"""Rules to determine the potential of renewables on a raster map and per administrative unit.

(1) Determine technical potential on a raster map.
(2) Define administrative units.
(3) Based on scenarios, restrict the potential, and allocate it to administrative units.

"""

root_dir = config["root-directory"] + "/" if config["root-directory"] not in ["", "."] else ""
script_dir = f"{root_dir}scripts/"

rule category_of_technical_eligibility:
    message:
        "Determine upper bound surface eligibility for renewables based on land cover, slope, bathymetry, and settlements."
    input:
        src = script_dir + "technical_eligibility.py",
        land_cover = rules.land_cover_in_europe.output[0],
        slope_pv = "build/slope-europe-pv.tif",
        slope_wind = "build/slope-europe-wind.tif",
        bathymetry = rules.bathymetry_in_europe.output[0],
        building_share = rules.settlements.output.buildings,
        urban_green_share = rules.settlements.output.urban_greens
    params:
        max_building_share = config["parameters"]["max-building-share"],
        max_urban_green_share = config["parameters"]["max-urban-green-share"],
        max_depth_offshore = config["parameters"]["max-depth-offshore"],
        slope_threshold = config["parameters"]["max-slope-pixel-fraction-threshold"]
    output:
        "build/technically-eligible-land.tif"
    conda: "../envs/default.yaml"
    script: "../scripts/technical_eligibility.py"


rule total_size_swiss_building_footprints_according_to_settlement_data:
    message: "Sum the size of building footprints from settlement data."
    input:
        src = script_dir + "swiss_building_footprints.py",
        building_footprints = rules.settlements.output.buildings,
        eligibility = "build/technically-eligible-land.tif",
        countries = rules.administrative_borders.output[0]
    output:
        "build/building-footprints-according-to-settlement-data-km2.txt"
    conda: "../envs/default.yaml"
    script: "../scripts/swiss_building_footprints.py"


rule correction_factor_building_footprint_to_available_rooftop:
    message: "Determine the factor that maps from building footprints to available rooftop area for CHE."
    input:
        rooftops = rules.total_size_swiss_rooftops_according_to_sonnendach_data.output[0],
        building_footprints = rules.total_size_swiss_building_footprints_according_to_settlement_data.output[0]
    output:
        "build/ratio-esm-available.txt"
    run:
        with open(input.rooftops, "r") as f_in:
            rooftops = float(f_in.read())
        with open(input.building_footprints, "r") as f_in:
            building_footprints = float(f_in.read())
        ratio = rooftops / building_footprints
        with open(output[0], "w") as f_out:
            f_out.write(f"{ratio:.3f}")


rule capacityfactor_of_technical_eligibility:
    message:
        "Determine capacityfactor of eligibility category."
    input:
        script = script_dir + "technically_eligible_capacityfactor.py",
        eligibility_categories = rules.category_of_technical_eligibility.output[0],
        rooftop_pv_cf = "build/capacityfactors/rooftop-pv-time-average.tif",
        open_field_pv_cf = "build/capacityfactors/open-field-pv-time-average.tif",
        wind_onshore_cf = "build/capacityfactors/wind-onshore-time-average.tif",
        wind_offshore_cf = "build/capacityfactors/wind-offshore-time-average.tif"
    params: availability = config["parameters"]["availability"]
    output:
        pv = "build/technically-eligible-capacityfactor-pv-prio.tif",
        wind = "build/technically-eligible-capacityfactor-wind-prio.tif"
    conda: "../envs/default.yaml"
    script: "../scripts/technically_eligible_capacityfactor.py"


rule area_of_technical_eligibility:
    message:
        "Quantify the area that is technically eligible for renewables."
    input:
        script = script_dir + "technically_eligible_area.py",
        eligibility_categories = rules.category_of_technical_eligibility.output[0],
        building_share = rules.settlements.output.buildings,
        rooftop_correction_factor = rules.correction_factor_building_footprint_to_available_rooftop.output[0]
    output:
        "build/technically-eligible-area-km2.tif"
    conda: "../envs/default.yaml"
    script: "../scripts/technically_eligible_area.py"


rule capacity_of_technical_eligibility:
    message:
        "Quantify the capacity that is technically eligible for renewables."
    input:
        script = script_dir + "technically_eligible_capacity.py",
        eligibility_categories = rules.category_of_technical_eligibility.output[0],
        eligible_areas = rules.area_of_technical_eligibility.output[0],
        statistical_roof_model = rules.sonnendach_statistics.output[0]
    params:
        maximum_installable_power_density = config["parameters"]["maximum-installable-power-density"]
    output:
        pv = "build/technically-eligible-capacity-pv-prio-mw.tif",
        wind = "build/technically-eligible-capacity-wind-prio-mw.tif",
    conda: "../envs/default.yaml"
    script: "../scripts/technically_eligible_capacity.py"


rule electricity_yield_of_technical_eligibility:
    message:
        "Quantify the max annual electricity yield that is technically eligible for renewables."
    input:
        script = script_dir + "technically_eligible_electricity_yield.py",
        eligibility_categories = rules.category_of_technical_eligibility.output[0],
        capacities_pv_prio = rules.capacity_of_technical_eligibility.output.pv,
        capacities_wind_prio = rules.capacity_of_technical_eligibility.output.wind,
        cf_pv_prio = rules.capacityfactor_of_technical_eligibility.output.pv,
        cf_wind_prio = rules.capacityfactor_of_technical_eligibility.output.wind
    output:
        pv = "build/technically-eligible-electricity-yield-pv-prio-twh.tif",
        wind = "build/technically-eligible-electricity-yield-wind-prio-twh.tif",
    conda: "../envs/default.yaml"
    script: "../scripts/technically_eligible_electricity_yield.py"


rule units:
    message: "Form units of layer {wildcards.layer} by remixing NUTS, LAU, and GADM."
    input:
        script = script_dir + "units.py",
        administrative_borders = rules.administrative_borders.output[0]
    params:
        layer_name = "{layer}",
        layer_config = lambda wildcards: config["layers"][wildcards.layer],
        countries = config["scope"]["countries"]
    output:
        "build/{layer}/units.geojson"
    conda: "../envs/default.yaml"
    script: "../scripts/units.py"


rule local_land_cover:
    message: "Land cover statistics per unit of layer {wildcards.layer}."
    input:
        script = script_dir + "land_cover_stats_to_csv.py",
        units = rules.units.output[0],
        land_cover = rules.land_cover_in_europe.output[0]
    params:
        attributes = [
            "lc_11", "lc_14", "lc_20", "lc_30", "lc_40", "lc_50",
            "lc_60", "lc_70", "lc_90", "lc_100", "lc_110", "lc_120",
            "lc_130", "lc_140", "lc_150", "lc_160", "lc_170", "lc_180",
            "lc_190", "lc_200", "lc_210", "lc_220", "lc_230"
        ]
    output:
        "build/{layer}/land-cover.csv"
    conda: "../envs/default.yaml"
    script: "../scripts/land_cover_stats_to_csv.py"


rule local_built_up_area:
    message: "Determine the built up area for administrative units in layer {wildcards.layer}."
    input:
        script = script_dir + "built_up_area.py",
        built_up_area = rules.settlements.output.built_up,
        units = rules.units.output[0]
    output:
        "build/{layer}/built-up-areas.csv"
    conda: "../envs/default.yaml"
    script: "../scripts/built_up_area.py"


rule shared_coast:
    message: "Determine share of coast length between eez and units of layer {wildcards.layer} using {threads} threads."
    input:
        script = script_dir + "shared_coast.py",
        units = rules.units.output[0],
        eez = rules.eez_in_europe.output[0]
    output:
        "build/{layer}/shared-coast.csv"
    threads: config["snakemake"]["max-threads"]
    conda: "../envs/default.yaml"
    script: "../scripts/shared_coast.py"


rule potentials:
    message:
        "Determine the constrained potentials for layer {wildcards.layer} in scenario {wildcards.scenario}."
    input:
        script = script_dir + "potentials.py",
        units = rules.units.output[0],
        eez = rules.eez_in_europe.output[0],
        shared_coast = rules.shared_coast.output[0],
        pv_yield = rules.electricity_yield_of_technical_eligibility.output.pv,
        wind_yield = rules.electricity_yield_of_technical_eligibility.output.wind,
        category = rules.category_of_technical_eligibility.output[0],
        land_cover = rules.land_cover_in_europe.output[0],
        protected_areas = rules.protected_areas_in_europe.output[0]
    params:
        scenario = lambda wildcards: config["scenarios"][wildcards.scenario],
        potential_metric = "electricity_yield"
    output:
        "build/{layer}/{scenario}/potentials.csv"
    conda: "../envs/default.yaml"
    script: "../scripts/potentials.py"


rule areas:
    message:
        "Determine eligible areas for layer {wildcards.layer} in scenario {wildcards.scenario}."
    input:
        script = script_dir + "areas.py",
        units = rules.units.output[0],
        eez = rules.eez_in_europe.output[0],
        shared_coast = rules.shared_coast.output[0],
        area = rules.area_of_technical_eligibility.output[0],
        category = rules.category_of_technical_eligibility.output[0],
        land_cover = rules.land_cover_in_europe.output[0],
        protected_areas = rules.protected_areas_in_europe.output[0]
    params:
        scenario = lambda wildcards: config["scenarios"][wildcards.scenario]
    output:
        "build/{layer}/{scenario}/areas.csv"
    conda: "../envs/default.yaml"
    script: "../scripts/areas.py"


rule capacities:
    message:
        "Determine installable capacities for layer {wildcards.layer} in scenario {wildcards.scenario}."
    input:
        script = script_dir + "potentials.py",
        units = rules.units.output[0],
        eez = rules.eez_in_europe.output[0],
        shared_coast = rules.shared_coast.output[0],
        pv_capacity = rules.capacity_of_technical_eligibility.output.pv,
        wind_capacity = rules.capacity_of_technical_eligibility.output.wind,
        pv_yield = rules.electricity_yield_of_technical_eligibility.output.pv,
        wind_yield = rules.electricity_yield_of_technical_eligibility.output.wind,
        category = rules.category_of_technical_eligibility.output[0],
        land_cover = rules.land_cover_in_europe.output[0],
        protected_areas = rules.protected_areas_in_europe.output[0]
    params:
        scenario = lambda wildcards: config["scenarios"][wildcards.scenario],
        potential_metric = "electricity_yield"
    output:
        "build/{layer}/{scenario}/capacities.csv"
    conda: "../envs/default.yaml"
    script: "../scripts/potentials.py"
