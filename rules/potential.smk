"""Rules to determine the potential of renewables on a raster map and per administrative unit.

(1) Determine technical potential on a raster map.
(2) Define administrative units.
(3) Based on scenarios, restrict the potential, and allocate it to administrative units.

"""

PYTHON = "PYTHONPATH=./ python"

root_dir = config["root-directory"] + "/" if config["root-directory"] not in ["", "."] else ""
script_dir = f"{root_dir}scripts/"

rule category_of_technical_eligibility:
    message:
        "Determine upper bound surface eligibility for renewables based on land cover, slope, bathymetry, and settlements."
    input:
        src = script_dir + "technical_eligibility.py",
        land_cover = rules.land_cover_in_europe.output,
        slope = rules.slope_in_europe.output,
        bathymetry = rules.bathymetry_in_europe.output,
        building_share = rules.settlements.output.buildings,
        urban_green_share = rules.settlements.output.urban_greens
    params:
        max_slope = config["parameters"]["max-slope"],
        max_building_share = config["parameters"]["max-building-share"],
        max_urban_green_share = config["parameters"]["max-urban-green-share"]
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
        countries = rules.administrative_borders.output
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
        eligibility_categories = rules.category_of_technical_eligibility.output,
        capacity_factors = expand(
            "build/capacityfactors/{technology}-time-average.tif",
            technology=["rooftop-pv", "open-field-pv", "wind-onshore", "wind-offshore"]
        )
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
        eligibility_categories = rules.category_of_technical_eligibility.output,
        building_share = rules.settlements.output.buildings,
        rooftop_correction_factor = rules.correction_factor_building_footprint_to_available_rooftop.output
    output:
        "build/technically-eligible-area-km2.tif"
    conda: "../envs/default.yaml"
    script: "../scripts/technically_eligible_area.py"


rule capacity_of_technical_eligibility:
    message:
        "Quantify the capacity that is technically eligible for renewables."
    input:
        script = script_dir + "technically_eligible_capacity.py",
        ligibility_categories = rules.category_of_technical_eligibility.output,
        eligible_areas = rules.area_of_technical_eligibility.output,
        statistical_roof_model = rules.sonnendach_statistics.output
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
        eligibility_categories = rules.category_of_technical_eligibility.output,
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
        administrative_borders = rules.administrative_borders.output
    params:
        layer_name = "{layer}",
        layer_config = lambda wildcards: config["shapes"][wildcards.layer],
        countries = config["scope"]["countries"]
    output:
        "build/{layer}/units.geojson"
    conda: "../envs/default.yaml"
    script: "../scripts/units.py"


rule local_land_cover:
    message: "Land cover statistics per unit of layer {wildcards.layer}."
    input:
        script = script_dir + "geojson_to_csv.py",
        units = rules.units.output,
        land_cover = rules.land_cover_in_europe.output
    output:
        "build/{layer}/land-cover.csv"
    conda: "../envs/default.yaml"
    shell:
        """
        fio cat {input.units} | \
        rio zonalstats -r {input.land_cover} --prefix 'lc_' --categorical | \
        {PYTHON} {input.script} -a id -a lc_11 -a lc_14 -a lc_20 -a lc_30 -a lc_40 \
        -a lc_50 -a lc_60 -a lc_70 -a lc_90 -a lc_100 -a lc_110 -a lc_120 -a lc_130 \
        -a lc_140 -a lc_150 -a lc_160 -a lc_170 -a lc_180 -a lc_190 -a lc_200 \
        -a lc_210 -a lc_220 -a lc_230 > \
        {output}
        """


rule local_built_up_area:
    message: "Determine the built up area for administrative units in layer {wildcards.layer}."
    input:
        script = script_dir + "built_up_area.py",
        built_up_area = rules.settlements.output.built_up,
        units = rules.units.output
    output:
        "build/{layer}/built-up-areas.csv"
    conda: "../envs/default.yaml"
    script: "../scripts/built_up_area.py"


# FIXME: where is eligibility_local.py??
rule eez_eligibility:
    message:
        "Allocate eligible land to exclusive economic zones using {threads} threads."
    input:
        src = script_dir + "eligibility_local.py",
        regions = rules.eez_in_europe.output,
        eligibility = rules.category_of_technical_eligibility.output
    output:
        "build/eez-eligibility.csv"
    threads: config["snakemake"]["max-threads"]
    conda: "../envs/default.yaml"
    shell:
        PYTHON + " {input.src} offshore {input.regions} {input.eligibility} {output} {threads}"


rule shared_coast:
    message: "Determine share of coast length between eez and units of layer {wildcards.layer} using {threads} threads."
    input:
        script = script_dir + "shared_coast.py",
        units = rules.units.output,
        eez = rules.eez_in_europe.output
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
        units = rules.units.output,
        eez = rules.eez_in_europe.output,
        shared_coast = rules.shared_coast.output,
        pv_yield = rules.electricity_yield_of_technical_eligibility.output.pv,
        wind_yield = rules.electricity_yield_of_technical_eligibility.output.wind,
        category = rules.category_of_technical_eligibility.output,
        land_cover = rules.land_cover_in_europe.output,
        protected_areas = rules.protected_areas_in_europe.output
    params:
        scenario = lambda wildcards: config["scenarios"][wildcards.scenario]
    output:
        "build/{layer}/{scenario}/potentials.csv"
    conda: "../envs/default.yaml"
    script: "../scripts/potentials.py"


rule areas:
    message:
        "Determine eligible areas for layer {wildcards.layer} in scenario {wildcards.scenario}."
    input:
        script = script_dir + "areas.py",
        units = rules.units.output,
        eez = rules.eez_in_europe.output,
        shared_coast = rules.shared_coast.output,
        area = rules.area_of_technical_eligibility.output,
        category = rules.category_of_technical_eligibility.output,
        land_cover = rules.land_cover_in_europe.output,
        protected_areas = rules.protected_areas_in_europe.output
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
        script = script_dir + "capacities.py",
        units = rules.units.output,
        eez = rules.eez_in_europe.output,
        shared_coast = rules.shared_coast.output,
        capacity_pv = rules.capacity_of_technical_eligibility.output.pv,
        capacity_wind = rules.capacity_of_technical_eligibility.output.wind,
        pv_yield = rules.electricity_yield_of_technical_eligibility.output.pv,
        wind_yield = rules.electricity_yield_of_technical_eligibility.output.wind,
        category = rules.category_of_technical_eligibility.output,
        land_cover = rules.land_cover_in_europe.output,
        protected_areas = rules.protected_areas_in_europe.output
    params:
        scenario = lambda wildcards: config["scenarios"][wildcards.scenario]
    output:
        "build/{layer}/{scenario}/capacities.csv"
    conda: "../envs/default.yaml"
    script: "../scripts/capacities.py"
