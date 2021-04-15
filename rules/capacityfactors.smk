"""Workflow to create spatially disaggregated capacity factors from renewables.ninja simulations.

This is based on simulations run based on input generated in the workflow `ninja-input.smk`. Beware
that renewables.ninja simulations are not in the loop, i.e. they are not run automatically but must
be run manually if they need to be altered.
"""

configfile: "./config/default.yaml"
root_dir = config["root-directory"] + "/" if config["root-directory"] not in ["", "."] else ""
script_dir = f"{root_dir}scripts/"


rule capacityfactor_timeseries:
    message: "Create index capacity factor timeseries of {wildcards.technology}."
    input:
        script = script_dir + "capacityfactors/timeseries.py",
        capacityfactor = lambda wildcards: config["data-sources"]["raw-capacity-factors"].format(technology=wildcards.technology)
    output:
        "build/capacityfactors/{technology}-timeseries.nc"
    conda: "../envs/default.yaml"
    script: "../scripts/capacityfactors/timeseries.py"


rule capacityfactor_id_map:
    message: "Create raster map of indices to time series for {wildcards.technology}."
    input:
        script = script_dir + "capacityfactors/id_map.py",
        timeseries = "build/capacityfactors/{technology}-timeseries.nc"
    output:
        temp("build/{technology}-ids-lowres.tif")
    shadow: "full"
    params:
        resolution = config["parameters"]["ninja"]["resolution-grid"]
    conda: "../envs/default.yaml"
    script: "../scripts/capacityfactors/id_map.py"


rule capacityfactor_id_map_warped_to_land_cover:
    message: "Warp raster map of indices for {wildcards.technology} to land cover map resolution."
    input:
        id_map = "build/{technology}-ids-lowres.tif",
        reference = "build/land-cover-europe.tif"
    output:
        "build/capacityfactors/{technology}-ids.tif"
    shadow: "full"
    conda: "../envs/default.yaml"
    shell: "rio warp {input.id_map} {output} --like {input.reference} --resampling nearest"


rule time_average_capacityfactor_map:
    message: "Create raster map of average capacity factors for {wildcards.technology}."
    input:
        scripts = script_dir + "capacityfactors/averages_map.py",
        id_map = "build/capacityfactors/{technology}-ids.tif",
        timeseries = rules.capacityfactor_timeseries.output
    output:
        "build/capacityfactors/{technology}-time-average.tif"
    conda: "../envs/default.yaml"
    script: "../scripts/capacityfactors/averages_map.py"
