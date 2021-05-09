"""Workflow to create simulation input for renewables.ninja.

We create a raster grid on top of a map of Europe in order of running one (wind)
or several (different roof configurations for pv) simulations per raster point.
"""

configfile: "config/default.yaml"

include: "../Snakefile"
include: "sonnendach.smk"
root_dir = config["root-directory"] + "/" if config["root-directory"] not in ["", "."] else ""
script_dir = f"{root_dir}scripts/"


rule ninja_simulation_input:
    message: "Create input files for renewable.ninja simulations."
    input:
        "build/capacityfactors/ninja-input-pv.csv",
        "build/capacityfactors/ninja-input-wind-onshore.csv",
        "build/capacityfactors/ninja-input-wind-offshore.csv"


rule pv_simulation_points:
    message: "Create locations and parameters of pv simulations for renewables.ninja."
    input:
        script = script_dir + "capacityfactors/ninja_input_pv.py",
        units = "build/continental/units.geojson",
        roof_categories = rules.sonnendach_statistics.output[0]
    params:
        bounds = config["scope"]["bounds"],
        ninja = config["parameters"]["ninja"],
        maximum_power_density = config["parameters"]["maximum-installable-power-density"]
    output:
        points = "build/capacityfactors/ninja-input-pv.csv"
    conda: "../envs/default.yaml"
    script: "../scripts/capacityfactors/ninja_input_pv.py"


rule wind_simulation_points:
    message: "Create locations and parameters of wind simulations for renewables.ninja."
    input:
        script = script_dir + "capacityfactors/ninja_input_wind.py",
        units = "build/continental/units.geojson",
        eez = "build/eez-in-europe.geojson"
    params:
        bounds = config["scope"]["bounds"],
        ninja = config["parameters"]["ninja"],
    output:
        points_onshore = "build/capacityfactors/ninja-input-wind-onshore.csv",
        points_offhore = "build/capacityfactors/ninja-input-wind-offshore.csv",
    conda: "../envs/default.yaml"
    script: "../scripts/capacityfactors/ninja_input_wind.py"
