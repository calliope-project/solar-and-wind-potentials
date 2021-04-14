from snakemake.utils import validate

configfile: "config/default.yaml"
validate(config, "config/schema.yaml")

include: "rules/data-preprocessing.smk"
include: "rules/sonnendach.smk"
include: "rules/capacityfactors.smk"
include: "rules/potential.smk"
include: "rules/sync.smk"

localrules: all, clean
root_dir = config["root-directory"] + "/" if config["root-directory"] not in ["", "."] else ""
__version__ = open(f"{root_dir}VERSION").readlines()[0].strip()
script_dir = f"{root_dir}scripts/"

wildcard_constraints:
    layer = "({layer_list})".format(layer_list="|".join((f"({layer})" for layer in config["layers"]))),
    scenario = "({scenario_list})".format(scenario_list="|".join((f"({scenario})" for scenario in config["scenarios"])))

onstart:
    shell("mkdir -p build/logs")
onsuccess:
    if "pushcut_secret" in config.keys():
        trigger_pushcut(event_name="snakemake_succeeded", secret=config["pushcut_secret"])
onerror:
    if "pushcut_secret" in config.keys():
        trigger_pushcut(event_name="snakemake_failed", secret=config["pushcut_secret"])


rule all:
    message: "Run entire analysis and compile report."
    input:
        "build/logs/test-report.html",
        expand(
            "build/{layer}/{scenario}/potentials.csv",
            scenario=config["scenarios"],
            layer=config["layers"]
        ),
        expand(
            "build/{layer}/{scenario}/capacities.csv",
            scenario=config["scenarios"],
            layer=config["layers"]
        ),
        expand(
            "build/{layer}/{scenario}/areas.csv",
            scenario=config["scenarios"],
            layer=config["layers"]
        )


rule clean: # removes all generated results
    shell:
        """
        rm -r ./build/*
        echo "Data downloaded to data/ has not been cleaned."
        """


rule test:
    message: "Run tests."
    input:
        #expand("build/{layer}/technical-potential/potentials.csv", layer=config["layers"]),
        #"build/technically-eligible-land.tif",
        #"build/technically-eligible-area-km2.tif",
        #"build/technically-eligible-electricity-yield-pv-prio-twh.tif",
        "build/administrative-borders.gpkg",
        #"data/automatic/sonnendach/total-rooftop-area-km2.txt",
        #"data/automatic/sonnendach/total-yield-twh.txt"
    output: "build/logs/test-report.html"
    conda: "envs/default.yaml"
    shell:
        "py.test --html={output} --self-contained-html"


def trigger_pushcut(event_name, secret):
    import requests
    response = requests.post(
            f'https://api.pushcut.io/{secret}/notifications/{event_name}'
    )
