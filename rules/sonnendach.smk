"""Retrieve statistics of the sonnendach.ch dataset."""

localrules: total_size_swiss_rooftops_according_to_sonnendach_data,
    total_swiss_yield_according_to_sonnendach_data, sonnendach_statistics


rule total_size_swiss_rooftops_according_to_sonnendach_data:
    message: "Download the size of rooftops from Sonnendach data."
    params: url = config["data-sources"]["sonnendach_statistics"]
    output: "data/automatic/sonnendach/total-rooftop-area-km2.txt"
    shell:
        "curl -sLo {output} '{params.url}'"


rule total_swiss_yield_according_to_sonnendach_data:
    message: "Download the yield of all available rooftops from Sonnendach data."
    params: url = config["data-sources"]["sonnendach_total_size"]
    output: "data/automatic/sonnendach/total-yield-twh.txt"
    shell:
        "curl -sLo {output} '{params.url}'"


rule sonnendach_statistics:
    message: "Download statistics of roofs in Switzerland."
    params: url = config["data-sources"]["sonnendach_total_yield"]
    output: "data/automatic/sonnendach/roof-statistics.csv",
    shell:
        "curl -sLo {output} '{params.url}'"
