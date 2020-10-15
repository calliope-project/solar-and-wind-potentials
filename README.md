# Possibility for renewable electricity autarky in Europe

Is your European region able to provide itself with 100% renewable electricity?

This repository contains the entire research project, including code and report. The philosophy behind this repository is that no intermediary results are included, but all results are computed from raw data and code.

[![article DOI](https://img.shields.io/badge/article-10.1016%2Fj.esr.2019.100388-blue)](https://doi.org/10.1016/j.esr.2019.100388)
[![data DOI](https://img.shields.io/badge/data-10.5281%2Fzenodo.3244985-blue)](https://doi.org/10.5281/zenodo.3244985)
[![code DOI](https://img.shields.io/badge/code-10.5281%2Fzenodo.3246303-blue)](https://doi.org/10.5281/zenodo.3246303)

## Getting ready

### Installation

The following dependencies are needed to set up an environment in which the analysis can be run and the paper be build:

* [conda](https://conda.io/docs/index.html)

When these dependencies are installed, you can create a conda environment from within you can run the analysis:

    conda env create -f environment.yaml

Don't forget to activate the environment. To see what you can do now, run:

    snakemake --list

### Data to be retrieved manually

Whenever possible, data is downloaded automatically. As this is not always possible, you will need to retrieve the following data sets manually:

* [European Settlement Map 2012, Release 2017, 100m](https://land.copernicus.eu/pan-european/GHSL/european-settlement-map), to be placed at `./data/esm-100m-2017/`
* [World Exclusive Economic Zones v10](http://www.marineregions.org/downloads.php), to be placed in `./data/World_EEZ_v10_20180221`
* capacity factors from renewable.ninja, to be placed in `./data/capacityfactors/{technology}.nc` for technology in ["wind-onshore", "wind-offshore", "rooftop-pv", "open-field-pv"] (where "open-field-pv" and "rooftop-pv" can be the same dataset and hence can be linked instead of copied)(to run simulations, see `Manual steps` below)

## Run the analysis

    snakemake --use-conda paper

This will run all analysis steps to reproduce results and eventually build the paper.

You can also run certain parts only by using other `snakemake` rules; to get a list of all rules run `snakemake --list`.

To generate a PDF of the dependency graph of all steps, run:

    snakemake --rulegraph | dot -Tpdf > dag.pdf

(needs `dot`: `conda install graphviz`).

## Run on Euler cluster

To run on Euler, use the following command:

    snakemake --use-conda --profile config/euler

If you want to run on another cluster, read [snakemake's documentation on cluster execution](https://snakemake.readthedocs.io/en/stable/executable.html#cluster-execution) and take `config/euler` as a starting point.

## Manual steps

At the moment, there is one manual step involved: running renewables.ninja simulations of wind and solar electricity. It is added to the automatic workflow as input data. Should you want to change the simulations, because you want to change parameters of the simulation (see `parameters.ninja` in the config), you can do that in three steps:

1) Create input files by first changing the config, then running `snakemake -s rules/ninja-input.smk`.
2) Run the simulations on renewables.ninja.
3) Update the data in `data/capacityfactors/{technology}`.

## Run the tests

    snakemake --use-conda test

## Repo structure

* `report`: contains all files necessary to build the paper; plots and result files are not in here but generated automatically
* `src`: contains the Python source code
* `tests`: contains the test code
* `config`: configurations used in the study
* `rules`: additional Snakemake rules and workflows
* `data`: place for raw data, whether retrieved manually and automatically
* `build`: will contain all results (does not exist initially)

## Citation

If you use this code or data in an academic publication, please see `CITATION.md`.

## License

The code in this repo is MIT licensed, see `./LICENSE.md`. This excludes the KlinicSlab font family (all files in `./report/fonts/`) which is copyright Lost Type.
