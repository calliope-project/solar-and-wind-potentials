"""Module to postprocess population data."""
import warnings

import click
import pandas as pd
import geopandas as gpd

from src.technical_eligibility import GlobCover
from src.conversion import area_in_squaremeters

WATER_THRESHOLD = 0.9 # units above this threshold are considered pure water bodies
POPULATION_THRESHOLD = 0.001 # share of population that can be removed

WATER = f"lc_{GlobCover.WATER_BODIES.value}"
NOT_WATER = [f"lc_{x.value}" for x in GlobCover if x is not GlobCover.WATER_BODIES]


@click.command()
@click.argument("path_to_units")
@click.argument("path_to_land_cover_data")
def postprocess_population(path_to_units, path_to_land_cover_data):
    """Postprocesses population data.

    (1) Removes nans, which can occur when units are too small to cover a population raster.
    (2) Removes population living in water bodies.
    (3) Calculates population density.

    Implemented as a filter which reads from stdin and writes to stdout.
    """
    population = pd.read_csv(click.get_text_stream('stdin'), index_col="id")
    population["population_sum"] = _fillna(population)
    population["population_sum"] = _remove_water_bodies(population, pd.read_csv(path_to_land_cover_data))
    population["density_p_per_km2"] = _calculate_density(population, gpd.read_file(path_to_units).set_index("id"))
    population[["population_sum", "density_p_per_km2"]].to_csv(click.get_text_stream('stdout'), header=True)


def _fillna(population):
    population["population_sum"] = pd.to_numeric(population["population_sum"], errors="coerce")
    return population["population_sum"].fillna(0.0)


def _remove_water_bodies(population, land_cover):
    data = land_cover.merge(population.reset_index(), on="id").fillna(0.0)
    data["water"] = data[WATER]
    data["not_water"] = data[NOT_WATER].sum(axis=1)
    data["rel_water"] = data["water"] / (data["not_water"] + data["water"])
    invalid_mask = (data["proper"] is False) & (data["rel_water"] > WATER_THRESHOLD)
    invalid_population = data.loc[invalid_mask, "population_sum"].sum()
    if invalid_population > 0:
        warnings.warn(f"Removing {invalid_population:.0f} people because they are in water bodies.")
    assert invalid_population < POPULATION_THRESHOLD * data["population_sum"].sum()
    data.loc[invalid_mask, "population_sum"] = 0.0
    return data.set_index("id")["population_sum"]


def _calculate_density(population, units):
    area_in_km2 = area_in_squaremeters(units) / 1e6
    return population["population_sum"] / area_in_km2


if __name__ == "__main__":
    postprocess_population()
