"""Determine potential of renewable electricity in each administrative unit.

* Take the (only technically restricted) raster data potentials,
* add restrictions based on scenario definitions,
* allocate the onshore potentials to the administrative units,
* allocate the offshore potentials to exclusive economic zones (EEZ),
* allocate the offshore potential of EEZ to units based on the fraction of shared coast.

This is in analogy to `areas.py` but for potentials [TWh/a] rather than areas [km2] .
"""

import numpy as np
import pandas as pd
import rasterio
from rasterstats import zonal_stats
import fiona

from renewablepotentialslib.eligibility import Eligibility, FOREST, FARM, OTHER


def potentials(path_to_units, path_to_eez, path_to_shared_coast,
               path_to_capacities_pv_prio, path_to_capacities_wind_prio,
               path_to_electricity_yield_pv_prio, path_to_electricity_yield_wind_prio,
               path_to_eligibility_categories, path_to_land_cover, path_to_protected_areas,
               scenario_config, metric, path_to_result):
    """Determine potential of renewable electricity in each administrative unit.

    * Take the (only technically restricted) raster data potentials,
    * add restrictions based on scenario definitions,
    * allocate the onshore potentials to the administrative units,
    * allocate the offshore potentials to exclusive economic zones (EEZ),
    * allocate the offshore potential of EEZ to units based on the fraction of shared coast.
    """
    with rasterio.open(path_to_eligibility_categories, "r") as src:
        eligibility_categories = src.read(1)
    with rasterio.open(path_to_electricity_yield_pv_prio, "r") as src:
        transform = src.transform
        electricity_yield_pv_prio = src.read(1)
    with rasterio.open(path_to_electricity_yield_wind_prio, "r") as src:
        electricity_yield_wind_prio = src.read(1)
    with rasterio.open(path_to_land_cover, "r") as src:
        land_cover = src.read(1)
    with rasterio.open(path_to_protected_areas, "r") as src:
        protected_areas = src.read(1)
    with fiona.open(path_to_units, "r") as src:
        unit_ids = [feature["properties"]["id"] for feature in src]
        unit_geometries = [feature["geometry"] for feature in src]
    with fiona.open(path_to_eez, "r") as src:
        eez_ids = [feature["properties"]["id"] for feature in src]
        eez_geometries = [feature["geometry"] for feature in src]
    shared_coasts = pd.read_csv(path_to_shared_coast, index_col=0)
    electricity_yield_pv_prio, electricity_yield_wind_prio = apply_scenario_config(
        potential_pv_prio=electricity_yield_pv_prio,
        potential_wind_prio=electricity_yield_wind_prio,
        categories=eligibility_categories,
        land_cover=land_cover,
        protected_areas=protected_areas,
        scenario_config=scenario_config
    )
    if metric == "capacity":
        with rasterio.open(path_to_capacities_pv_prio, "r") as src:
            transform = src.transform
            capacities_pv_prio = src.read(1)
        with rasterio.open(path_to_capacities_wind_prio, "r") as src:
            capacities_wind_prio = src.read(1)

        potential_pv_prio, potential_wind_prio = apply_scenario_config(
            potential_pv_prio=capacities_pv_prio,
            potential_wind_prio=capacities_wind_prio,
            categories=eligibility_categories,
            land_cover=land_cover,
            protected_areas=protected_areas,
            scenario_config=scenario_config
        )
    elif metric == "electricity_yield":
        potential_pv_prio = electricity_yield_pv_prio
        potential_wind_prio = electricity_yield_wind_prio

    capacities_pv_prio, capacities_wind_prio = decide_between_pv_and_wind(
        potential_pv_prio=potential_pv_prio,
        potential_wind_prio=potential_wind_prio,
        electricity_yield_pv_prio=electricity_yield_pv_prio,
        electricity_yield_wind_prio=electricity_yield_wind_prio,
        eligibility_categories=eligibility_categories
    )

    onshore_potentials = pd.DataFrame(
        index=unit_ids,
        data={
            getattr(potential, f"{metric}_name"): potentials_per_shape(
                eligibilities=potential.eligible_on,
                potential_map=(capacities_pv_prio if "pv" in str(potential).lower()
                               else capacities_wind_prio),
                eligibility_categories=eligibility_categories,
                shapes=unit_geometries,
                transform=transform
            )
            for potential in Potential.onshore()
        }
    )
    offshore_eez_potentials = pd.DataFrame(
        index=eez_ids,
        data={
            getattr(potential, f"{metric}_name"): potentials_per_shape(
                eligibilities=potential.eligible_on,
                potential_map=(capacities_pv_prio if "pv" in str(potential).lower()
                               else capacities_wind_prio),
                eligibility_categories=eligibility_categories,
                shapes=eez_geometries,
                transform=transform
            )
            for potential in Potential.offshore()
        }
    )
    offshore_potentials = pd.DataFrame(
        data=shared_coasts.dot(offshore_eez_potentials),
        columns=[getattr(potential, f"{metric}_name") for potential in Potential.offshore()]
    )
    potentials = pd.concat([onshore_potentials, offshore_potentials], axis=1)
    potentials.index.name = "id"
    potentials.to_csv(
        path_to_result,
        header=True,
        index=True
    )


def apply_scenario_config(potential_pv_prio, potential_wind_prio, categories,
                          land_cover, protected_areas, scenario_config):
    """Limit potential in each pixel based on scenario config."""

    # share-rooftops-used
    share_rooftops_used = scenario_config["share-rooftops-used"]
    mask = categories == Eligibility.ROOFTOP_PV
    potential_pv_prio[mask] = potential_pv_prio[mask] * share_rooftops_used
    potential_wind_prio[mask] = potential_wind_prio[mask] * share_rooftops_used

    # share-forest-used-for-wind
    share_forest_used_for_wind = scenario_config["share-forest-used-for-wind"]
    mask = np.isin(land_cover, FOREST) & (categories != Eligibility.ROOFTOP_PV)
    potential_pv_prio[mask] = potential_pv_prio[mask] * share_forest_used_for_wind
    potential_wind_prio[mask] = potential_wind_prio[mask] * share_forest_used_for_wind

    # share-other-land-used
    share_other_land_used = scenario_config["share-other-land-used"]
    mask = np.isin(land_cover, OTHER) & (categories != Eligibility.ROOFTOP_PV)
    potential_pv_prio[mask] = potential_pv_prio[mask] * share_other_land_used
    potential_wind_prio[mask] = potential_wind_prio[mask] * share_other_land_used

    # share-farmland-used
    share_farmland_used = scenario_config["share-farmland-used"]
    mask = np.isin(land_cover, FARM) & (categories != Eligibility.ROOFTOP_PV)
    potential_pv_prio[mask] = potential_pv_prio[mask] * share_farmland_used
    potential_wind_prio[mask] = potential_wind_prio[mask] * share_farmland_used

    # share-offshore-used
    share_offshore_used = scenario_config["share-offshore-used"]
    mask = categories == Eligibility.OFFSHORE_WIND
    potential_pv_prio[mask] = potential_pv_prio[mask] * share_offshore_used
    potential_wind_prio[mask] = potential_wind_prio[mask] * share_offshore_used

    # pv-on-farmland
    pv_on_farmland = scenario_config["pv-on-farmland"]
    if not pv_on_farmland:
        mask = np.isin(land_cover, FARM) & (categories == Eligibility.ONSHORE_WIND_AND_PV)
        potential_pv_prio[mask] = 0

    # share-protected-areas-used
    use_protected_areas = scenario_config["use-protected-areas"]
    if not use_protected_areas:
        mask = (protected_areas == ProtectedArea.PROTECTED) & (categories != Eligibility.ROOFTOP_PV)
        potential_pv_prio[mask] = 0
        potential_wind_prio[mask] = 0

    return potential_pv_prio, potential_wind_prio


def decide_between_pv_and_wind(potential_pv_prio, potential_wind_prio,
                               electricity_yield_pv_prio, electricity_yield_wind_prio,
                               eligibility_categories):
    """When both are possible, choose PV when its electricity yield is higher, or vice versa."""
    pv_and_wind_possible = eligibility_categories == Eligibility.ONSHORE_WIND_AND_PV
    higher_wind_yield = electricity_yield_pv_prio <= electricity_yield_wind_prio

    potential_pv_prio[pv_and_wind_possible & higher_wind_yield] = 0
    potential_wind_prio[pv_and_wind_possible & ~higher_wind_yield] = 0

    return potential_pv_prio, potential_wind_prio


def potentials_per_shape(eligibilities, potential_map, eligibility_categories, shapes, transform):
    """Determine potential of one eligibility category per shape."""
    potential_map = potential_map.copy()
    potential_map[~np.isin(eligibility_categories, eligibilities)] = 0
    potentials = zonal_stats(
        shapes,
        potential_map,
        affine=transform,
        stats="sum",
        nodata=-999
    )
    return [stat["sum"] for stat in potentials]


if __name__ == "__main__":
    potentials(
        path_to_units=snakemake.input.units,
        path_to_eez=snakemake.input.eez,
        path_to_shared_coast=snakemake.input.shared_coast,
        path_to_capacities_pv_prio=snakemake.input.get("pv_capacity", None),
        path_to_capacities_wind_prio=snakemake.input.get("wind_capacity", None),
        path_to_electricity_yield_pv_prio=snakemake.input.pv_yield,
        path_to_electricity_yield_wind_prio=snakemake.input.wind_yield,
        path_to_eligibility_categories=snakemake.input.category,
        path_to_land_cover=snakemake.input.land_cover,
        path_to_protected_areas=snakemake.input.protected_areas,
        scenario_config=snakemake.params.scenario,
        metric=snakemake.params.potential_metric,
        path_to_result=snakemake.output[0]
    )
