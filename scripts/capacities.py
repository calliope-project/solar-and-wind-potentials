"""Determine potential of renewable electricity in each administrative unit.

* Take the (only technically restricted) raster data potentials,
* add restrictions based on scenario definitions,
* allocate the onshore potentials to the administrative units,
* allocate the offshore potentials to exclusive economic zones (EEZ),
* allocate the offshore potential of EEZ to units based on the fraction of shared coast.

This is in analogy to `potentials.py` and `areas.py` but for installable capacities [MW]
rather than areas [km2] or [TWh/a].
"""
import pandas as pd
import rasterio
import fiona

from potentials import Potential, apply_scenario_config, decide_between_pv_and_wind, potentials_per_shape


def potentials(path_to_units, path_to_eez, path_to_shared_coast,
               path_to_capacities_pv_prio, path_to_capacities_wind_prio,
               path_to_electricity_yield_pv_prio, path_to_electricity_yield_wind_prio,
               path_to_eligibility_categories, path_to_land_cover, path_to_protected_areas,
               scenario_config, path_to_result):
    """Determine potential of renewable electricity in each administrative unit.

    * Take the (only technically restricted) raster data potentials,
    * add restrictions based on scenario definitions,
    * allocate the onshore potentials to the administrative units,
    * allocate the offshore potentials to exclusive economic zones (EEZ),
    * allocate the offshore potential of EEZ to units based on the fraction of shared coast.
    """
    with rasterio.open(path_to_eligibility_categories, "r") as src:
        eligibility_categories = src.read(1)
    with rasterio.open(path_to_capacities_pv_prio, "r") as src:
        transform = src.transform
        capacities_pv_prio = src.read(1)
    with rasterio.open(path_to_capacities_wind_prio, "r") as src:
        capacities_wind_prio = src.read(1)
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

    capacities_pv_prio, capacities_wind_prio = apply_scenario_config(
        potential_pv_prio=capacities_pv_prio,
        potential_wind_prio=capacities_wind_prio,
        categories=eligibility_categories,
        land_cover=land_cover,
        protected_areas=protected_areas,
        scenario_config=scenario_config
    )
    electricity_yield_pv_prio, electricity_yield_wind_prio = apply_scenario_config(
        potential_pv_prio=electricity_yield_pv_prio,
        potential_wind_prio=electricity_yield_wind_prio,
        categories=eligibility_categories,
        land_cover=land_cover,
        protected_areas=protected_areas,
        scenario_config=scenario_config
    )
    capacities_pv_prio, capacities_wind_prio = decide_between_pv_and_wind(
        potential_pv_prio=capacities_pv_prio,
        potential_wind_prio=capacities_wind_prio,
        electricity_yield_pv_prio=electricity_yield_pv_prio,
        electricity_yield_wind_prio=electricity_yield_wind_prio,
        eligibility_categories=eligibility_categories
    )

    onshore_potentials = pd.DataFrame(
        index=unit_ids,
        data={
            potential.capacity_name: potentials_per_shape(
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
            potential.capacity_name: potentials_per_shape(
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
        columns=[potential.capacity_name for potential in Potential.offshore()]
    )
    potentials = pd.concat([onshore_potentials, offshore_potentials], axis=1)
    potentials.index.name = "id"
    potentials.to_csv(
        path_to_result,
        header=True,
        index=True
    )


if __name__ == "__main__":
    potentials(
        path_to_units=snakemake.input.units,
        path_to_eez=snakemake.input.eez,
        path_to_shared_coast=snakemake.input.shared_coast,
        path_to_capacities_pv_prio=snakemake.input.capacity_pv,
        path_to_capacities_wind_prio=snakemake.input.capacity_wind,
        path_to_electricity_yield_pv_prio=snakemake.input.pv_yield,
        path_to_electricity_yield_wind_prio=snakemake.input.wind_yield,
        path_to_eligibility_categories=snakemake.input.category,
        path_to_land_cover=snakemake.input.land_cover,
        path_to_protected_areas=snakemake.input.protected_areas,
        scenario_config=snakemake.params.scenario,
        path_to_result=snakemake.output[0]
    )
