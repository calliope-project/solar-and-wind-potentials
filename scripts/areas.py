"""Determine available area of renewable electricity in each administrative unit.

* Take the (only technically restricted) raster data potentials,
* add restrictions based on scenario definitions,
* allocate the onshore areas to the administrative units,
* allocate the offshore areas to exclusive economic zones (EEZ),
* allocate the offshore areas of EEZ to units based on the fraction of shared coast.

This is in analogy to `potentials.py` but for areas [km2] instead of potentials [TWh/a].
"""
import numpy as np
import pandas as pd
import rasterio
from rasterstats import zonal_stats
import fiona

from src.technical_eligibility import Eligibility, FOREST, FARM, OTHER
from src.potentials import ProtectedArea


def areas(path_to_units, path_to_eez, path_to_shared_coast, path_to_eligible_area,
          path_to_eligibility_categories, path_to_land_cover, path_to_protected_areas,
          scenario_config, path_to_result):
    """Determine available area of renewable electricity in each administrative unit.

    * Take the (only technically restricted) raster data potentials,
    * add restrictions based on scenario definitions,
    * allocate the onshore areas to the administrative units,
    * allocate the offshore areas to exclusive economic zones (EEZ),
    * allocate the offshore areas of EEZ to units based on the fraction of shared coast.
    """
    with rasterio.open(path_to_eligibility_categories, "r") as src:
        category_map = src.read(1)
    with rasterio.open(path_to_eligible_area, "r") as src:
        transform = src.transform
        area_map = src.read(1)
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

    area_map = apply_scenario_config_to_areas(
        area_map=area_map,
        category_map=category_map,
        land_cover=land_cover,
        protected_areas=protected_areas,
        scenario_config=scenario_config
    )
    category_map = apply_scenario_config_to_categories(
        category_map=category_map,
        land_cover=land_cover,
        protected_areas=protected_areas,
        scenario_config=scenario_config
    )
    onshore_areas = pd.DataFrame(
        index=unit_ids,
        data={
            eligibility_category.area_column_name: _area(
                eligibility_category=eligibility_category,
                area_map=area_map,
                category_map=category_map,
                shapes=unit_geometries,
                transform=transform
            )
            for eligibility_category in Eligibility.onshore()
        }
    )
    offshore_eez_areas = pd.DataFrame(
        index=eez_ids,
        data={
            eligibility_category.area_column_name: _area(
                eligibility_category=eligibility_category,
                area_map=area_map,
                category_map=category_map,
                shapes=eez_geometries,
                transform=transform
            )
            for eligibility_category in Eligibility.offshore()
        }
    )
    offshore_areas = pd.DataFrame(
        data=shared_coasts.dot(offshore_eez_areas),
        columns=[cat.area_column_name for cat in Eligibility.offshore()]
    )
    areas = pd.concat([onshore_areas, offshore_areas], axis=1)
    areas.index.name = "id"
    areas.to_csv(
        path_to_result,
        header=True,
        index=True
    )


def apply_scenario_config_to_areas(area_map, category_map,
                                   land_cover, protected_areas, scenario_config):
    """Limit eligible area of each pixel based on scenario config."""

    # share-rooftops-used
    share_rooftops_used = scenario_config["share-rooftops-used"]
    mask = category_map == Eligibility.ROOFTOP_PV
    area_map[mask] = area_map[mask] * share_rooftops_used

    # share-forest-used-for-wind
    share_forest_used_for_wind = scenario_config["share-forest-used-for-wind"]
    mask = np.isin(land_cover, FOREST) & (category_map != Eligibility.ROOFTOP_PV)
    area_map[mask] = area_map[mask] * share_forest_used_for_wind

    # share-other-land-used
    share_other_land_used = scenario_config["share-other-land-used"]
    mask = np.isin(land_cover, OTHER) & (category_map != Eligibility.ROOFTOP_PV)
    area_map[mask] = area_map[mask] * share_other_land_used

    # share-farmland-used
    share_farmland_used = scenario_config["share-farmland-used"]
    mask = np.isin(land_cover, FARM) & (category_map != Eligibility.ROOFTOP_PV)
    area_map[mask] = area_map[mask] * share_farmland_used

    # share-offshore-used
    share_offshore_used = scenario_config["share-offshore-used"]
    mask = category_map == Eligibility.OFFSHORE_WIND
    area_map[mask] = area_map[mask] * share_offshore_used

    return area_map


def apply_scenario_config_to_categories(category_map, land_cover, protected_areas, scenario_config):
    """Change categories of each pixel based on scenario config."""

    # pv-on-farmland
    pv_on_farmland = scenario_config["pv-on-farmland"]
    if not pv_on_farmland:
        mask = np.isin(land_cover, FARM) & (category_map == Eligibility.ONSHORE_WIND_AND_PV)
        category_map[mask] = Eligibility.ONSHORE_WIND

    # protected-areas-used
    use_protected_areas = scenario_config["use-protected-areas"]
    if not use_protected_areas:
        mask = (protected_areas == ProtectedArea.PROTECTED) & (category_map != Eligibility.ROOFTOP_PV)
        category_map[mask] = Eligibility.NOT_ELIGIBLE

    return category_map


def _area(eligibility_category, area_map, category_map, shapes, transform):
    """Determine eligible area of one eligibility category per shape."""
    area_map = area_map.copy()
    area_map[category_map != eligibility_category] = 0
    potentials = zonal_stats(
        shapes,
        area_map,
        affine=transform,
        stats="sum",
        nodata=-999
    )
    return [stat["sum"] for stat in potentials]


if __name__ == "__main__":
    areas(
        path_to_units=snakemake.input.units,
        path_to_eez=snakemake.input.eez,
        path_to_shared_coast=snakemake.input.shared_coast,
        path_to_eligible_area=snakemake.input.area,
        path_to_eligibility_categories=snakemake.input.category,
        path_to_land_cover=snakemake.input.land_cover,
        path_to_protected_areas=snakemake.input.protected_areas,
        scenario_config=snakemake.params.scenario,
        path_to_result=snakemake.output[0]
    )
