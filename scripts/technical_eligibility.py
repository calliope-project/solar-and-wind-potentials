"""This module determines an upper bound of land eligibility for renewable generation based on geospatial data.

In here, we only exclude areas based on technical restrictions.
"""

import numpy as np
import rasterio

DATATYPE = np.uint8

from renewablepotentialslib.eligibility import Eligibility, FARM, FOREST, OTHER, WATER


def determine_eligibility(path_to_land_cover, path_to_slope, path_to_bathymetry,
                          path_to_building_share, path_to_urban_green_share, path_to_result,
                          max_slope, max_building_share, max_urban_green_share, max_depth_offshore):
    """Determines eligibility of land for renewables."""
    with rasterio.open(path_to_land_cover) as src:
        transform = src.transform
        land_cover = src.read(1)
        crs = src.crs
    with rasterio.open(path_to_slope) as src:
        slope = src.read(1)
    with rasterio.open(path_to_bathymetry) as src:
        bathymetry = src.read(1)
    with rasterio.open(path_to_building_share) as src:
        building_share = src.read(1)
    with rasterio.open(path_to_urban_green_share) as src:
        urban_green_share = src.read(1)
    eligibility = _determine_eligibility(
        land_cover=land_cover,
        slope=slope,
        bathymetry=bathymetry,
        building_share=building_share,
        urban_green_share=urban_green_share,
        max_slope=max_slope,
        max_building_share=max_building_share,
        max_urban_green_share=max_urban_green_share,
        max_depth_offshore=max_depth_offshore
    )
    with rasterio.open(path_to_result, 'w', driver='GTiff', height=eligibility.shape[0],
                       width=eligibility.shape[1], count=1, dtype=DATATYPE,
                       crs=crs, transform=transform) as new_geotiff:
        new_geotiff.write(eligibility, 1)


def _determine_eligibility(
    land_cover, slope, bathymetry, building_share, urban_green_share,
    max_slope, max_building_share, max_urban_green_share, max_depth_offshore
):
    # parameters
    max_slope_pv = max_slope["pv"]
    max_slope_wind = max_slope["wind"]
    assert max_slope_pv <= max_slope_wind # wind can be built whereever pv can be built

    # prepare masks
    settlements = (building_share > max_building_share) | (urban_green_share > max_urban_green_share)
    farm = np.isin(land_cover, FARM)
    forest = np.isin(land_cover, FOREST)
    other = np.isin(land_cover, OTHER)
    water = np.isin(land_cover, WATER)
    pv = (slope <= max_slope_pv) & ~settlements & (farm | other)
    wind = (slope <= max_slope_wind) & ~settlements & (farm | forest | other)
    offshore = (bathymetry > max_depth_offshore) & ~settlements & water

    # allocate eligibility
    land = np.ones_like(land_cover, dtype=DATATYPE) * Eligibility.NOT_ELIGIBLE
    _add_eligibility(land, Eligibility.ROOFTOP_PV, settlements)
    _add_eligibility(land, Eligibility.ONSHORE_WIND_AND_PV, wind & pv)
    _add_eligibility(land, Eligibility.ONSHORE_WIND, wind & ~pv)
    _add_eligibility(land, Eligibility.OFFSHORE_WIND, offshore)
    return land


def _add_eligibility(land, eligibility, mask):
    assert all(land[mask] == Eligibility.NOT_ELIGIBLE), f"Overwriting other eligibility with {eligibility}."
    land[mask] = eligibility


if __name__ == "__main__":
    determine_eligibility(
        path_to_land_cover=snakemake.input.land_cover,
        path_to_slope=snakemake.input.slope,
        path_to_bathymetry=snakemake.input.bathymetry,
        path_to_building_share=snakemake.input.building_share,
        path_to_urban_green_share=snakemake.input.urban_green_share,
        max_slope=snakemake.params.max_slope,
        max_building_share=snakemake.params.max_building_share,
        max_urban_green_share=snakemake.params.max_urban_green_share,
        max_depth_offshore=snakemake.params.max_depth_offshore,
        path_to_result=snakemake.output[0]
    )
