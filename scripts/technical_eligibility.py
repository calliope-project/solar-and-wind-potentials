"""This module determines an upper bound of land eligibility for renewable generation based on geospatial data.

In here, we only exclude areas based on technical restrictions.
"""
from enum import IntEnum

import numpy as np
import rasterio

DATATYPE = np.uint8


class Eligibility(IntEnum):
    """Categories defining land eligibility for renewable power."""
    NOT_ELIGIBLE = 0
    ROOFTOP_PV = 250
    ONSHORE_WIND_AND_PV = 180
    ONSHORE_WIND = 110
    OFFSHORE_WIND = 40

    @property
    def area_column_name(self):
        return "eligibility_{}_km2".format(self.name.lower())

    @staticmethod
    def onshore():
        """Returns all onshore eligibilities."""
        return [
            Eligibility.NOT_ELIGIBLE,
            Eligibility.ROOFTOP_PV,
            Eligibility.ONSHORE_WIND_AND_PV,
            Eligibility.ONSHORE_WIND,
        ]

    @staticmethod
    def offshore():
        """Returns all offshore eligibilities."""
        return [
            Eligibility.OFFSHORE_WIND
        ]


class GlobCover(IntEnum):
    """Original categories taken from GlobCover 2009 land cover."""
    POST_FLOODING = 11
    RAINFED_CROPLANDS = 14
    MOSAIC_CROPLAND = 20
    MOSAIC_VEGETATION = 30
    CLOSED_TO_OPEN_BROADLEAVED_FOREST = 40
    CLOSED_BROADLEAVED_FOREST = 50
    OPEN_BROADLEAVED_FOREST = 60
    CLOSED_NEEDLELEAVED_FOREST = 70
    OPEN_NEEDLELEAVED_FOREST = 90
    CLOSED_TO_OPEN_MIXED_FOREST = 100
    MOSAIC_FOREST = 110
    MOSAIC_GRASSLAND = 120
    CLOSED_TO_OPEN_SHRUBLAND = 130
    CLOSED_TO_OPEN_HERBS = 140
    SPARSE_VEGETATION = 150
    CLOSED_TO_OPEN_REGULARLY_FLOODED_FOREST = 160 # doesn't exist in Europe
    CLOSED_REGULARLY_FLOODED_FOREST = 170 # doesn't exist in Europe
    CLOSED_TO_OPEN_REGULARLY_FLOODED_GRASSLAND = 180 # roughly 2.3% of land in Europe
    ARTIFICAL_SURFACES_AND_URBAN_AREAS = 190
    BARE_AREAS = 200
    WATER_BODIES = 210
    PERMANENT_SNOW = 220
    NO_DATA = 230


FARM = [GlobCover.POST_FLOODING, GlobCover.RAINFED_CROPLANDS,
        GlobCover.MOSAIC_CROPLAND, GlobCover.MOSAIC_VEGETATION]
FOREST = [GlobCover.CLOSED_TO_OPEN_BROADLEAVED_FOREST, GlobCover.CLOSED_BROADLEAVED_FOREST,
          GlobCover.OPEN_BROADLEAVED_FOREST, GlobCover.CLOSED_NEEDLELEAVED_FOREST,
          GlobCover.OPEN_NEEDLELEAVED_FOREST, GlobCover.CLOSED_TO_OPEN_MIXED_FOREST,
          GlobCover.MOSAIC_FOREST, GlobCover.CLOSED_TO_OPEN_REGULARLY_FLOODED_FOREST,
          GlobCover.CLOSED_REGULARLY_FLOODED_FOREST]
VEGETATION = [GlobCover.MOSAIC_GRASSLAND, GlobCover.CLOSED_TO_OPEN_SHRUBLAND,
              GlobCover.CLOSED_TO_OPEN_HERBS, GlobCover.SPARSE_VEGETATION,
              GlobCover.CLOSED_TO_OPEN_REGULARLY_FLOODED_GRASSLAND]
BARE = [GlobCover.BARE_AREAS]
OTHER = VEGETATION + BARE
URBAN = [GlobCover.ARTIFICAL_SURFACES_AND_URBAN_AREAS]
WATER = [GlobCover.WATER_BODIES]


def determine_eligibility(path_to_land_cover, path_to_slope,
                          path_to_bathymetry, path_to_building_share, path_to_urban_green_share,
                          path_to_result, max_slope, max_building_share, max_urban_green_share):
    """Determines eligibility of land for renewables."""
    with rasterio.open(str(path_to_land_cover)) as src:
        transform = src.transform
        land_cover = src.read(1)
        crs = src.crs
    with rasterio.open(str(path_to_slope)) as src:
        slope = src.read(1)
    with rasterio.open(str(path_to_bathymetry)) as src:
        bathymetry = src.read(1)
    with rasterio.open(str(path_to_building_share)) as src:
        building_share = src.read(1)
    with rasterio.open(str(path_to_urban_green_share)) as src:
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
    )
    with rasterio.open(str(path_to_result), 'w', driver='GTiff', height=eligibility.shape[0],
                       width=eligibility.shape[1], count=1, dtype=DATATYPE,
                       crs=crs, transform=transform) as new_geotiff:
        new_geotiff.write(eligibility, 1)


def _determine_eligibility(
    land_cover, slope, bathymetry, building_share, urban_green_share, max_slope, max_building_share, max_urban_green_share
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
    offshore = (bathymetry > config["parameters"]["max-depth-offshore"]) & ~settlements & water

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
        path_to_result=snakemake.output
    )
