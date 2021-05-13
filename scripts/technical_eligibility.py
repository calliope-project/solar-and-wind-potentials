"""This module determines an upper bound of land eligibility for renewable generation based on geospatial data.

In here, we only exclude areas based on technical restrictions.
"""

import numpy as np
import rasterio

from renewablepotentialslib.eligibility import eligibility_land_mask, DATATYPE


def determine_eligibility(path_to_land_cover, path_to_integer_slope, path_to_bathymetry,
                          path_to_building_share, path_to_urban_green_share, path_to_result,
                          max_slope, max_building_share, max_urban_green_share, max_depth_offshore):
    """Determines eligibility of land for renewables."""
    with rasterio.open(path_to_land_cover) as src:
        transform = src.transform
        land_cover = src.read(1)
        crs = src.crs
    with rasterio.open(path_to_integer_slope) as src:
        # see https://land.copernicus.eu/user-corner/technical-library/slope-conversion-table
        slope = np.degrees(np.arccos(src.read(1) / 250))
    with rasterio.open(path_to_bathymetry) as src:
        bathymetry = src.read(1)
    with rasterio.open(path_to_building_share) as src:
        building_share = src.read(1)
    with rasterio.open(path_to_urban_green_share) as src:
        urban_green_share = src.read(1)
    eligibility = eligibility_land_mask(
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


if __name__ == "__main__":
    determine_eligibility(
        path_to_land_cover=snakemake.input.land_cover,
        path_to_slope=snakemake.input.slope,
        path_to_bathymetry=snakemake.input.bathymetry,
        path_to_building_share=snakemake.input.building_share,
        path_to_urban_green_share=snakemake.input.urban_green_share,
        max_integer_slope=snakemake.params.max_integer_slope,
        max_building_share=snakemake.params.max_building_share,
        max_urban_green_share=snakemake.params.max_urban_green_share,
        max_depth_offshore=snakemake.params.max_depth_offshore,
        path_to_result=snakemake.output[0]
    )
