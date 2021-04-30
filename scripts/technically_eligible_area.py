"""Determines area of technically eligible land for renewables."""
import numpy as np
import rasterio

from renewablepotentialslib.geo.shape_utils import determine_pixel_areas
from renewablepotentialslib.eligibility import Eligibility

DATATYPE = np.float32


def determine_area(path_to_eligibility_categories, path_to_building_share,
                   path_to_rooftop_correction_factor, path_to_result):
    """Determines area of technically eligible land for renewables.

    For all eligibility categories other than rooftop PV, this is simply the pixel/cell size.
    For rooftop PV, we reduce the area to the building footprints, and we furthermore apply a
    correction factor to map from building footprint to available rooftop space.
    """
    with rasterio.open(path_to_eligibility_categories) as src:
        eligibility_categories = src.read(1)
        meta = src.meta
        bounds = src.bounds
        resolution = src.res[0]
    pixel_area = determine_pixel_areas(meta["crs"], bounds, resolution).astype(DATATYPE)
    areas_of_eligibility = pixel_area.copy()
    rooftop_area = determine_rooftop_areas(pixel_area, path_to_building_share, path_to_rooftop_correction_factor)
    pv_rooftop_mask = eligibility_categories == Eligibility.ROOFTOP_PV
    areas_of_eligibility[pv_rooftop_mask] = rooftop_area[pv_rooftop_mask]
    write_to_file(areas_of_eligibility, path_to_result, meta)


def determine_rooftop_areas(pixel_areas, path_to_building_share, path_to_rooftop_correction_factor):
    """Returns a raster in which the value corresponds to the rooftop area in the pixel."""
    with rasterio.open(path_to_building_share) as f_building_share, \
            open(path_to_rooftop_correction_factor, "r") as f_factor:
        factor = float(f_factor.readline())
        building_share = f_building_share.read(1)
    return pixel_areas * building_share * factor


def write_to_file(areas_of_eligibility, path_to_result, meta):
    meta.update(dtype=DATATYPE)
    with rasterio.open(path_to_result, 'w', **meta) as new_geotiff:
        new_geotiff.write(areas_of_eligibility, 1)


if __name__ == "__main__":
    determine_area(
        path_to_eligibility_categories=snakemake.input.eligibility_categories,
        path_to_building_share=snakemake.input.building_share,
        path_to_rooftop_correction_factor=snakemake.input.rooftop_correction_factor,
        path_to_result=snakemake.output[0]
    )
