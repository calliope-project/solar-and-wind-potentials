"""Determine the built up area in administrative units."""
import fiona
import rasterio
from rasterstats import zonal_stats
import pandas as pd

from renewablepotentialslib.geo.shape_utils import determine_pixel_areas


def built_up_areas(path_to_built_up_share, path_to_units, path_to_result):
    """Determine the built up area in administrative units."""
    with rasterio.open(path_to_built_up_share) as src:
        built_up_share = src.read(1)
        crs = src.crs
        transform = src.transform
        bounds = src.bounds
        resolution = src.res[0]
    with fiona.open(path_to_units, "r") as src:
        unit_ids = [feature["properties"]["id"] for feature in src]
        unit_geometries = [feature["geometry"] for feature in src]

    pixel_area = determine_pixel_areas(crs, bounds, resolution)
    built_up_stats = pd.DataFrame(
        index=unit_ids,
        data={
            "built_up_km2": _stats(unit_geometries, built_up_share * pixel_area, transform),
            "non_built_up_km2": _stats(unit_geometries, (1 - built_up_share) * pixel_area, transform)
        }
    )
    built_up_stats["built_up_share"] = (built_up_stats["built_up_km2"] /
                                        (built_up_stats["built_up_km2"] + built_up_stats["non_built_up_km2"]))
    built_up_stats.index.name = "id"
    built_up_stats.to_csv(
        path_to_result,
        header=True,
        index=True
    )


def _stats(unit_geometries, raster_area, transform):
    area = zonal_stats(
        unit_geometries,
        raster_area,
        affine=transform,
        stats="sum",
        nodata=-999
    )
    return [stat["sum"] for stat in area]


if __name__ == "__main__":
    built_up_areas(
        path_to_built_up_share=snakemake.input.built_up_area,
        path_to_units=snakemake.input.units,
        path_to_result=snakemake.output[0]
    )
