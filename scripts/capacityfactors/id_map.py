"""Create maps of ids to capacity factor timeseries of renewables."""
import math

import numpy as np
import geopandas as gpd
import shapely
import rasterio
from rasterio.transform import from_origin
import xarray as xr

DTYPE = np.uint16
NO_DATA_VALUE = 64001
INDEX_EPSILON = 10e-3

from renewablepotentialslib.geo import EPSG_3035, EPSG_3035_PROJ4, WGS84, WGS84_PROJ4


def id_map(path_to_timeseries, path_to_map, resolution_km):
    """Create maps of ids to capacity factor timeseries of renewables.

    Each point on the map links to a timeseries of capacity factors of renewables. Together with the
    timeseries, both files form the spatio-temporal data format used in this study.
    """
    ds = xr.open_dataset(path_to_timeseries)
    pv_config = ds[["lat", "lon"]].to_dataframe()
    points = gpd.GeoDataFrame(
        geometry=[shapely.geometry.Point(row.lon, row.lat) for _, row in pv_config.iterrows()],
        data={
            "site_id": pv_config.index
        },
        crs=WGS84_PROJ4
    ).to_crs(EPSG_3035_PROJ4)
    x_min = min([point.x for point in points.geometry])
    x_max = max([point.x for point in points.geometry])
    y_min = min([point.y for point in points.geometry])
    y_max = max([point.y for point in points.geometry])
    resolution_m = resolution_km * 1000
    width = (x_max - x_min) / resolution_m + 1
    height = (y_max - y_min) / resolution_m + 1
    assert isclose(round(width), width) # diff is purely numerics
    assert isclose(round(height), height) # diff is purely numerics
    width = round(width)
    height = round(height)
    raster = np.ones(shape=(height, width), dtype=DTYPE) * NO_DATA_VALUE
    for _, point in points.iterrows():
        index_x = (point.geometry.x - x_min) / resolution_m
        index_y = (y_max - point.geometry.y) / resolution_m
        assert isclose(round(index_x), index_x) # diff is purely numerics
        assert isclose(round(index_y), index_y) # diff is purely numerics
        int_index_x = round(index_x)
        int_index_y = round(index_y)
        raster[int_index_y, int_index_x] = point.site_id
    transform = from_origin(
        west=x_min - resolution_m / 2,
        north=y_max + resolution_m / 2,
        xsize=resolution_m,
        ysize=resolution_m
    )
    with rasterio.open(path_to_map, 'w', driver='GTiff', height=height, width=width,
                       count=1, dtype=DTYPE, crs=EPSG_3035, transform=transform,
                       nodata=NO_DATA_VALUE) as f_map:
        f_map.write(raster, 1)


def isclose(a, b):
    return math.isclose(a, b, abs_tol=INDEX_EPSILON, rel_tol=0)


if __name__ == "__main__":
    id_map(
        path_to_timeseries=snakemake.input.timeseries,
        resolution_km=snakemake.params.resolution,
        path_to_map=snakemake.output[0]
    )
