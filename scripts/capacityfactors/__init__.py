import math

import numpy as np
import geopandas as gpd
import shapely

from src.conversion import transform_bounds

# from https://epsg.io/3035
EPSG_3035 = "EPSG:3035"
EPSG_3035_PROJ4 = "+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs "
# from https://epsg.io/4326
WGS84 = "EPSG:4326"
WGS84_PROJ4 = "+proj=longlat +datum=WGS84 +no_defs "


def point_raster_on_shapes(bounds_wgs84, resolution_km2, shapes):
    """Creates a point raster with given resolution on a set of shapes.

    Extends (=buffers) the shapes, so that whenever a raster cell is touched by any shape,
    a point is created for that cell.

    Parameters:
        * bounds_wgs84: the bounds of the point raster, given in WGS84
        * resolution_km2: the resolution of the point raster, given in km2
        * shapes: GeoDataFrame containing the shapes
    Returns:
        * point raster in WGS84 with given resolution, filtered by the shapes
    """
    x_min, y_min, x_max, y_max = transform_bounds(
        bounds_wgs84["x_min"], bounds_wgs84["y_min"], bounds_wgs84["x_max"], bounds_wgs84["y_max"],
        from_epsg=WGS84,
        to_epsg=EPSG_3035
    )
    all_points = [
        shapely.geometry.Point(x, y)
        for x in np.arange(start=x_min, stop=x_max, step=resolution_km2 * 1000)
        for y in np.arange(start=y_min, stop=y_max, step=resolution_km2 * 1000)
    ]
    simplification_strength = resolution_km2 * 1000 / 20
    buffer_size = math.sqrt(resolution_km2 ** 2 + resolution_km2 ** 2) / 2 * 1000
    surface_areas = (shapes.to_crs(EPSG_3035_PROJ4)
                           .simplify(simplification_strength)
                           .buffer(buffer_size))
    prepared_polygons = [shapely.prepared.prep(polygon) for polygon in surface_areas.geometry]
    return gpd.GeoSeries(
        list(filter(
            lambda point: any(polygon.intersects(point) for polygon in prepared_polygons),
            all_points
        )),
        crs=EPSG_3035_PROJ4
    ).to_crs(WGS84_PROJ4)
