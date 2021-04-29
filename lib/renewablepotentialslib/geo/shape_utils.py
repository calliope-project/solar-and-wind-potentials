"""Module containing utilities."""
import math

import numpy as np
import geopandas as gpd
import shapely
import pycountry
import pyproj

from renewablepotentialslib.geo.conversion import transform_bounds, eu_country_code_to_iso3
from renewablepotentialslib.geo import EPSG_3035, EPSG_3035_PROJ4, WGS84, WGS84_PROJ4


def determine_pixel_areas(crs, bounds, resolution):
    """Returns a raster in which the value corresponds to the area [km2] of the pixel.

    This assumes the data comprises square pixel in WGS84.

    Parameters:
        crs: the coordinate reference system of the data (must be WGS84)
        bounds: an object with attributes left/right/top/bottom given in degrees
        resolution: the scalar resolution (remember: square pixels) given in degrees
    """
    # the following is based on https://gis.stackexchange.com/a/288034/77760
    # and assumes the data to be in WGS84
    assert pyproj.crs.CRS(crs) == pyproj.crs.CRS(WGS84)
    width = int((bounds.right - bounds.left) / resolution)
    height = int((bounds.top - bounds.bottom) / resolution)
    latitudes = np.linspace(
        start=bounds.top,
        stop=bounds.bottom,
        num=height,
        endpoint=True,
        dtype=np.float64
    )
    varea_of_pixel = np.vectorize(lambda lat: _area_of_pixel(resolution, lat))
    pixel_area = varea_of_pixel(latitudes) # vector
    return pixel_area.repeat(width).reshape(height, width).astype(np.float64)


def _area_of_pixel(pixel_size, center_lat):
    """Calculate km^2 area of a wgs84 square pixel.

    Adapted from: https://gis.stackexchange.com/a/127327/2397

    Parameters:
        pixel_size (float): length of side of pixel in degrees.
        center_lat (float): latitude of the center of the pixel. Note this
            value +/- half the `pixel-size` must not exceed 90/-90 degrees
            latitude or an invalid area will be calculated.

    Returns:
        Area of square pixel of side length `pixel_size` centered at
        `center_lat` in km^2.

    """
    a = 6378137  # meters
    b = 6356752.3142  # meters
    e = math.sqrt(1 - (b / a)**2)
    area_list = []
    for f in [center_lat + pixel_size / 2, center_lat - pixel_size / 2]:
        zm = 1 - e * math.sin(math.radians(f))
        zp = 1 + e * math.sin(math.radians(f))
        area_list.append(
            math.pi * b**2 * (
                math.log(zp / zm) / (2 * e) +
                math.sin(math.radians(f)) / (zp * zm)))
    return pixel_size / 360. * (area_list[0] - area_list[1]) / 1e6


def buffer_if_necessary(shape):
    """Fix shapes which are invalid.

    Following the advice given here:
    https://github.com/Toblerity/Shapely/issues/344
    """
    if shape.is_valid:
        return shape

    new_shape = shape.buffer(0.0)
    assert new_shape.is_valid
    assert np.isclose(new_shape.area, shape.area, rtol=1e-5)

    return new_shape


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


def update_features(gdf, src):
    if src == 'nuts':
        gdf["CNTR_CODE"] = gdf.CNTR_CODE.apply(eu_country_code_to_iso3)
        gdf = gdf.rename(columns={"NUTS_NAME": "name", "CNTR_CODE": "country_code"})
        gdf["type"] = gdf.LEVL_CODE.map({0: "country"})
        gdf["proper"] = True
        gdf = gdf.drop(columns=["FID"])

        # Country IDs should have three letters instead of two
        gdf.loc[gdf.LEVL_CODE == 0, "id"] = gdf.loc[gdf.LEVL_CODE == 0, "country_code"]
        # Many country names are broken or missing
        gdf.loc[gdf.LEVL_CODE == 0, "name"] = gdf.loc[gdf.LEVL_CODE == 0, "id"].apply(
            lambda x: pycountry.countries.lookup(x).name
        )

        gdf["level"] = 'nuts' + gdf.LEVL_CODE.astype(str)

    elif src == 'gadm':
        gdf["level"] = gdf.source_ds_lyr.str.rsplit('_', 1, expand=True)[1].astype(int)
        for lvl in gdf.level.unique():
            lvl_mask = gdf.level == lvl
            gdf.loc[lvl_mask, "country_code"] = gdf.loc[lvl_mask, "GID_0"]
            gdf.loc[lvl_mask, "id"] = gdf.loc[lvl_mask, f"GID_{lvl}"]
            gdf.loc[lvl_mask, "name"] = gdf.loc[lvl_mask, f"NAME_{lvl}"]
            if lvl > 0:
                gdf.loc[lvl_mask, "type"] = gdf.loc[lvl_mask, f"ENGTYPE_{lvl}"]
            else:
                gdf.loc[lvl_mask, "type"] = "country"
        gdf["proper"] = True
        gdf["level"] = 'gadm' + gdf.level.astype(str)

    elif src == 'lau':
        gdf["CNTR_CODE"] = gdf.COMM_ID.str[:2].apply(eu_country_code_to_iso3)
        gdf = gdf.rename(
            columns={"COMM_ID": "id", "NAME_LATN": "name", "CNTR_CODE": "country_code"}
        )
        gdf["type"] = "commune"
        gdf["proper"] = gdf["TRUE_COMM_"] == "T"
        gdf["level"] = "lau2"

    return gdf


def drop_countries(gdf, scope_config):
    countries = [pycountry.countries.lookup(i).alpha_3 for i in scope_config["countries"]]
    _not_in = set(gdf.country_code).difference(countries)
    if _not_in:
        print(f"Removing {_not_in} as they are outside of study area.")

    return gdf[gdf.country_code.isin(countries)]


def drop_geoms_completely_outside_study_area(gdf, scope_config):
    _study_area = study_area(scope_config)
    completely_in = gdf.intersects(_study_area)
    for row_index, row in gdf[~completely_in].iterrows():
        print(
            "Removing {} ({}, country={}) as they are outside of study area."
            .format(*row[["name", "level", "country_code"]])
        )
    gdf = gdf[completely_in]

    return gdf


def drop_parts_of_geoms_completely_outside_study_area(gdf, scope_config):
    gdf = gdf.copy()
    _study_area = study_area(scope_config)
    all_geoms = gdf.explode()
    geoms_within_study_area = all_geoms.within(_study_area)
    geoms_partially_out = ~geoms_within_study_area.all(level=0)

    # work only with geoms which have some polygons within the study area and some out
    geoms_to_update = geoms_within_study_area.mul(geoms_partially_out, level=0)
    if gdf.loc[geoms_to_update.any(level=0)].empty:
        return gdf

    for row_index, row in gdf.loc[geoms_to_update.any(level=0)].iterrows():
        print(
            "Removing parts of {} ({}, country={}) as they are outside of study area."
            .format(*row[["name", "level", "country_code"]])
        )
    # Unlike groupby, dissolve can only operate on columns, not multiindex levels

    new_geoms = (
        all_geoms[geoms_to_update]
        .geometry
        .map(buffer_if_necessary)
        .groupby(level=0)
        .apply(lambda x: x.unary_union)
        .map(to_multi_polygon)
    )
    gdf.loc[new_geoms.index, 'geometry'] = new_geoms

    return gdf


def to_multi_polygon(geometry):
    if isinstance(geometry, dict):
        geometry = shapely.geometry.shape(geometry)
    if isinstance(geometry, shapely.geometry.polygon.Polygon):
        return shapely.geometry.MultiPolygon(polygons=[geometry])
    else:
        return geometry


def study_area(scope_config):
    """
    Create a bounding box for the study area, and cut out holes for all defined
    exclusion zones. For plotting purposes, exclusion zones and the bounding box are
    defined in opposite orientations, see https://github.com/geopandas/geopandas/issues/951
    """
    holes = [
        (
            (exclusion_zone["x_max"], exclusion_zone["y_min"]),
            (exclusion_zone["x_max"], exclusion_zone["y_max"]),
            (exclusion_zone["x_min"], exclusion_zone["y_max"]),
            (exclusion_zone["x_min"], exclusion_zone["y_min"])
        )
        for exclusion_zone in scope_config.get("exclusion_zones", {}).values()
    ]

    _study_area = shapely.geometry.Polygon(
        ((scope_config["bounds"]["x_min"], scope_config["bounds"]["y_min"]),
        (scope_config["bounds"]["x_min"], scope_config["bounds"]["y_max"]),
        (scope_config["bounds"]["x_max"], scope_config["bounds"]["y_max"]),
        (scope_config["bounds"]["x_max"], scope_config["bounds"]["y_min"])),
        holes=holes
    )
    _study_area = buffer_if_necessary(_study_area)

    return _study_area
