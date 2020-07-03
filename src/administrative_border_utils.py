import geopandas as gpd
import shapely.geometry
import shapely.errors
import pycountry

SCHEMA = {
    "properties": {
        "country_code": "str", # id of the country to which the unit belongs
        "id": "str", # a unique id of this unit
        "name": "str", # the name of the unit, not necessarily unique
        "type": "str", # the type of the unit
        "proper": "bool" # flag indicating proper administrative unit (not the case for water bodies e.g.)
    },
    "geometry": "MultiPolygon"
}


def drop_countries(gdf, config):
    countries = [pycountry.countries.lookup(i).alpha_3 for i in config["scope"]["countries"]]
    _not_in = set(gdf.country_code).difference(countries)
    print(f"Removing {_not_in} as they are outside of study area.")

    return gdf[gdf.country_code.isin(countries)]


def drop_geoms(gdf, config):
    study_area = _study_area(config)
    completely_in = gdf.intersects(study_area)
    for i in gdf[~completely_in].iterrows():
        print(
            "Removing {} ({}) as they are outside of study area."
            .format(*i[1][["name", "level"]])
        )
    gdf = gdf[completely_in]

    all_geoms = gdf.explode()
    partially_in = all_geoms.within(study_area)
    partially_out = ~partially_in.groupby(level=0).min()
    for i in gdf.loc[partially_out].iterrows():
        print(
            "Removing parts of {} ({}) as they are outside of study area."
            .format(*i[1][["name", "level"]])
        )
    # Unlike groupby, dissolve can only operate on columns, not multiindex levels
    new_geoms = all_geoms[partially_in.mul(partially_out, level=0)].reset_index().dissolve('level_0')
    gdf.update(new_geoms.geometry.map(_to_multi_polygon).to_frame('geometry'))

    return gdf


def to_file(gdf, path_to_output, output_driver):
    assert gdf.id.duplicated().sum() == 0
    allowed_cols = list(SCHEMA["properties"].keys()) + ['geometry']
    if output_driver.lower() == 'gpkg':
        for lvl in gdf.level.unique():
            gdf.loc[gdf.level == lvl, allowed_cols].to_file(
                path_to_output, schema=SCHEMA, layer=lvl, driver=output_driver
            )
    else:
        gdf[allowed_cols].to_file(path_to_output, schema=SCHEMA, driver=output_driver)


def _to_multi_polygon(geometry):
    if isinstance(geometry, dict):
        geometry = shapely.geometry.shape(geometry)
    if isinstance(geometry, shapely.geometry.polygon.Polygon):
        return shapely.geometry.MultiPolygon(polygons=[geometry])
    else:
        return geometry


def _study_area(config):
    """
    Create a bounding box for the study area, and cut out holes for all defined
    exclusion zones. For plotting purposes, exclusion zones and the bounding box are
    defined in opposite orientations, see https://github.com/geopandas/geopandas/issues/951
    """
    if config["scope"].get("exclusion_zones", {}) and isinstance(config["scope"]["exclusion_zones"], dict):
        holes = [
            (
                (exclusion_zone["x_max"], exclusion_zone["y_min"]),
                (exclusion_zone["x_max"], exclusion_zone["y_max"]),
                (exclusion_zone["x_min"], exclusion_zone["y_max"]),
                (exclusion_zone["x_min"], exclusion_zone["y_min"])
            )
            for exclusion_zone in config["scope"]["exclusion_zones"].values()
        ]
    else:
        holes = []

    study_area = shapely.geometry.Polygon(
        ((config["scope"]["bounds"]["x_min"], config["scope"]["bounds"]["y_min"]),
        (config["scope"]["bounds"]["x_min"], config["scope"]["bounds"]["y_max"]),
        (config["scope"]["bounds"]["x_max"], config["scope"]["bounds"]["y_max"]),
        (config["scope"]["bounds"]["x_max"], config["scope"]["bounds"]["y_min"])),
        holes=holes
    )
    if study_area.is_valid is False:
        study_area = study_area.buffer(0)

    return study_area