import click
import geopandas as gpd
import shapely.geometry
import shapely.errors
import pycountry

from conversion import eu_country_code_to_iso3
from utils import Config

OUTPUT_DRIVER = 'GPKG'
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


@click.command()
@click.argument("path_to_nuts")
@click.argument("path_to_gadm")
@click.argument("path_to_lau")
@click.argument("path_to_output")
@click.argument("config", type=Config())
def normalise_admin_borders(path_to_nuts, path_to_gadm, path_to_lau, path_to_output, config):
    """Normalises raw administrative boundary data and places it in one, layered geodatapackage."""

    #lau
    for _src, _path in {
        'nuts': path_to_nuts, 'gadm': path_to_gadm, 'lau': path_to_lau
    }.items():
        gdf = gpd.read_file(_path)
        gdf = gdf.to_crs(config["crs"])
        gdf.geometry = gdf.geometry.buffer(0).map(_to_multi_polygon)
        gdf = _update_features(gdf, _src)
        gdf = _drop_countries(gdf, config)
        gdf = _drop_geoms(gdf, config)

        assert gdf.id.duplicated().sum() == 0

        allowed_cols = list(SCHEMA["properties"].keys()) + ['geometry']

        for lvl in gdf.level.unique():
            gdf.loc[gdf.level == lvl, allowed_cols].to_file(
                path_to_output, schema=SCHEMA, layer=lvl, driver=OUTPUT_DRIVER
            )


def _update_features(gdf, src):
    if src == 'nuts':
        gdf["CNTR_CODE"] = gdf.CNTR_CODE.apply(eu_country_code_to_iso3)
        gdf = gdf.rename(
            columns={"NUTS_ID": "id", "NUTS_NAME": "name", "CNTR_CODE": "country_code"}
        )
        gdf["type"] = gdf.LEVL_CODE.map({0: "country"})
        gdf["proper"] = True
        gdf = gdf.drop(columns=["FID"])

        # Country IDs should have three letters instead of two
        gdf.loc[gdf.LEVL_CODE == 0, "id"] = gdf.loc[gdf.LEVL_CODE == 0, "country_code"]
        # Many country names are broken or missing
        gdf.loc[gdf.LEVL_CODE == 0, "name"] = gdf.loc[gdf.LEVL_CODE == 0, "id"].apply(
            lambda x: pycountry.countries.lookup(x).name
        )

        gdf["level"] = 'nuts' + gdf.LEVL_CODE

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
        gdf["level"] = 'gadm' + gdf.level

    elif src == 'lau':
        gdf["CNTR_CODE"] = gdf.COMM_ID.str[:2].apply(eu_country_code_to_iso3)
        gdf = gdf.rename(
            columns={"COMM_ID": "id", "NAME_LATN": "name", "CNTR_CODE": "country_code"}
        )
        gdf["type"] = "commune"
        gdf["proper"] = gdf["TRUE_COMM_"] == "T"
        gdf["level"] = "lau2"

    return gdf


def _drop_countries(gdf, config):
    countries = [pycountry.countries.lookup(i).alpha_3 for i in config["scope"]["countries"]]
    _not_in = set(gdf.country_code).difference(countries)
    print(f"Removing {_not_in} as they are outside of study area.")

    return gdf[gdf.country_code.isin(countries)]


def _drop_geoms(gdf, config):
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


if __name__ == "__main__":
    normalise_admin_borders()