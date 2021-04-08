import click
import geopandas as gpd
import shapely.geometry
import shapely.errors
import pycountry

from src.conversion import eu_country_code_to_iso3
from src.utils import Config, buffer_if_necessary

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
        gdf.geometry = gdf.geometry.map(buffer_if_necessary).map(_to_multi_polygon)
        gdf = _update_features(gdf, _src)
        gdf = _drop_countries(gdf, config)
        gdf = _drop_geoms_completely_outside_study_area(gdf, config)
        gdf = _drop_parts_of_geoms_completely_outside_study_area(gdf, config)

        assert gdf.id.duplicated().sum() == 0

        allowed_cols = list(SCHEMA["properties"].keys()) + ['geometry']

        for lvl in gdf.level.unique():
            gdf.loc[gdf.level == lvl, allowed_cols].to_file(
                path_to_output, schema=SCHEMA, layer=lvl, driver=OUTPUT_DRIVER
            )


def _update_features(gdf, src):
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


def _drop_countries(gdf, config, print_dropped=True):
    countries = [pycountry.countries.lookup(i).alpha_3 for i in config["scope"]["countries"]]
    _not_in = set(gdf.country_code).difference(countries)
    if print_dropped:
        print(f"Removing {_not_in} as they are outside of study area.")

    return gdf[gdf.country_code.isin(countries)]


def _drop_geoms_completely_outside_study_area(gdf, config, print_dropped=True):
    study_area = _study_area(config)
    completely_in = gdf.intersects(study_area)
    if print_dropped:
        for row_index, row in gdf[~completely_in].iterrows():
            print(
                "Removing {} ({}, country={}) as they are outside of study area."
                .format(*row[["name", "level", "country_code"]])
            )
    gdf = gdf[completely_in]

    return gdf


def _drop_parts_of_geoms_completely_outside_study_area(gdf, config, print_dropped=True):
    gdf = gdf.copy()
    study_area = _study_area(config)
    all_geoms = gdf.explode()
    geoms_within_study_area = all_geoms.within(study_area)
    geoms_partially_out = ~geoms_within_study_area.all(level=0)

    # work only with geoms which have some polygons within the study area and some out
    geoms_to_update = geoms_within_study_area.mul(geoms_partially_out, level=0)
    if gdf.loc[geoms_to_update.any(level=0)].empty:
        return gdf
    if print_dropped:
        for row_index, row in gdf.loc[geoms_to_update.any(level=0)].iterrows():
            print(
                "Removing parts of {} ({}, country={}) as they are outside of study area."
                .format(*row[["name", "level", "country_code"]])
            )

    new_geoms = (
        all_geoms[geoms_to_update]
        .geometry
        .map(buffer_if_necessary)
        .groupby(level=0)
        .apply(lambda x: x.unary_union)
        .map(_to_multi_polygon)
    )
    gdf.loc[new_geoms.index, 'geometry'] = new_geoms

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
    holes = [
        (
            (exclusion_zone["x_max"], exclusion_zone["y_min"]),
            (exclusion_zone["x_max"], exclusion_zone["y_max"]),
            (exclusion_zone["x_min"], exclusion_zone["y_max"]),
            (exclusion_zone["x_min"], exclusion_zone["y_min"])
        )
        for exclusion_zone in config["scope"].get("exclusion_zones", {}).values()
    ]

    study_area = shapely.geometry.Polygon(
        ((config["scope"]["bounds"]["x_min"], config["scope"]["bounds"]["y_min"]),
        (config["scope"]["bounds"]["x_min"], config["scope"]["bounds"]["y_max"]),
        (config["scope"]["bounds"]["x_max"], config["scope"]["bounds"]["y_max"]),
        (config["scope"]["bounds"]["x_max"], config["scope"]["bounds"]["y_min"])),
        holes=holes
    )
    study_area = buffer_if_necessary(study_area)

    return study_area


if __name__ == "__main__":
    normalise_admin_borders()
