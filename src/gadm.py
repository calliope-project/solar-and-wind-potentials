"""Module to merge and preprocess GADM administrative borders."""
from itertools import chain

import click
import fiona
import fiona.transform
import geopandas as gpd
import shapely.geometry
import shapely.errors
from shapely.prepared import prep

from src.utils import Config

LAYER_NAME = "gadm{layer_id}"
SCHEMA = {
    "properties": {
        "country_code": "str", # id of the country to which the unit belongs
        "id": "str", # a unique id of this unit
        "name": "str", # the name of the unit, not necessarily unqique
        "type": "str", # the type of the unit
        "proper": "int" # flag indicating proper administrative unit (not the case for water bodies e.g.)
        # FIXME should be bool, but GeoPandas doesn't support bool: https://github.com/geopandas/geopandas/issues/437
        # will be fixed here: https://github.com/geopandas/geopandas/pull/855
    },
    "geometry": "MultiPolygon"
}


@click.command()
@click.argument("path_to_countries", nargs=-1, metavar="COUNTRIES...")
@click.argument("max_layer_depths", type=click.INT)
@click.argument("path_to_output")
@click.argument("config", type=Config())
def retrieve_administrative_borders(path_to_countries, max_layer_depths, path_to_output, config):
    study_area = _study_area(config)
    study_area = prep(study_area) # improves performance
    with fiona.open(path_to_countries[0], "r", layer=0) as first_country:
        src_crs = first_country.crs
        src_driver = first_country.driver
    for layer_id in range(max_layer_depths + 1):
        print(f"Merging layer {layer_id}...")
        with fiona.open(path_to_output,
                        "w",
                        crs=config["crs"],
                        schema=SCHEMA,
                        driver=src_driver,
                        layer=LAYER_NAME.format(layer_id=layer_id)) as merged_file:
            merged_file.writerecords(
                [_reproject(feature, src_crs, config["crs"]) for feature in chain(
                    *[_country_features(path_to_country, layer_id, study_area)
                      for path_to_country in path_to_countries])
                 ]
            )
    _test_id_uniqueness(path_to_output)


def _country_features(path_to_file, layer_id, study_area):
    max_layer_id = int(sorted(fiona.listlayers(path_to_file))[-1][-1])
    layer_id = min(layer_id, max_layer_id)
    layer_name = fiona.listlayers(path_to_file)[0][:-1] + str(layer_id)
    with fiona.open(path_to_file, "r", layer=layer_name) as country_file:
        for feature in filter(_in_study_area(study_area), country_file):
            new_feature = {}
            new_feature["properties"] = {}
            new_feature["properties"]["country_code"] = feature["properties"]["GID_0"]
            new_feature["properties"]["id"] = feature["properties"][f"GID_{layer_id}"]
            new_feature["properties"]["name"] = feature["properties"][f"NAME_{layer_id}"]
            new_feature["properties"]["type"] = (
                feature["properties"][f"ENGTYPE_{layer_id}"] if layer_id > 0
                else "country"
            )
            new_feature["properties"]["proper"] = True
            new_feature["geometry"] = _all_parts_in_study_area(feature, study_area)
            yield new_feature


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
        raise shapely.errors.TopologicalError(
            "Invalid study area geometry. "
            "Ensure that exclusion zones do not share a border with the study bounds."
        )
    else:
        return study_area


def _in_study_area(study_area):
    def _in_study_area(feature):
        unit = shapely.geometry.shape(feature["geometry"])
        if study_area.contains(unit) or study_area.intersects(unit):
            return True
        else:
            print("Removing {} as it is outside of study area.".format(_feature_name(feature)))
            return False
    return _in_study_area


def _all_parts_in_study_area(feature, study_area):
    unit = _to_multi_polygon(feature["geometry"])
    if not study_area.contains(unit):
        print("Removing parts of {} outside of study area.".format(_feature_name(feature)))
        new_unit = shapely.geometry.MultiPolygon([polygon for polygon in unit.geoms
                                                  if study_area.contains(polygon)])
        unit = new_unit
    return shapely.geometry.mapping(unit)


def _feature_name(feature):
    # brute force way of finding name
    # the problem is that the name of the feature depends on the layer
    for property_name in ["NAME_7", "NAME_6", "NAME_5", "NAME_4", "NAME_3", "NAME_2", "NAME_1",
                          "NAME_0"]:
        try:
            name = feature["properties"][property_name]
            break
        except KeyError:
            pass # nothing to do here
    return name


def _reproject(feature, src_crs, dst_crs):
    feature["geometry"] = fiona.transform.transform_geom(
        src_crs=src_crs,
        dst_crs=dst_crs,
        geom=feature["geometry"]
    )
    return feature


def _to_multi_polygon(geometry):
    if isinstance(geometry, dict):
        geometry = shapely.geometry.shape(geometry)
    if isinstance(geometry, shapely.geometry.polygon.Polygon):
        return shapely.geometry.MultiPolygon(polygons=[geometry])
    else:
        return geometry


def _test_id_uniqueness(path_to_file):
    for layer_name in fiona.listlayers(path_to_file):
        assert not gpd.read_file(path_to_file, layer=layer_name).id.duplicated().any()


if __name__ == "__main__":
    retrieve_administrative_borders()
