"""Preprocessing of raw LAU2 data to bring it into normalised form."""
import click
import fiona
import fiona.transform
import shapely.geometry
import geopandas as gpd
import pandas as pd
import pycountry

from gadm import SCHEMA, _test_id_uniqueness, _study_area
from nuts import _to_multi_polygon
from conversion import eu_country_code_to_iso3
from utils import Config

OUTPUT_DRIVER = "GeoJSON"
KOSOVO_MUNICIPALITIES = [f"RS{x:02d}" for x in range(1, 38)]


@click.group()
def lau():
    pass


@lau.command()
@click.argument("path_to_shapes")
@click.argument("path_to_attributes")
@click.argument("path_to_output")
def merge(path_to_shapes, path_to_attributes, path_to_output):
    """Merge LAU shapes with attributes."""
    shapes = gpd.read_file(path_to_shapes)
    shapes.geometry = shapes.geometry.map(_to_multi_polygon)
    attributes = gpd.read_file(path_to_attributes)
    attributes = pd.DataFrame(attributes) # to be able to remove the geo information
    del attributes["geometry"]
    shapes.merge(attributes, on="COMM_ID", how="left").to_file(path_to_output, driver=OUTPUT_DRIVER)


@lau.command()
@click.argument("path_to_shapes")
@click.argument("path_to_output")
def identify(path_to_shapes, path_to_output):
    """Identify and remove municipalities in Kosovo.

    Those municipalities must be removed as we do not have load data and pycountry
    cannot handle them at the moment (as of 2018, Kosovo does not have a standardised
    country code).
    """
    gpd.read_file(path_to_shapes).set_index("COMM_ID").drop(KOSOVO_MUNICIPALITIES).reset_index().to_file(
        path_to_output,
        driver=OUTPUT_DRIVER
    )


@lau.command()
@click.argument("path_to_lau")
@click.argument("path_to_degurba")
@click.argument("path_to_output")
def degurba(path_to_lau, path_to_degurba, path_to_output):
    """Merge LAU2 units with DEGURBA data."""
    lau2 = gpd.read_file(path_to_lau)
    degurba = gpd.read_file(path_to_degurba)
    degurba_codes = lau2.merge(degurba, how="left", on=["CNTR_CODE", "NSI_CODE"])
    degurba_codes = pd.DataFrame(degurba_codes) # to be able to remove the geo information
    degurba_codes.rename(columns={
        "COMM_ID": "id",
        "DGURBA_CLA": "urbanisation_class"
    })[["id", "urbanisation_class"]].to_csv(path_to_output, header=True, index=False)


@lau.command()
@click.argument("path_to_lau")
@click.argument("path_to_output")
@click.argument("config", type=Config())
def normalise(path_to_lau, path_to_output, config):
    """Normalises raw LAU2 data."""
    with fiona.open(path_to_lau, "r") as lau_file, fiona.open(path_to_output,
                                                              "w",
                                                              crs=config["crs"],
                                                              schema=SCHEMA,
                                                              driver=OUTPUT_DRIVER) as result_file:
        result_file.writerecords(_layer_features(lau_file, config))
    _test_id_uniqueness(path_to_output)


def _layer_features(lau_file, config):
    for feature in filter(lambda feature: _in_study_area(config, feature), lau_file):
        new_feature = {}
        new_feature["properties"] = {}
        new_feature["properties"]["country_code"] = eu_country_code_to_iso3(feature["properties"]["COMM_ID"][:2])
        new_feature["properties"]["id"] = feature["properties"]["COMM_ID"]
        new_feature["properties"]["name"] = feature["properties"]["NAME_LATN"]
        new_feature["properties"]["type"] = "commune"
        new_feature["properties"]["proper"] = True if feature["properties"]["TRUE_COMM_"] == "T" else False
        new_feature["geometry"] = _all_parts_in_study_area_and_crs(feature, lau_file.crs, config)
        yield new_feature


def _all_parts_in_study_area_and_crs(feature, src_crs, config):
    study_area = _study_area(config)
    unit = _to_multi_polygon(feature["geometry"])
    if not study_area.contains(unit):
        print("Removing parts of {} outside of study area.".format(feature["properties"]["COMM_ID"]))
        new_unit = shapely.geometry.MultiPolygon([polygon for polygon in unit.geoms
                                                  if study_area.contains(polygon)])
        unit = new_unit
    geometry = shapely.geometry.mapping(unit)
    return fiona.transform.transform_geom(
        src_crs=src_crs,
        dst_crs=config["crs"],
        geom=geometry
    )


def _in_study_area(config, feature):
    study_area = _study_area(config)
    countries = [pycountry.countries.lookup(country) for country in config["scope"]["countries"]]
    unit = shapely.geometry.shape(feature["geometry"])
    country = pycountry.countries.lookup(eu_country_code_to_iso3(feature["properties"]["COMM_ID"][:2]))
    if (country in countries) and (study_area.contains(unit) or study_area.intersects(unit)):
        return True
    else:
        print("Removing {} as it is outside of study area.".format(feature["properties"]["COMM_ID"]))
        return False


if __name__ == "__main__":
    lau()
