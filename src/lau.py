"""Preprocessing of raw LAU2 data to bring it into normalised form."""
import click
import geopandas as gpd
import pandas as pd
import pycountry

import administrative_border_utils
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
    shapes.geometry = shapes.geometry.map(administrative_border_utils._to_multi_polygon)
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

    gdf = gpd.read_file(path_to_lau)
    gdf = gdf.to_crs(config["crs"])
    gdf.geometry = gdf.geometry.buffer(0).map(administrative_border_utils._to_multi_polygon)
    gdf = _update_features(gdf)
    gdf = administrative_border_utils.drop_countries(gdf, config)
    gdf = administrative_border_utils.drop_geoms(gdf, config)

    administrative_border_utils.to_file(gdf, path_to_output, OUTPUT_DRIVER)


def _update_features(gdf):
    gdf["CNTR_CODE"] = gdf.COMM_ID.str[:2].apply(eu_country_code_to_iso3)
    gdf = gdf.rename(
        columns={"COMM_ID": "id", "NAME_LATN": "name", "CNTR_CODE": "country_code"}
    )
    gdf["type"] = "commune"
    gdf["proper"] = gdf["TRUE_COMM_"] == "T"
    gdf["level"] = "lau2"

    return gdf


if __name__ == "__main__":
    lau()
