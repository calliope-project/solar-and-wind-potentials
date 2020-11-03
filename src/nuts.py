"""Preprocessing of raw NUTS data to bring it into normalised form."""
import click
import geopandas as gpd
import pycountry

import administrative_border_utils
from conversion import eu_country_code_to_iso3
from utils import Config, buffer_if_necessary

OUTPUT_DRIVER = "GPKG"
LAYER_NAME = "nuts{}"


@click.command()
@click.argument("path_to_nuts")
@click.argument("path_to_output")
@click.argument("config", type=Config())
def normalise_nuts(path_to_nuts, path_to_output, config):
    """Normalises raw NUTS data.

    Raw data contains all NUTS layers in one layer of one shapefile. The output
    of this function corresponds to the form the data is used in this analysis,
    where each geographical layer is stored in one layer of a GeoPackage.
    """
    gdf = gpd.read_file(path_to_nuts)
    gdf = gdf.to_crs(config["crs"])
    gdf.geometry = gdf.geometry.buffer(0).map(administrative_border_utils._to_multi_polygon)
    gdf = _update_features(gdf)
    gdf = administrative_border_utils.drop_countries(gdf, config)
    gdf = administrative_border_utils.drop_geoms(gdf, config)

    administrative_border_utils.to_file(gdf, path_to_output, OUTPUT_DRIVER)


def _update_features(gdf):
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

    gdf["level"] = gdf.LEVL_CODE.apply(lambda i: LAYER_NAME.format(i))

    return gdf


if __name__ == "__main__":
    normalise_nuts()
