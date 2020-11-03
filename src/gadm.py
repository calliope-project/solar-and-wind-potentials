"""Module to merge and preprocess GADM administrative borders."""
from itertools import chain

import click
import fiona
import fiona.transform
import geopandas as gpd
import shapely.geometry
import shapely.errors
from shapely.prepared import prep

import administrative_border_utils
from utils import Config

OUTPUT_DRIVER = "GPKG"
LAYER_NAME = "gadm{}"


@click.command()
@click.argument("path_to_gadm")
@click.argument("path_to_output")
@click.argument("config", type=Config())
def normalise_gadm(path_to_gadm, path_to_output, config):
    gdf = gpd.read_file(path_to_gadm)
    gdf.geometry = gdf.geometry.buffer(0).map(administrative_border_utils._to_multi_polygon)
    gdf = _update_features(gdf)
    gdf = administrative_border_utils.drop_countries(gdf, config)
    gdf = administrative_border_utils.drop_geoms(gdf, config)

    administrative_border_utils.to_file(gdf, path_to_output, OUTPUT_DRIVER)


def _update_features(gdf):
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
    gdf["level"] = gdf.level.apply(lambda i: LAYER_NAME.format(i))

    return gdf


if __name__ == "__main__":
    normalise_gadm()
