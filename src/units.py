import click
import pandas as pd
import geopandas as gpd
import pycountry

from utils import Config

DRIVER = "GeoJSON"


@click.command()
@click.argument("path_to_borders")
@click.argument("path_to_output")
@click.argument("layer_name")
@click.argument("config", type=Config())
def remix_units(path_to_borders, path_to_output, layer_name, config):
    """Remixes NUTS, LAU, and GADM data to form the units of the analysis."""
    source_layers = _read_source_layers(path_to_borders, config["layers"][layer_name])
    _validate_source_layers(source_layers)
    _validate_layer_config(config, layer_name)
    layer = _build_layer(config["layers"][layer_name], source_layers)
    _validate_layer(layer, layer_name, config["scope"]["countries"])
    if layer_name == "continental": # treat special case
        layer = _continental_layer(layer)
    _write_layer(layer, path_to_output)


def _read_source_layers(path_to_borders, layers):
    source_layers = {
        layer_name: gpd.read_file(path_to_borders, layer=layer_name)
        for layer_name in set(layers.values())
    }
    return source_layers


def _validate_source_layers(source_layers):
    crs = [layer.crs for layer in source_layers.values()]
    assert not crs or crs.count(crs[0]) == len(crs), "Source layers have different crs. They must match."


def _validate_layer_config(config, layer_name):
    country_scope = config["scope"]["countries"]
    layer = config["layers"][layer_name]
    assert all(country in layer.keys() for country in country_scope), ("Layer {} is not correctly "
                                                                       "defined.".format(layer_name))


def _build_layer(country_to_source_map, source_layers):
    crs = [layer.crs for layer in source_layers.values()][0]
    layer = pd.concat([
        source_layers[source_layer][source_layers[source_layer].country_code == _iso3(country)]
        for country, source_layer in country_to_source_map.items()
    ])
    assert isinstance(layer, pd.DataFrame)
    return gpd.GeoDataFrame(layer, crs=crs)


def _validate_layer(layer, layer_name, countries):
    assert all(_iso3(country) in layer.country_code.unique()
               for country in countries), (f"Countries are missing in layer {layer_name}.")


def _iso3(country_name):
    return pycountry.countries.lookup(country_name).alpha_3


def _continental_layer(layer):
    # special case all Europe
    layer = layer.dissolve(by=[1 for idx in layer.index])
    index = layer.index[0]
    layer.loc[index, "id"] = "EUR"
    layer.loc[index, "country_code"] = "EUR"
    layer.loc[index, "name"] = "Europe"
    layer.loc[index, "type"] = "continent"
    layer.loc[index, "proper"] = 1
    return layer


def _write_layer(gdf, path_to_file):
    gdf.to_file(
        path_to_file,
        driver=DRIVER
    )


if __name__ == "__main__":
    remix_units()
