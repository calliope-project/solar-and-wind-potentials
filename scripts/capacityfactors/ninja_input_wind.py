"""Create wind simulation input for renewables.ninja."""
import click
import pandas as pd
import geopandas as gpd

from src.utils import Config
from src.capacityfactors import point_raster_on_shapes


@click.command()
@click.argument("path_to_shapes_of_land_surface")
@click.argument("path_to_shapes_of_water_surface")
@click.argument("path_to_onshore_output")
@click.argument("path_to_offshore_output")
@click.argument("config", type=Config())
def wind(path_to_shapes_of_land_surface, path_to_shapes_of_water_surface, path_to_onshore_output,
         path_to_offshore_output, config):
    """Create wind on- and offshore simulation input for renewables.ninja."""
    write_parameters(
        bounds=config["scope"]["bounds"],
        resolution=config["parameters"]["ninja"]["resolution-grid"],
        path_to_shapes=path_to_shapes_of_land_surface,
        hub_height=config["parameters"]["ninja"]["hub-height"]["onshore"],
        turbine=config["parameters"]["ninja"]["turbine"]["onshore"],
        path_to_output=path_to_onshore_output
    )
    write_parameters(
        bounds=config["scope"]["bounds"],
        resolution=config["parameters"]["ninja"]["resolution-grid"],
        path_to_shapes=path_to_shapes_of_water_surface,
        hub_height=config["parameters"]["ninja"]["hub-height"]["offshore"],
        turbine=config["parameters"]["ninja"]["turbine"]["offshore"],
        path_to_output=path_to_offshore_output
    )


def write_parameters(bounds, resolution, path_to_shapes, hub_height, turbine, path_to_output):
    points = point_raster_on_shapes(
        bounds_wgs84=bounds,
        shapes=gpd.read_file(path_to_shapes),
        resolution_km2=resolution
    )
    parameters = pd.DataFrame(
        data={
            "lat": [point.y for point in points.geometry],
            "long": [point.x for point in points.geometry],
            "hub_height": hub_height,
            "turbine": turbine
        }
    )
    parameters["sim_id"] = parameters.index
    parameters["site_id"] = parameters.index
    parameters["weight"] = 1
    parameters[["sim_id", "weight", "site_id", "lat", "long", "hub_height", "turbine"]].to_csv(
        path_to_output,
        header=True,
        index=False
    )


if __name__ == "__main__":
    wind()
