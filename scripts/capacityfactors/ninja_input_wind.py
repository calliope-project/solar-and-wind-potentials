"""Create wind simulation input for renewables.ninja."""
import pandas as pd
import geopandas as gpd

from renewablepotentialslib.geo.shape_utils import point_raster_on_shapes


def wind(path_to_shapes_of_land_surface, path_to_shapes_of_water_surface, bounds, ninja,
         path_to_onshore_output, path_to_offshore_output):
    """Create wind on- and offshore simulation input for renewables.ninja."""
    write_parameters(
        bounds=bounds,
        resolution=ninja["resolution-grid"],
        path_to_shapes=path_to_shapes_of_land_surface,
        hub_height=ninja["hub-height"]["onshore"],
        turbine=ninja["turbine"]["onshore"],
        path_to_output=path_to_onshore_output
    )
    write_parameters(
        bounds=bounds,
        resolution=ninja["resolution-grid"],
        path_to_shapes=path_to_shapes_of_water_surface,
        hub_height=ninja["hub-height"]["offshore"],
        turbine=ninja["turbine"]["offshore"],
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
    wind(
        path_to_shapes_of_land_surface=snakemake.input.units,
        path_to_shapes_of_water_surface=snakemake.input.eez,
        bounds=snakemake.params.bounds,
        ninja=snakemake.params.ninja,
        path_to_onshore_output=snakemake.output.points_onshore,
        path_to_offshore_output=snakemake.output.points_offhore
    )
