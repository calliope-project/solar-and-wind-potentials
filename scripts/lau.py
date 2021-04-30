"""Preprocessing of raw LAU2 data to bring it into normalised form."""
import geopandas as gpd
import pandas as pd

from renewablepotentialslib.geo.shape_utils import to_multi_polygon

OUTPUT_DRIVER = "GeoJSON"
KOSOVO_MUNICIPALITIES = [f"RS{x:02d}" for x in range(1, 38)]


def merge_lau(path_to_shapes, path_to_attributes, path_to_output):
    """Merge LAU shapes with attributes."""
    shapes = gpd.read_file(path_to_shapes)
    shapes.geometry = shapes.geometry.map(to_multi_polygon)
    attributes = gpd.read_file(path_to_attributes)
    attributes = pd.DataFrame(attributes) # to be able to remove the geo information
    del attributes["geometry"]
    all_shapes = shapes.merge(attributes, on="COMM_ID", how="left")
    all_shapes_no_kosovo = _remove_kosovo(all_shapes)
    all_shapes_no_kosovo.to_file(path_to_output, driver=OUTPUT_DRIVER)


def _remove_kosovo(shapes):
    """Identify and remove municipalities in Kosovo.

    Those municipalities must be removed as we do not have load data and pycountry
    cannot handle them at the moment (as of 2018, Kosovo does not have a standardised
    country code).
    """
    return shapes.set_index("COMM_ID").drop(KOSOVO_MUNICIPALITIES).reset_index()


if __name__ == "__main__":
    merge_lau(
        path_to_shapes=snakemake.input.shapes,
        path_to_attributes=snakemake.input.attributes,
        path_to_output=snakemake.output[0]
    )
