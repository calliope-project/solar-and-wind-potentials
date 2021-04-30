"""This module estimates the shape of protected areas, for which only centroids are known.

This procedure is applied by the provider of the database, UNEP-WCMC, as well. See:
https://www.protectedplanet.net/c/calculating-protected-area-coverage
or the manual of the database for further information.
"""
import math

import geopandas as gpd
import pycountry

from renewablepotentialslib.geo import EPSG_3035_PROJ4


def estimate_shapes(path_to_input, scope_config, path_to_output):
    """Estimates the shap of protected areas for which only centroids are known."""
    points = gpd.read_file(path_to_input)
    points_in_scope = filter_points(points, scope_config)
    original_crs = points_in_scope.crs
    # convert points to circles
    points_in_scope = points_in_scope.to_crs(EPSG_3035_PROJ4)
    points_in_scope.geometry = [rec[1].geometry.buffer(radius_meter(rec[1]["REP_AREA"]))
                                for rec in points_in_scope.iterrows()]
    test_area_size(points_in_scope)
    points_in_scope.to_crs(original_crs).to_file(path_to_output, driver="GeoJSON")


def filter_points(points, scope_config):
    x_min, x_max, y_min, y_max = [scope_config["bounds"][z]
                                  for z in ["x_min", "x_max", "y_min", "y_max"]]
    countries = [pycountry.countries.lookup(country).alpha_3
                 for country in scope_config["countries"]]
    return points.cx[x_min:x_max, y_min:y_max].loc[
        (points.ISO3.isin(countries)) &
        (points.REP_AREA > 0)
    ].copy()


def radius_meter(area_squarekilometer):
    area_squaremeter = area_squarekilometer * 1e6
    return math.sqrt(area_squaremeter / math.pi)


def test_area_size(points):
    area_size_calculated = points.area.sum() / 1e6
    area_size_reported = points.REP_AREA.sum()
    assert abs(area_size_calculated - area_size_reported) < (area_size_reported / 100)


if __name__ == "__main__":
    estimate_shapes(
        path_to_input=snakemake.input.protected_areas,
        scope_config=snakemake.params.scope,
        path_to_output=snakemake.output[0]
    )
