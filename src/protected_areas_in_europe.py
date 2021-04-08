import math

import glob
import os
import click
import geopandas as gpd

from src.utils import Config

from src.administrative_borders import _to_multi_polygon, _study_area, _drop_countries, _drop_geoms_completely_outside_study_area, _drop_parts_of_geoms_completely_outside_study_area

OUTPUT_DRIVER = "GeoJSON"

@click.command()
@click.argument("path_to_shape_dir")
@click.argument("path_to_output")
@click.argument("shapes_to_include")
@click.argument("config", type=Config())
def wdpa(path_to_shape_dir, path_to_output, shapes_to_include, config):
    all_shapes = gpd.GeoSeries()
    list_of_shapes_to_include = shapes_to_include.split(",")
    if not isinstance(list_of_shapes_to_include, list):
        list_of_shapes_to_include = [list_of_shapes_to_include]
    for geodata_type in list_of_shapes_to_include:
        shapes = _get_shapes(path_to_shape_dir, geodata_type, config)
        all_shapes = all_shapes.append(shapes.geometry)
    all_shapes.to_file(path_to_output, driver=OUTPUT_DRIVER)


def _get_shapes(path_to_shapes, geodata_type, config):
    shapefiles = glob.glob(os.path.join(path_to_shapes, "**", f"*-{geodata_type}.shp"))
    shapes = gpd.GeoDataFrame()
    for shapefile in shapefiles:
        print(f"Reading shape data from {shapefile}")
        shapes = shapes.append(gpd.read_file(shapefile))

    print(f"Number of shapes pre-filter = {len(shapes)}")
    shapes = _filter_shapes(shapes, config)
    print(f"Number of shapes post-filter = {len(shapes)}")

    if geodata_type == "points":
        shapes = _estimate_polygons_from_points(shapes)

    shapes.geometry = shapes.geometry.map(_to_multi_polygon)

    return shapes


def _filter_shapes(shapes, config):
    shapes = shapes[
        (shapes.REP_AREA.astype(float) > 0) &
        # The filter is in accordance to the way UNEP-WCMC calculates statistics:
        # https://www.protectedplanet.net/c/calculating-protected-area-coverage
        (shapes.STATUS.isin(['Designated', 'Inscribed', 'Established'])) &
        (shapes.DESIG_ENG != 'UNESCO-MAB Biosphere Reserve')
    ]
    shapes = _drop_countries(shapes.rename(columns={"ISO3": "country_code"}), config, print_dropped=False)
    shapes = _drop_geoms_completely_outside_study_area(shapes, config, print_dropped=False)

    return shapes


def _estimate_polygons_from_points(points):
    """Estimates the shape of protected areas for which only centroids are known."""
    def __radius_meter(area_squarekilometer):
        area_squaremeter = area_squarekilometer * 1e6
        return math.sqrt(area_squaremeter / math.pi)

    def __test_area_size(points):
        area_size_calculated = points.area.sum() / 1e6
        area_size_reported = points.REP_AREA.sum()
        assert abs(area_size_calculated - area_size_reported) < (area_size_reported / 100)

    original_crs = points.crs
    # convert points to circles
    points_in_metres = points.to_crs("epsg:3035")
    points_in_metres.geometry = [
        point_data.geometry.buffer(__radius_meter(point_data["REP_AREA"]))
        for _point_idx, point_data in points_in_metres.iterrows()
    ]
    __test_area_size(points_in_metres)
    return points_in_metres.to_crs(original_crs)





if __name__ == "__main__":
    wdpa()
