import glob
import os
import geopandas as gpd

from renewablepotentialslib.geo.shape_utils import (
    to_multi_polygon,
    drop_countries,
    drop_geoms_completely_outside_study_area,
    estimate_polygons_from_points
)

OUTPUT_DRIVER = "GeoJSON"


def wdpa(path_to_shape_dir, path_to_output, shapes_to_include, scope):
    all_shapes = gpd.GeoSeries()
    list_of_shapes_to_include = shapes_to_include.split(",")
    if not isinstance(list_of_shapes_to_include, list):
        list_of_shapes_to_include = [list_of_shapes_to_include]
    for geodata_type in list_of_shapes_to_include:
        shapes = _get_shapes(path_to_shape_dir, geodata_type, scope)
        all_shapes = all_shapes.append(shapes.geometry)
    all_shapes.to_file(path_to_output, driver=OUTPUT_DRIVER)


def _get_shapes(path_to_shapes, geodata_type, scope):
    shapefiles = glob.glob(os.path.join(path_to_shapes, "**", f"*-{geodata_type}.shp"))
    shapes = gpd.GeoDataFrame()
    for shapefile in shapefiles:
        print(f"Reading shape data from {shapefile}")
        shapes = shapes.append(gpd.read_file(shapefile))

    print(f"Number of shapes pre-filter = {len(shapes)}")
    shapes = _filter_shapes(shapes, scope)
    print(f"Number of shapes post-filter = {len(shapes)}")

    if geodata_type == "points":
        shapes = estimate_polygons_from_points(shapes, "REP_AREA")

    shapes.geometry = shapes.geometry.map(to_multi_polygon)

    return shapes


def _filter_shapes(shapes, scope):
    shapes = shapes[
        (shapes.REP_AREA.astype(float) > 0) &
        # The filter is in accordance to the way UNEP-WCMC calculates statistics:
        # https://www.protectedplanet.net/c/calculating-protected-area-coverage
        (shapes.STATUS.isin(['Designated', 'Inscribed', 'Established'])) &
        (shapes.DESIG_ENG != 'UNESCO-MAB Biosphere Reserve')
    ]
    shapes = drop_countries(shapes.rename(columns={"ISO3": "country_code"}), scope, print_dropped=False)
    shapes = drop_geoms_completely_outside_study_area(shapes, scope, print_dropped=False)

    return shapes


if __name__ == "__main__":
    wdpa(
        path_to_shape_dir=snakemake.input.shapes,
        shapes_to_include=snakemake.params.shapes_to_include,
        scope=snakemake.params.scope_config,
        path_to_output=snakemake.output[0]
    )
