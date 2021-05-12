import glob
import os
import geopandas as gpd
import pandas as pd

from renewablepotentialslib.shape_utils import (
    to_multi_polygon,
    drop_countries,
    estimate_polygons_from_points,
    study_area
)

OUTPUT_DRIVER = "GeoJSON"


def wdpa(path_to_shape_dir, path_to_output, shapefile_prefix, shapes_to_include, scope):
    all_shapes = pd.concat([
        _get_shapes(path_to_shape_dir, shapefile_prefix, geodata_type, scope).geometry
        for geodata_type in shapes_to_include
    ])

    all_shapes.to_file(path_to_output, driver=OUTPUT_DRIVER)


def _get_shapes(path_to_shape_dir, shapefile_prefix, geodata_type, scope):
    shapefile_zips = glob.glob(os.path.join(path_to_shape_dir, "*.zip"))
    shapefile_name = f"zip://{{}}!{shapefile_prefix}-{geodata_type}.shp"
    bbox = study_area(scope)
    shapes = pd.concat([
        gpd.read_file(shapefile_name.format(zipfile), bbox=bbox)
        for zipfile in shapefile_zips
    ])

    print(f"Number of WDPA shapes pre-filter = {len(shapes)}")
    shapes = _filter_shapes(shapes, scope)
    print(f"Number of WDPA shapes post-filter = {len(shapes)}")

    if geodata_type == "points":
        shapes = estimate_polygons_from_points(shapes, "REP_AREA")
    print(f"Processed points to polys")
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

    return shapes


if __name__ == "__main__":
    wdpa(
        path_to_shape_dir=snakemake.input.shape_dir,
        shapefile_prefix=snakemake.params.shapefile_prefix,
        shapes_to_include=snakemake.params.shapes_to_include,
        scope=snakemake.params.scope_config,
        path_to_output=snakemake.output[0]
    )
