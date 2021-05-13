import geopandas as gpd

from renewablepotentialslib.geo.shape_utils import (
    buffer_if_necessary,
    to_multi_polygon,
    drop_countries,
    drop_geoms_completely_outside_study_area,
    drop_parts_of_geoms_completely_outside_study_area,
    update_features
)

OUTPUT_DRIVER = 'GPKG'
SCHEMA = {
    "properties": {
        "country_code": "str", # id of the country to which the unit belongs
        "id": "str", # a unique id of this unit
        "name": "str", # the name of the unit, not necessarily unique
        "type": "str", # the type of the unit
        "proper": "bool" # flag indicating proper administrative unit (not the case for water bodies e.g.)
    },
    "geometry": "MultiPolygon"
}


def normalise_admin_borders(crs, scope_config, path_to_output, **shape_dirs):
    """Normalises raw administrative boundary data and places it in one, layered geodatapackage."""

    for _src, _path in shape_dirs.items():
        gdf = gpd.read_file(_path)
        gdf = gdf.to_crs(crs)
        gdf.geometry = gdf.geometry.map(buffer_if_necessary).map(to_multi_polygon)
        gdf = update_features(gdf, _src)
        gdf = drop_countries(gdf, scope_config)
        gdf = drop_geoms_completely_outside_study_area(gdf, scope_config)
        gdf = drop_parts_of_geoms_completely_outside_study_area(gdf, scope_config)

        assert gdf.id.duplicated().sum() == 0

        allowed_cols = list(SCHEMA["properties"].keys()) + ['geometry']

        for lvl in gdf.level.unique():
            gdf.loc[gdf.level == lvl, allowed_cols].to_file(
                path_to_output, schema=SCHEMA, layer=lvl, driver=OUTPUT_DRIVER
            )


if __name__ == "__main__":
    shape_dirs = {
        source.replace("SHAPEINPUTS_", ""): source_dir
        for source, source_dir in snakemake.input.items() if source.startswith("SHAPEINPUTS")
    }
    normalise_admin_borders(
        crs=snakemake.params.crs,
        scope_config=snakemake.params.scope,
        path_to_output=snakemake.output[0],
        **shape_dirs
    )
