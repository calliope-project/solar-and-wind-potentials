import geopandas as gpd

from renewablepotentialslib.shape_utils import (
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


def normalise_admin_borders(path_to_nuts, path_to_gadm, path_to_lau, crs, scope_config, path_to_output):
    """Normalises raw administrative boundary data and places it in one, layered geodatapackage."""

    for _src, _path in {
        'nuts': path_to_nuts, 'gadm': path_to_gadm, 'lau': path_to_lau
    }.items():
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
    normalise_admin_borders(
        path_to_nuts=snakemake.input.nuts_geojson,
        path_to_gadm=snakemake.input.gadm_gpkg,
        path_to_lau=snakemake.input.lau_gpkg,
        crs=snakemake.params.crs,
        scope_config=snakemake.params.scope,
        path_to_output=snakemake.output[0]
    )
