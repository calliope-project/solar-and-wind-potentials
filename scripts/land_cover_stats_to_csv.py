"""Script to extract attributes of features and write to csv."""
import geopandas as gpd
from rasterstats import zonal_stats
import pandas as pd


def land_cover_stats_to_csv(path_to_units, path_to_land_cover, attributes, path_to_results_csv):
    """Extract attributes of features and write to csv."""
    #ipdb.set_trace()
    units = gpd.read_file(path_to_units).set_index('id')

    cmap = {int(attribute.split("_")[1]): attribute for attribute in attributes}
    stats = zonal_stats(path_to_units, path_to_land_cover, categorical=True, category_map=cmap)
    stat_df = pd.DataFrame(data=stats, index=units.index)
    stat_df.to_csv(path_to_results_csv)


if __name__ == "__main__":
    land_cover_stats_to_csv(
        path_to_units=snakemake.input.units,
        path_to_land_cover=snakemake.input.land_cover,
        attributes=snakemake.params.attributes,
        path_to_results_csv=snakemake.output[0]
    )
