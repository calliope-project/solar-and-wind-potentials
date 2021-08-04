"""Determine total Swiss building footprint of ESM data."""
import rasterio
import fiona
from rasterstats import zonal_stats
import pandas as pd
import geopandas as gpd

from renewablepotentialslib.eligibility import Eligibility
from renewablepotentialslib.geo.conversion import area_in_squaremeters


def swiss_building_footprint(path_to_building_footprint, path_to_eligibility,
                             path_to_countries, path_to_output):
    with rasterio.open(path_to_eligibility, "r") as f_eligibility:
        eligibility = f_eligibility.read(1)
    with rasterio.open(path_to_building_footprint, "r") as f_building_share:
        building_share = f_building_share.read(1)
        transform = f_building_share.transform
    building_share[eligibility != Eligibility.ROOFTOP_PV] = 0

    with fiona.open(path_to_countries, "r", layer="nuts0") as src:
        zs = zonal_stats(
            vectors=src,
            raster=building_share,
            affine=transform,
            stats="mean",
            nodata=-999
        )
        building_share = pd.Series(
            index=[feat["properties"]["id"] for feat in src],
            data=[stat["mean"] for stat in zs]
        )
    building_footprint_km2 = (
        area_in_squaremeters(gpd.read_file(path_to_countries).set_index("id"))
        .div(1e6)
        .mul(building_share)
    )
    swiss_building_footprint = building_footprint_km2.loc["CHE"]
    with open(path_to_output, "w") as f_out:
        f_out.write(f"{swiss_building_footprint}")


if __name__ == "__main__":
    swiss_building_footprint(
        path_to_building_footprint=snakemake.input.building_footprints,
        path_to_eligibility=snakemake.input.eligibility,
        path_to_countries=snakemake.input.countries,
        path_to_output=snakemake.output[0]
    )
