"""Module to Determine share of shared coast between eez and administrative units."""
from textwrap import dedent
from multiprocessing import Pool
from itertools import cycle

import numpy as np
import pandas as pd
import geopandas as gpd
from shapely.prepared import prep

from renewablepotentialslib.geo.shape_utils import buffer_if_necessary

DRIVER = "GeoJSON"


def allocate_eezs(path_to_units, path_to_eezs, path_to_output, threads):
    """Determine share of shared coast between eez and administrative units."""
    units = gpd.read_file(path_to_units)
    units.set_index("id", inplace=True)
    eezs = gpd.read_file(path_to_eezs)
    eezs.geometry = eezs.geometry.map(buffer_if_necessary)
    with Pool(int(threads)) as pool:
        share_of_coast_length = pool.map(
            _share_of_coast_length,
            zip((eez[1] for eez in eezs.iterrows()), cycle([units]))
        )
    share = pd.DataFrame(
        index=units.index,
        data=dict(zip(eezs["id"].values, share_of_coast_length))
    )
    assert (
        ((share.sum() > 0.99) & (share.sum() < 1.01)) |
        (share.sum() == 0.0)
    ).all(), share.sum()
    share.to_csv(path_to_output, header=True)


def _share_of_coast_length(args):
    # How to determine the length of the shared coast?
    # I intersect eez with the unit and determine the length of the resulting polygon.
    # This approach is fairly rough, but accurate enough for this analysis.
    eez = args[0]
    units = args[1]
    length_of_shared_coast = pd.Series(data=0.0, index=units.index, dtype=np.float32)
    prep_eez = prep(eez["geometry"]) # increase performance
    intersection_mask = ((units["country_code"].isin([eez["ISO_Ter1"], "EUR"])) &
                         (units["geometry"].map(lambda unit: prep_eez.intersects(unit))))
    if intersection_mask.sum() == 0:
        msg = dedent("""No shared coast found for {}.
        Ignoring eez with area {} km^2.""".format(
            eez["GeoName"],
            eez["Area_km2"]
        ))
        print(msg)
        share = length_of_shared_coast.copy()
    elif intersection_mask.sum() == 1:
        # performance improvement in cases where only one unit matches
        share = length_of_shared_coast.copy()
        share[intersection_mask] = 1
    else:
        length_of_shared_coast[intersection_mask] = units.loc[intersection_mask, "geometry"].map(
            lambda unit: eez["geometry"].intersection(unit).length
        )
        share = length_of_shared_coast / length_of_shared_coast.sum()
    return share


if __name__ == "__main__":
    allocate_eezs(
        path_to_units=snakemake.input.units,
        path_to_eezs=snakemake.input.eez,
        path_to_output=snakemake.output[0],
        threads=snakemake.threads
    )
