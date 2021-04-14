"""Create maps of time averaged capacitfy factors of renewables."""
import numpy as np
import rasterio
import xarray as xr

from timeseries import CAPACITY_FACTOR_VAR

DTYPE = np.float32
NODATA = -1


def averages_map(path_to_id_map, path_to_timeseries, path_to_output):
    """Create maps of time averaged capacitfy factors of renewables."""
    with rasterio.open(str(path_to_id_map), "r") as f_ids:
        ids = f_ids.read(1)
        meta = f_ids.meta
    averages = map_id_to_average_capacity_factor(ids, path_to_timeseries, meta["nodata"])
    meta["dtype"] = DTYPE
    meta["nodata"] = NODATA
    with rasterio.open(str(path_to_output), "w", **meta) as f_avg:
        f_avg.write(averages, 1)


def map_id_to_average_capacity_factor(ids, path_to_timeseries, nodata_id):
    average_capacity_factors = xr.open_dataset(str(path_to_timeseries)).mean("time")[CAPACITY_FACTOR_VAR].to_dataframe()
    average_capacity_factors.index = average_capacity_factors.index.astype(np.int32)
    average_capacity_factors = average_capacity_factors.to_dict()[CAPACITY_FACTOR_VAR]
    average_capacity_factors[nodata_id] = NODATA
    mapping_function = np.vectorize(
        lambda site_id: average_capacity_factors[site_id],
        otypes=[DTYPE]
    )
    return mapping_function(ids)


if __name__ == "__main__":
    averages_map(
        path_to_id_map=snakemake.input.id_map,
        path_to_timeseries=snakemake.input.timeseries,
        path_to_output=snakemake.output
    )
