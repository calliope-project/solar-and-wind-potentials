"""Determines maximal electricity yield for renewables."""
from datetime import timedelta

import numpy as np
import rasterio

from renewablepotentialslib.geo.conversion import watt_to_watthours


def determine_electricity_yield(path_to_eligibility_categories, path_to_capacities_pv_prio,
                                path_to_capacities_wind_prio, path_to_cf_pv_prio,
                                path_to_cf_wind_prio,
                                path_to_pv_prio_result, path_to_wind_prio_result):
    """Determines maximal electricity yield for renewables."""
    with rasterio.open(path_to_capacities_pv_prio) as src:
        meta = src.meta
        capacity_pv_prio_mw = src.read(1)
    with rasterio.open(path_to_capacities_wind_prio) as src:
        capacity_wind_prio_mw = src.read(1)
    with rasterio.open(path_to_eligibility_categories) as src:
        eligibility_categories = src.read(1)
    with rasterio.open(path_to_cf_pv_prio) as src:
        no_cf = src.nodata
        cf_pv_prio = src.read(1)
    with rasterio.open(path_to_cf_wind_prio) as src:
        cf_wind_prio = src.read(1)
    electricity_yield_pv_prio = _determine_electricity_yield(
        capacity_mw=capacity_pv_prio_mw,
        eligibility_category=eligibility_categories,
        cf=cf_pv_prio,
        data_mask=cf_pv_prio != no_cf
    )
    electricity_yield_wind_prio = _determine_electricity_yield(
        capacity_mw=capacity_wind_prio_mw,
        eligibility_category=eligibility_categories,
        cf=cf_wind_prio,
        data_mask=cf_wind_prio != no_cf
    )
    _write_to_file(path_to_pv_prio_result, electricity_yield_pv_prio, meta)
    _write_to_file(path_to_wind_prio_result, electricity_yield_wind_prio, meta)


def _determine_electricity_yield(capacity_mw, eligibility_category, cf, data_mask):
    electricity_yield_twh = np.zeros_like(capacity_mw)
    cf[data_mask] = watt_to_watthours(cf[data_mask], duration=timedelta(days=365))
    electricity_yield_twh[data_mask] = (capacity_mw * cf)[data_mask] / 1e6
    return electricity_yield_twh


def _write_to_file(path_to_file, electricity_yield, meta):
    with rasterio.open(path_to_file, 'w', **meta) as new_geotiff:
        new_geotiff.write(electricity_yield, 1)


if __name__ == "__main__":
    determine_electricity_yield(
        path_to_eligibility_categories=snakemake.input.eligibility_categories,
        path_to_capacities_pv_prio=snakemake.input.capacities_pv_prio,
        path_to_capacities_wind_prio=snakemake.input.capacities_wind_prio,
        path_to_cf_pv_prio=snakemake.input.cf_pv_prio,
        path_to_cf_wind_prio=snakemake.input.cf_wind_prio,
        path_to_pv_prio_result=snakemake.output.pv,
        path_to_wind_prio_result=snakemake.output.wind
    )
