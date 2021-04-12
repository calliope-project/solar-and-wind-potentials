"""Determines capacity factors for each eligibility category on a map."""
import numpy as np
import rasterio

from src.technical_eligibility import Eligibility
def determine_capacityfactor(path_to_eligibility_categories, path_to_cfs, availability,
                             path_to_output_pv_prio, path_to_output_wind_prio):
    """Determines capacity factors for each eligibility category on a map."""
    with rasterio.open(path_to_eligibility_categories) as src:
        eligibility_categories = src.read(1)
    with rasterio.open(path_to_cfs[0]) as src:
        meta = src.meta
        rooftop_pv_cf = src.read(1)
    with rasterio.open(path_to_cfs[1]) as src:
        open_field_pv_cf = src.read(1)
    with rasterio.open(path_to_cfs[2]) as src:
        wind_onshore_cf = src.read(1)
        valid = wind_onshore_cf != meta["nodata"]
        wind_onshore_cf[valid] = wind_onshore_cf[valid] * availability["wind-onshore"]
    with rasterio.open(path_to_cfs[3]) as src:
        wind_offshore_cf = src.read(1)
        valid = wind_offshore_cf != meta["nodata"]
        wind_offshore_cf[valid] = wind_offshore_cf[valid] * availability["wind-offshore"]
    capacityfactor_pv_prio = _determine_capacityfactor(
        eligibility_category=eligibility_categories,
        rooftop_pv_cf=rooftop_pv_cf,
        open_field_pv_cf=open_field_pv_cf,
        wind_onshore_cf=wind_onshore_cf,
        wind_offshore_cf=wind_offshore_cf,
        pv_prio=True,
        nodata=meta["nodata"]
    )
    capacityfactor_wind_prio = _determine_capacityfactor(
        eligibility_category=eligibility_categories,
        rooftop_pv_cf=rooftop_pv_cf,
        open_field_pv_cf=open_field_pv_cf,
        wind_onshore_cf=wind_onshore_cf,
        wind_offshore_cf=wind_offshore_cf,
        pv_prio=False,
        nodata=meta["nodata"]
    )
    _write_to_file(path_to_output_pv_prio, capacityfactor_pv_prio, meta)
    _write_to_file(path_to_output_wind_prio, capacityfactor_wind_prio, meta)


def _determine_capacityfactor(eligibility_category, rooftop_pv_cf, open_field_pv_cf,
                              wind_onshore_cf, wind_offshore_cf, pv_prio, nodata):
    cf_map = np.zeros_like(rooftop_pv_cf)
    for eligibility in Eligibility:
        cf = _capacity_factor(
            eligibility=eligibility,
            pv_prio=pv_prio,
            rooftop_pv_cf=rooftop_pv_cf,
            open_field_pv_cf=open_field_pv_cf,
            wind_onshore_cf=wind_onshore_cf,
            wind_offshore_cf=wind_offshore_cf,
            nodata=nodata
        )
        mask = eligibility_category == eligibility
        cf_map[mask] = cf[mask]
    return cf_map


def _capacity_factor(eligibility, pv_prio, rooftop_pv_cf, open_field_pv_cf,
                     wind_onshore_cf, wind_offshore_cf, nodata):
    return {
        Eligibility.NOT_ELIGIBLE: np.ones_like(rooftop_pv_cf) * nodata,
        Eligibility.ROOFTOP_PV: rooftop_pv_cf,
        Eligibility.ONSHORE_WIND_AND_PV: open_field_pv_cf if pv_prio else wind_onshore_cf,
        Eligibility.ONSHORE_WIND: wind_onshore_cf,
        Eligibility.OFFSHORE_WIND: wind_offshore_cf
    }[eligibility]


def _write_to_file(path_to_file, data, meta):
    with rasterio.open(path_to_file, 'w', **meta) as new_geotiff:
        new_geotiff.write(data, 1)


if __name__ == "__main__":
    determine_capacityfactor(
        path_to_eligibility_categories=snakemake.input.eligibility_categories,
        path_to_cfs=snakemake.input.capacity_factors,
        availability=snakemake.params.availability,
        path_to_output_pv_prio=snakemake.output.pv,
        path_to_output_wind_prio=snakemake.output.wind,
    )
