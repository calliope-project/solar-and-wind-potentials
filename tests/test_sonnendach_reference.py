"""Test whether our estimations are close to the ones from sonnendach.ch"""
import os
from pathlib import Path

import pytest
import rasterio
import rasterio.mask
from rasterstats import zonal_stats
import fiona

from renewablepotentialslib.eligibility import Eligibility

ROOT_DIR = Path(os.path.abspath(__file__)).parent.parent
PATH_TO_CATEGORIES = ROOT_DIR / "build" / "technically-eligible-land.tif"
PATH_TO_AREAS = ROOT_DIR / "build" / "technically-eligible-area-km2.tif"
PATH_TO_ENERGY_YIELD = ROOT_DIR / "build" / "technically-eligible-electricity-yield-pv-prio-twh.tif"
PATH_TO_NUTS = ROOT_DIR / "build" / "administrative-borders.gpkg"
PATH_TO_SONNENDACH_AREA_ESTIMATE = ROOT_DIR / "data" / "automatic" / "sonnendach" /\
    "total-rooftop-area-km2.txt"
PATH_TO_SONNENDACH_YIELD_ESTIMATE = ROOT_DIR / "data" / "automatic" / "sonnendach" /\
    "total-yield-twh.txt"


@pytest.mark.skipif(not PATH_TO_AREAS.exists(), reason="Eligible area raster data not available.")
@pytest.mark.skipif(not PATH_TO_NUTS.exists(), reason="Switzerland shape not available.")
@pytest.mark.skipif(not PATH_TO_SONNENDACH_AREA_ESTIMATE.exists(), reason="Sonnendach area estimation not available.")
def test_switzerland_rooftop_area():
    with open(PATH_TO_SONNENDACH_AREA_ESTIMATE, "r") as f_sonnendach_estimate:
        sonnendach_estimate = float(f_sonnendach_estimate.readline())
    with fiona.open(PATH_TO_NUTS.as_posix(), "r", layer="nuts0") as shapefile:
        switzerland = [feature["geometry"] for feature in shapefile if feature["properties"]["country_code"] == "CHE"]
        assert len(switzerland) == 1
    with rasterio.open(PATH_TO_AREAS.as_posix()) as src:
        transform = src.transform
        areas = src.read(1)
    with rasterio.open(PATH_TO_CATEGORIES.as_posix()) as src:
        categories = src.read(1)
    areas[categories != Eligibility.ROOFTOP_PV] = 0
    zs = zonal_stats(switzerland, areas, affine=transform, stats="sum", nodata=-999)
    our_estimate = zs[0]["sum"]
    assert our_estimate == pytest.approx(sonnendach_estimate, 0.02) # 2% tolerance


@pytest.mark.skipif(not PATH_TO_ENERGY_YIELD.exists(), reason="Eligible energy yield raster data not available.")
@pytest.mark.skipif(not PATH_TO_NUTS.exists(), reason="Switzerland shape not available.")
@pytest.mark.skipif(
    not PATH_TO_SONNENDACH_YIELD_ESTIMATE.exists(),
    reason="Sonnendach yield estimation not available.")
def test_switzerland_energy_yield():
    with open(PATH_TO_SONNENDACH_YIELD_ESTIMATE, "r") as f_sonnendach_estimate:
        sonnendach_estimate = float(f_sonnendach_estimate.readline())
    with fiona.open(PATH_TO_NUTS.as_posix(), "r", layer="nuts0") as shapefile:
        switzerland = [feature["geometry"] for feature in shapefile if feature["properties"]["country_code"] == "CHE"]
        assert len(switzerland) == 1
    with rasterio.open(PATH_TO_ENERGY_YIELD.as_posix()) as src:
        transform = src.transform
        energy_yield = src.read(1)
    with rasterio.open(PATH_TO_CATEGORIES.as_posix()) as src:
        categories = src.read(1)
    energy_yield[categories != Eligibility.ROOFTOP_PV] = 0
    zs = zonal_stats(switzerland, energy_yield, affine=transform, stats="sum", nodata=-999)
    our_estimate = zs[0]["sum"]
    assert our_estimate <= sonnendach_estimate
    assert our_estimate == pytest.approx(sonnendach_estimate, 0.10) # 10% tolerance
