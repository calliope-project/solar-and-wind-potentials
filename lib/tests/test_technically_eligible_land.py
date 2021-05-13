import pytest
import numpy as np

from renewablepotentialslib.eligibility import eligibility_land_mask, Eligibility, GlobCover

@pytest.fixture
def config():
    return {
        "max-depth-offshore": -50,
        "max-building-share": 0.1,
        "max-urban-green-share": 0.1,
        "max-slope-pixel-fraction-threshold": 0.9
    }


@pytest.mark.parametrize(
    "land_cover,slope_pv,slope_wind,bathymetry,building_share,urban_green_share,expected", [
        (GlobCover.RAINFED_CROPLANDS, 1, 1, 0, 0, 0, Eligibility.ONSHORE_WIND_AND_PV),
        (GlobCover.RAINFED_CROPLANDS, 0, 1, 0, 0, 0, Eligibility.ONSHORE_WIND),
        (GlobCover.RAINFED_CROPLANDS, 0.8, 0.8, 0, 0, 0, Eligibility.NOT_ELIGIBLE),
        (GlobCover.RAINFED_CROPLANDS, 1, 1, 0, 0.11, 0, Eligibility.ROOFTOP_PV),
        (GlobCover.RAINFED_CROPLANDS, 1, 1, 0, 0, 0.11, Eligibility.ROOFTOP_PV),
        (GlobCover.RAINFED_CROPLANDS, 1, 1, 0, 0, 0, Eligibility.ONSHORE_WIND_AND_PV),
        (GlobCover.MOSAIC_FOREST, 1, 1, 0, 0, 0, Eligibility.ONSHORE_WIND),
        (GlobCover.MOSAIC_FOREST, 0.8, 0.8, 0, 0, 0, Eligibility.NOT_ELIGIBLE),
        (GlobCover.MOSAIC_FOREST, 1, 1, 0, 0.11, 0, Eligibility.ROOFTOP_PV),
        (GlobCover.MOSAIC_FOREST, 1, 1, 0, 0, 0.11, Eligibility.ROOFTOP_PV),
        (GlobCover.MOSAIC_FOREST, 1, 1, 0, 0, 0, Eligibility.ONSHORE_WIND),
        (GlobCover.MOSAIC_GRASSLAND, 1, 1, 0, 0, 0, Eligibility.ONSHORE_WIND_AND_PV),
        (GlobCover.MOSAIC_GRASSLAND, 0.8, 1, 0, 0, 0, Eligibility.ONSHORE_WIND),
        (GlobCover.WATER_BODIES, 1, 1, 0, 0, 0, Eligibility.OFFSHORE_WIND),
        (GlobCover.WATER_BODIES, 1, 1, -51, 0, 0, Eligibility.NOT_ELIGIBLE),
        (GlobCover.WATER_BODIES, 1, 1, 0, 0, 0, Eligibility.OFFSHORE_WIND),
        (GlobCover.ARTIFICAL_SURFACES_AND_URBAN_AREAS, 1, 1, 0, 0, 0, Eligibility.NOT_ELIGIBLE)

    ]
)
def test_eligibility(land_cover, slope_pv, slope_wind, bathymetry, building_share, urban_green_share,
                     expected, config):
    max_depth_offshore = config["max-depth-offshore"]
    max_building_share = config["max-building-share"]
    max_urban_green_share = config["max-urban-green-share"]
    slope_threshold = config["max-slope-pixel-fraction-threshold"]
    result = eligibility_land_mask(
        land_cover=np.array([land_cover]),
        slope_pv=np.array([slope_pv]),
        slope_wind=np.array([slope_wind]),
        bathymetry=np.array([bathymetry]),
        building_share=np.array([building_share]),
        urban_green_share=np.array([urban_green_share]),
        max_depth_offshore=max_depth_offshore,
        max_building_share=max_building_share,
        max_urban_green_share=max_urban_green_share,
        slope_threshold=slope_threshold
    )
    assert Eligibility(result[0]) == expected
