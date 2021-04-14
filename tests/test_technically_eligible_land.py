import pytest
import numpy as np

from scripts.technical_eligibility import Eligibility, _determine_eligibility, GlobCover


@pytest.fixture
def config():
    return {
        "parameters": {
            "max-slope": {
                "pv": 3,
                "wind": 20
            },
            "max-depth-offshore": -50,
            "max-building-share": 0.1,
            "max-urban-green-share": 0.1
        }
    }


@pytest.mark.parametrize(
    "land_cover,slope,bathymetry,building_share,urban_green_share,expected", [
        (GlobCover.RAINFED_CROPLANDS, 0, 0, 0, 0, Eligibility.ONSHORE_WIND_AND_PV),
        (GlobCover.RAINFED_CROPLANDS, 4, 0, 0, 0, Eligibility.ONSHORE_WIND),
        (GlobCover.RAINFED_CROPLANDS, 21, 0, 0, 0, Eligibility.NOT_ELIGIBLE),
        (GlobCover.RAINFED_CROPLANDS, 0, 0, 0.11, 0, Eligibility.ROOFTOP_PV),
        (GlobCover.RAINFED_CROPLANDS, 0, 0, 0, 0.11, Eligibility.ROOFTOP_PV),
        (GlobCover.RAINFED_CROPLANDS, 0, 0, 0, 0, Eligibility.ONSHORE_WIND_AND_PV),
        (GlobCover.MOSAIC_FOREST, 0, 0, 0, 0, Eligibility.ONSHORE_WIND),
        (GlobCover.MOSAIC_FOREST, 21, 0, 0, 0, Eligibility.NOT_ELIGIBLE),
        (GlobCover.MOSAIC_FOREST, 0, 0, 0.11, 0, Eligibility.ROOFTOP_PV),
        (GlobCover.MOSAIC_FOREST, 0, 0, 0, 0.11, Eligibility.ROOFTOP_PV),
        (GlobCover.MOSAIC_FOREST, 0, 0, 0, 0, Eligibility.ONSHORE_WIND),
        (GlobCover.MOSAIC_GRASSLAND, 0, 0, 0, 0, Eligibility.ONSHORE_WIND_AND_PV),
        (GlobCover.MOSAIC_GRASSLAND, 4, 0, 0, 0, Eligibility.ONSHORE_WIND),
        (GlobCover.WATER_BODIES, 0, 0, 0, 0, Eligibility.OFFSHORE_WIND),
        (GlobCover.WATER_BODIES, 0, -51, 0, 0, Eligibility.NOT_ELIGIBLE),
        (GlobCover.WATER_BODIES, 0, 0, 0, 0, Eligibility.OFFSHORE_WIND),
        (GlobCover.ARTIFICAL_SURFACES_AND_URBAN_AREAS, 0, 0, 0, 0, Eligibility.NOT_ELIGIBLE)

    ]
)
def test_eligibility(land_cover, slope, bathymetry, building_share, urban_green_share,
                     expected, config):
    result = _determine_eligibility(
        land_cover=np.array([land_cover]),
        slope=np.array([slope]),
        bathymetry=np.array([bathymetry]),
        building_share=np.array([building_share]),
        urban_green_share=np.array([urban_green_share]),
        config=config
    )
    assert Eligibility(result[0]) == expected
