import io

import pytest
import pandas as pd
from pandas.testing import assert_frame_equal

from scripts.capacityfactors.ninja_input_pv import area_to_capacity

ROOF_MODEL = """orientation,average_tilt,share_of_roof_areas
E, 18.155579, 0.049090
E, 25.863758, 0.039782
E, 32.876361, 0.036700
E, 43.523447, 0.040453
N, 17.312256, 0.048285
N, 24.879743, 0.041521
N, 32.361540, 0.046410
N, 43.655379, 0.045527
S, 18.063436, 0.055544
S, 25.412273, 0.036332
S, 32.368793, 0.047489
S, 43.059819, 0.042767
W, 18.107352, 0.051856
W, 25.376952, 0.034674
W, 32.340545, 0.041847
W, 43.504116, 0.039763
flat, 0.000000, 0.301960
"""


@pytest.fixture()
def roof_model_area_based():
    model = pd.read_csv(io.StringIO(ROOF_MODEL)).set_index(["orientation", "average_tilt"])
    assert model.sum().sum() == pytest.approx(1.0)
    return model


def test_capacity_based_roof_model_sums_to_one(roof_model_area_based):
    roof_model_capacity_based = area_to_capacity(
        roof_model_area_based,
        power_density_flat=1,
        power_density_tilted=2
    )
    assert roof_model_capacity_based["share_of_roof_areas"].sum() == pytest.approx(1.0)


def test_capacity_based_roof_model_equals_area_based_for_equal_power_density(roof_model_area_based):
    roof_model_capacity_based = area_to_capacity(
        roof_model_area_based,
        power_density_flat=1,
        power_density_tilted=1
    )
    assert_frame_equal(roof_model_capacity_based, roof_model_area_based)


def test_weight_of_flat_reduced_for_lower_power_density(roof_model_area_based):
    roof_model_capacity_based = area_to_capacity(
        roof_model_area_based,
        power_density_flat=1,
        power_density_tilted=2
    )
    capacity_weight = float(roof_model_capacity_based.loc[("flat", 0.0)])
    area_weight = float(roof_model_area_based.loc[("flat", 0.0)])
    assert capacity_weight < area_weight


def test_weight_of_flat_increased_for_higher_power_density(roof_model_area_based):
    roof_model_capacity_based = area_to_capacity(
        roof_model_area_based,
        power_density_flat=2,
        power_density_tilted=1
    )
    capacity_weight = float(roof_model_capacity_based.loc[("flat", 0.0)])
    area_weight = float(roof_model_area_based.loc[("flat", 0.0)])
    assert capacity_weight > area_weight
