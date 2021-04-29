"""Test whether potential estimations between layers are similar."""
import os
from pathlib import Path

import pytest
import pandas as pd

from renewablepotentialslib.eligibility import Potential

TOLERANCE = 0.005 # 0.5%

BUILD_DIR = Path(os.path.abspath(__file__)).parent.parent / "build"
PATH_TO_CONTINENTAL_POTENTIALS = BUILD_DIR / "continental" / "technical-potential" / "potentials.csv"
PATH_TO_NATIONAL_POTENTIALS = BUILD_DIR / "national" / "technical-potential" / "potentials.csv"
PATH_TO_REGIONAL_POTENTIALS = BUILD_DIR / "regional" / "technical-potential" / "potentials.csv"
PATH_TO_MUNICIPAL_POTENTIALS = BUILD_DIR / "municipal" / "technical-potential" / "potentials.csv"


@pytest.mark.skipif(not PATH_TO_CONTINENTAL_POTENTIALS.exists(), reason="Continental potentials not available.")
@pytest.mark.skipif(not PATH_TO_NATIONAL_POTENTIALS.exists(), reason="National potentials not available.")
@pytest.mark.parametrize(
    "potential", Potential
)
def test_continental_to_national(potential):
    continental = pd.read_csv(PATH_TO_CONTINENTAL_POTENTIALS, index_col=0).sum()
    national = pd.read_csv(PATH_TO_NATIONAL_POTENTIALS, index_col=0).sum()

    assert continental[str(potential)] == pytest.approx(national[str(potential)], TOLERANCE)


@pytest.mark.skipif(not PATH_TO_CONTINENTAL_POTENTIALS.exists(), reason="Continental potentials not available.")
@pytest.mark.skipif(not PATH_TO_REGIONAL_POTENTIALS.exists(), reason="Regional potentials not available.")
@pytest.mark.parametrize(
    "potential", Potential
)
def test_continental_to_regional(potential):
    continental = pd.read_csv(PATH_TO_CONTINENTAL_POTENTIALS, index_col=0).sum()
    regional = pd.read_csv(PATH_TO_REGIONAL_POTENTIALS, index_col=0).sum()

    assert continental[str(potential)] == pytest.approx(regional[str(potential)], TOLERANCE)


@pytest.mark.skipif(not PATH_TO_CONTINENTAL_POTENTIALS.exists(), reason="Continental potentials not available.")
@pytest.mark.skipif(not PATH_TO_MUNICIPAL_POTENTIALS.exists(), reason="Municipal potentials not available.")
@pytest.mark.parametrize(
    "potential", Potential
)
def test_continental_to_municipal(potential):
    continental = pd.read_csv(PATH_TO_CONTINENTAL_POTENTIALS, index_col=0).sum()
    municipal = pd.read_csv(PATH_TO_MUNICIPAL_POTENTIALS, index_col=0).sum()

    assert continental[str(potential)] == pytest.approx(municipal[str(potential)], TOLERANCE)
