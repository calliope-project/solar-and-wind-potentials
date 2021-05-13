from pathlib import Path
import os

import pytest
import fiona
import yaml

ROOT_DIR = Path(os.path.abspath(__file__)).parent.parent

with open('config/default.yaml', 'r') as src:
    CONFIG_DEFAULT = yaml.safe_load(src)

PATH_TO_BORDERS = ROOT_DIR / "build" / "administrative-borders.gpkg"


@pytest.mark.skipif(not PATH_TO_BORDERS.exists(), reason="Consolidated administrative border shapefile not available.")
def test_administrative_border_layers():
    layers = fiona.listlayers(str(PATH_TO_BORDERS))
    default_layers = set(
        v for countries in CONFIG_DEFAULT['layers'].values() for k, v in countries.items()
    )
    assert len(default_layers.difference(layers)) == 0
