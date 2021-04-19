
from pathlib import Path
import os

import pytest
import fiona
import yaml

ROOT_DIR = Path(os.path.abspath(__file__)).parent.parent

with open('config/default.yaml', 'r') as src:
    CONFIG_DEFAULT = yaml.safe_load(src)


BUILD_DIR = Path(os.path.abspath(__file__)).parent.parent / "build"
PATH_TO_BORDERS = BUILD_DIR / "administrative-borders.gpkg"
@pytest.mark.skipif(not PATH_TO_BORDERS.exists(), reason="Consolidated administrative border shapefile not available.")
class TestFinalBorders():
    def test_administrative_border_layers(self):
        layers = fiona.listlayers(self.PATH_TO_BORDERS)
        default_layers = set(
            v for countries in CONFIG_DEFAULT['layers'].values() for k, v in countries.items()
        )
        assert len(default_layers.difference(layers)) == 0
