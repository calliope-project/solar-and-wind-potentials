import yaml

import geopandas as gpd
import pytest
import shapely.geometry

from src.administrative_borders import _study_area, _to_multi_polygon
# Loading the 20M file for its smaller size (1M used in actual workflow)
URL_NUTS = "https://ec.europa.eu/eurostat/cache/GISCO/distribution/v2/nuts/shp/NUTS_RG_20M_{}_4326.shp.zip"

with open('config/schema.yaml', 'r') as src:
    config_schema = yaml.safe_load(src)

@pytest.fixture
def config(exclusions=0):
    _config = {
        "scope": {
            'bounds': {
                'x_min': 0,
                'x_max': 1,
                'y_min': 50,
                'y_max': 51,
            }
        }
    }
    if exclusions > 0:
        _config["scope"]["exclusion_zones"] = {
            'foo': {
                'x_min': 0.5,
                'x_max': 0.6,
                'y_min': 50.2,
                'y_max': 50.5,
            }
        }
    if exclusions == 2:
        _config["scope"]["exclusion_zones"]["bar"] = {
            'x_min': 0,
            'x_max': 0.1,
            'y_min': 50,
            'y_max': 50.1,
        }
    return _config

@pytest.fixture
def europe_gdf():
    return gpd.read_file(
        URL_NUTS.format(config_schema["properties"]["parameters"]["properties"]["nuts-year"]["default"])
    )

@pytest.mark.parametrize("year", config_schema["properties"]["parameters"]["properties"]["nuts-year"]["enum"])
def test_nuts_years(year):
    gdf = gpd.read_file(URL_NUTS.format(year))
    cols = ['FID', 'NUTS_ID', 'LEVL_CODE', 'NUTS_NAME']
    assert all(i in gdf.columns for i in cols)


def test_bounding_box(config):
    area = _study_area(config())
    assert area.bounds == (0.0, 50.0, 1.0, 51.0)

@pytest.mark.parametrize(
    "exclusions,not_in,is_in",
    [(1, (0.55, 50.3), (0.7, 50.7)), (2, (0.05, 50.05), (0.05, 50.3))]
)
def test_exclusion_zone(exclusions, not_in, is_in, config):
    no_exclusions = _study_area(config())
    area = _study_area(config(exclusions))

    assert no_exclusions.area > area.area
    assert area.area == no_exclusions.area
    assert area.bounds == no_exclusions.bounds

    assert not shapely.geometry.Point(not_in).within(area)
    assert shapely.geometry.Point(is_in).within(area)

def test_multipolygon(europe_gdf):
    gdf = europe_gdf()
    assert any([isinstance(i, shapely.geometry.Polygon) for i in gdf.geometry])

    gdf = gdf.map(_to_multi_polygon)

    assert all([isinstance(i, shapely.geometry.MultiPolygon) for i in gdf.geometry])