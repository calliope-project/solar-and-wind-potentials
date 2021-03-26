import os
import yaml
from pathlib import Path

import geopandas as gpd
import pytest
import shapely.geometry
import pycountry
import fiona

from src.administrative_borders import _study_area, _to_multi_polygon, _drop_countries, _drop_geoms_completely_outside_study_area, _drop_parts_of_geoms_completely_outside_study_area
from src.conversion import eu_country_code_to_iso3
# Loading the 60M file for its smaller size (1M used in actual workflow)
URL_NUTS = "https://ec.europa.eu/eurostat/cache/GISCO/distribution/v2/nuts/shp/NUTS_RG_60M_{}_4326.geojson"
ROOT_DIR = Path(os.path.abspath(__file__)).parent.parent
PATH_TO_BORDERS = ROOT_DIR / "build" / "administrative-borders.gpkg"

with open('config/schema.yaml', 'r') as src:
    config_schema = yaml.safe_load(src)

with open('config/default.yaml', 'r') as src:
    config_default = yaml.safe_load(src)

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
        URL_NUTS.format(config_default["parameters"]["nuts-year"])
    )

@pytest.mark.parametrize("year", config_schema["properties"]["parameters"]["properties"]["nuts-year"]["enum"])
def test_nuts_years(year):
    gdf = gpd.read_file(URL_NUTS.format(year))
    cols = ['id', 'FID', 'LEVL_CODE', 'NUTS_NAME']
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
    assert not all([isinstance(i, shapely.geometry.Polygon) for i in gdf.geometry])

    gdf = gdf.map(_to_multi_polygon)

    assert all([isinstance(i, shapely.geometry.MultiPolygon) for i in gdf.geometry])


def test_drop_countries(europe_gdf):
    config = {
        "scope": {"countries": ["Germany", "France", "Greece"]}
    }
    gdf = europe_gdf().rename(columns={'CNTR_CODE': 'country_code'})
    gdf["country_code"] = gdf["country_code"].apply(eu_country_code_to_iso3)

    assert set(gdf.country_code) != set(
        [pycountry.countries.lookup(i).alpha_3 for i in config["scope"]["countries"]]
    )

    gdf = _drop_countries(gdf, config)

    assert set(gdf.country_code) == set(
        [pycountry.countries.lookup(i).alpha_3 for i in config["scope"]["countries"]]
    )


def test_drop_geoms(capsys, europe_gdf):
    config = {
        "scope": {
            'bounds': {  # IRE and GBR, not including shetland and some of cornwall
                'x_min': -11.12,
                'x_max': 2.24,
                'y_min': 51.11,
                'y_max': 59.72,
            }
        }
    }

    with capsys.readouterr() as captured:
        gdf = _drop_geoms_completely_outside_study_area(europe_gdf, config)
    not_dropped = ['IRL', 'North West (England)', 'United Kingdom', 'Ireland']
    for i in not_dropped:
        assert i not in captured.out
    dropped = [
        "Removing France (nuts0, country=FRA) as they are outside of study area",
        "(FRA)",
        "(GBR)",
    ]
    for i in dropped:
        assert i in captured.out

    with capsys.readouterr() as captured:
        gdf = _drop_parts_of_geoms_completely_outside_study_area(gdf, config)
    for i in not_dropped:
        assert i not in captured.out
    dropped = [
        "Removing parts of Shetland Islands (nuts3) as they are outside of study area",
        "(GBR)",
    ]
    for i in dropped:
        assert i in captured.out
    assert sorted(gdf.country_code.unique()) == sorted(['IRL', 'GBR'])

def test_administrative_border_layers(config_default):
    layers = fiona.listlayers(PATH_TO_BORDERS)
    default_layers = set(
        v for countries in config_default['layers'].values() for k, v in countries.items()
    )
    assert len(default_layers.difference(layers)) == 0
