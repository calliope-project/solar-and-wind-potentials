import yaml

import geopandas as gpd
import pytest

# Loading the 20M file for its smaller size (1M used in actual workflow)
URL_NUTS = "https://ec.europa.eu/eurostat/cache/GISCO/distribution/v2/nuts/shp/NUTS_RG_20M_{}_4326.shp.zip"

with open('config/schema.yaml', 'r') as src:
    config_schema = yaml.safe_load(src)

@pytest.mark.parametrize("year", config_schema["properties"]["parameters"]["properties"]["nuts-year"]["enum"])
def test_nuts_years(year):
    gdf = gpd.read_file(URL_NUTS.format(year))
    cols = ['FID', 'NUTS_ID', 'LEVL_CODE', 'NUTS_NAME']
    assert all(i in gdf.columns for i in cols)
