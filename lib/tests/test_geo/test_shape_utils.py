import geopandas as gpd
import pytest
import shapely.geometry
import pycountry

from renewablepotentialslib.geo.shape_utils import (
    study_area,
    to_multi_polygon,
    drop_countries,
    drop_geoms_completely_outside_study_area,
    drop_parts_of_geoms_completely_outside_study_area,
    update_features
)

# Loading the 60M file for its smaller size (1M used in actual workflow)
URL_NUTS = "https://ec.europa.eu/eurostat/cache/GISCO/distribution/v2/nuts/geojson/NUTS_RG_60M_{}_4326.geojson"

@pytest.mark.parametrize("year", [2006, 2010, 2013, 2016, 2021])
def test_nuts_years(year):
    gdf = gpd.read_file(URL_NUTS.format(year))
    cols = ['id', 'FID', 'LEVL_CODE', 'NUTS_NAME']
    assert all(i in gdf.columns for i in cols)

class TestStudyArea():
    def config(self, exclusions=0):
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

    def test_bounding_box(self):
        area = study_area(self.config()["scope"])
        assert area.bounds == (0.0, 50.0, 1.0, 51.0)


    @pytest.mark.parametrize(
        "exclusions,not_in,is_in",
        [(1, (0.55, 50.3), (0.7, 50.7)), (2, (0.05, 50.05), (0.05, 50.3))]
    )
    def test_exclusion_zone(self, exclusions, not_in, is_in):
        no_exclusions = study_area(self.config()["scope"])
        with_exclusions = study_area(self.config(exclusions)["scope"])

        assert no_exclusions.area > with_exclusions.area
        assert with_exclusions.bounds == no_exclusions.bounds

        assert not shapely.geometry.Point(not_in).within(with_exclusions)
        assert shapely.geometry.Point(is_in).within(with_exclusions)


class TestGeomManipulation():
    CONFIG = {
        'bounds': {  # IRE and GBR, not including shetland and some of cornwall
            'x_min': -10,
            'x_max': 1,
            'y_min': 51.11,
            'y_max': 60.1,
        }
    }
    NOT_DROPPED = ['IRL', 'North West (England)']
    FULL_DROPPED = [
        "Removing France (nuts0, country=FRA) as they are outside of study area",
    ]
    PARTIAL_DROPPED = [
        "Removing parts of Highlands and Islands (nuts2, country=GBR)",
        "Removing parts of United Kingdom (nuts0, country=GBR)",
    ]

    @pytest.fixture
    def europe_gdf(self):
        gdf = gpd.read_file(
            URL_NUTS.format(2016)
        )
        return update_features(gdf, "nuts")

    @pytest.mark.filterwarnings("ignore:Geometry is in a geographic CRS:UserWarning")
    def test_multipolygon(self, europe_gdf):
        assert not all([isinstance(i, shapely.geometry.Polygon) for i in europe_gdf.geometry])

        polys = europe_gdf.geometry.map(to_multi_polygon)

        assert all([isinstance(i, shapely.geometry.MultiPolygon) for i in polys])
        assert all(polys.area == europe_gdf.area)

    def test_drop_countries(self, europe_gdf):
        config = {"countries": ["Germany", "France", "Greece"]}

        assert set(europe_gdf.country_code) != set(
            [pycountry.countries.lookup(i).alpha_3 for i in config["countries"]]
        )

        gdf = drop_countries(europe_gdf, config)

        assert set(gdf.country_code) == set(
            [pycountry.countries.lookup(i).alpha_3 for i in config["countries"]]
        )

    def test_drop_geoms(self, capsys, europe_gdf):
        gdf = drop_geoms_completely_outside_study_area(europe_gdf, self.CONFIG)
        out, err = capsys.readouterr()

        for i in self.NOT_DROPPED + self.PARTIAL_DROPPED:
            assert i not in out

        for i in self.FULL_DROPPED:
            assert i in out
        assert sorted(gdf.country_code.unique()) == sorted(['IRL', 'GBR'])

    def test_drop_parts_of_geoms(self, capsys, europe_gdf):
        gdf = drop_parts_of_geoms_completely_outside_study_area(europe_gdf, self.CONFIG)
        out, err = capsys.readouterr()
        for i in self.NOT_DROPPED + self.FULL_DROPPED:
            assert i not in out
        for i in self.PARTIAL_DROPPED:
            assert i in out
        assert (
            europe_gdf[europe_gdf.name == 'Highlands and Islands'].area
            > gdf[gdf.name == 'Highlands and Islands'].area
        ).all()
