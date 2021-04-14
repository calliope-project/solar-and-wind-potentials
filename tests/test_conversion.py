from datetime import timedelta
import math

import pytest

from scripts.conversion import watt_to_watthours, eu_country_code_to_iso3, coordinate_string_to_decimal,\
    transform_coordinates


@pytest.mark.parametrize("watt,duration,expected_watthour", [
    (10, timedelta(minutes=60), 10),
    (6, timedelta(minutes=30), 3),
    (1 / 8760, timedelta(days=365), 1)
])
def test_watt_to_watthour_conversion(watt, duration, expected_watthour):
    assert watt_to_watthours(watt=watt, duration=duration) == expected_watthour


@pytest.mark.parametrize(
    "eu_country_code,iso3",
    [("DE", "DEU"),
     ("EL", "GRC"),
     ("UK", "GBR")]
)
def test_eu_country_code(eu_country_code, iso3):
    assert eu_country_code_to_iso3(eu_country_code) == iso3


@pytest.mark.parametrize(
    "arcminutes,expected_easting,expected_northing",
    [("""48°18'N 14°17'E""", 14.283333, 48.300000),
     ("""54°35'20.0"N 1°11'15.0"W""", -1.187500, 54.588889)]
)
def test_coordinates_to_decimal(arcminutes, expected_easting, expected_northing):
    easting, northing = coordinate_string_to_decimal(arcminutes)
    assert math.isclose(easting, expected_easting, abs_tol=0.00001)
    assert math.isclose(northing, expected_northing, abs_tol=0.00001)


@pytest.mark.parametrize(
    "arcminutes,expected_easting,expected_northing",
    [("""48°18'N 14°17'O""", 14.283333, 48.300000),
     ("""48°18'N 14°17' O""", 14.283333, 48.300000),
     ("""48°18'N, 14°17' O""", 14.283333, 48.300000),
     ("""48.300000 N, 14.283333O""", 14.283333, 48.300000),
     ("""48°18′N 14°17′E""", 14.283333, 48.300000),
     ("""48°18′0.0″N 14°17′0.0″E""", 14.283333, 48.300000)]
)
def test_coordinates_to_decimal_edgecases(arcminutes, expected_easting, expected_northing):
    easting, northing = coordinate_string_to_decimal(arcminutes)
    assert math.isclose(easting, expected_easting, abs_tol=0.00001)
    assert math.isclose(northing, expected_northing, abs_tol=0.00001)


@pytest.mark.parametrize(
    "from_epsg,from_x,from_y,to_x,to_y",
    [("EPSG:4326", 8.55, 47.36, 4211389.55, 2695117.37), # values from epsg.io
     ("EPSG:4326", 33.87, 89.44, 4347749.36, 7315609.95)] # values from epsg.io
)
def test_transform_coordinates_to_epsg3035(from_epsg, from_x, from_y, to_x, to_y):
    x, y = transform_coordinates(
        from_epsg=from_epsg,
        to_epsg="EPSG:3035",
        x=from_x,
        y=from_y
    )
    assert math.isclose(x, to_x, abs_tol=0.1) # tolerance 0.1m
    assert math.isclose(y, to_y, abs_tol=0.1)


@pytest.mark.parametrize(
    "from_epsg,from_x,from_y,to_x,to_y",
    [("EPSG:4326", 92.4, 32.8, 8315488.22495176, 3969803.31307849), # values from epsg.io
     ("EPSG:4326", 33.87, 89.44, 0.00000000, 9020047.84807365)] # values from epsg.io
)
@pytest.mark.xfail(reason="ESRI:54009 may not be the right identifier.")
def test_transform_coordinates_to_esri54009(from_epsg, from_x, from_y, to_x, to_y):
    x, y = transform_coordinates(
        from_epsg=from_epsg,
        to_epsg="ESRI:54009",
        x=from_x,
        y=from_y
    )
    assert math.isclose(x, to_x, abs_tol=0.1) # tolerance 0.1m
    assert math.isclose(y, to_y, abs_tol=0.1)


@pytest.mark.parametrize(
    "from_epsg,to_x,to_y,from_x,from_y",
    [("EPSG:3035", 8.55, 47.36, 4211389.55, 2695117.37), # values from epsg.io
     ("EPSG:3035", 33.87, 89.44, 4347749.36, 7315609.95)] # values from epsg.io
)
def test_transform_coordinates_to_epsg4326(from_epsg, to_x, to_y, from_x, from_y):
    x, y = transform_coordinates(
        from_epsg=from_epsg,
        to_epsg="EPSG:4326",
        x=from_x,
        y=from_y
    )
    assert math.isclose(x, to_x, abs_tol=0.01)
    assert math.isclose(y, to_y, abs_tol=0.01)
