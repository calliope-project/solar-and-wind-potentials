"""Functions to convert units."""
def _set_proj_lib():
    # FIXME This is fragile and should not be necessary.
    import os
    from pathlib import Path
    path_to_projlib = Path(os.environ['CONDA_PREFIX']) / "share" / "proj"
    if not path_to_projlib.exists():
        msg = "Could not locate proj library. Path works on Unix only."
        raise ImportError(msg)
    os.environ["PROJ_LIB"] = path_to_projlib.as_posix()
_set_proj_lib()
from functools import partial
from itertools import product

import pycountry
import shapely.geometry
import shapely.ops
import pyproj

# from https://epsg.io/3035
EPSG_3035_PROJ4 = "+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs "


def watt_to_watthours(watt, duration):
    """Converts from [k|M|G|T]Watt to [k|M|G|T]WattHour."""
    return watt * duration.total_seconds() / 3600


def area_in_squaremeters(geodataframe):
    """Calculates the area sizes of a geo dataframe in square meters.

    Following https://gis.stackexchange.com/a/20056/77760 I am choosing equal-area projections
    to receive a most accurate determination of the size of polygons in the geo dataframe.
    Instead of Gall-Peters, as suggested in the answer, I am using EPSG_3035 which is
    particularly usefull for Europe.

    Returns a pandas series of area sizes in square meters.
    """
    return geodataframe.to_crs(EPSG_3035_PROJ4).area


def eu_country_code_to_iso3(eu_country_code):
    """Converts EU country code to ISO 3166 alpha 3.

    The European Union uses its own country codes, which often but not always match ISO 3166.
    """
    assert len(eu_country_code) == 2, "EU country codes are of length 2, yours is '{}'.".format(eu_country_code)
    if eu_country_code.lower() == "el":
        iso2 = "gr"
    elif eu_country_code.lower() == "uk":
        iso2 = "gb"
    else:
        iso2 = eu_country_code
    return pycountry.countries.lookup(iso2).alpha_3


def coordinate_string_to_decimal(coordinate_string):
    """Converts a coordinate string to decimal coordinates in degrees (easting, northing).

    A coordinate string looks something like this "48°18'N 14°17'E" hence is given in degrees,
    arcminutes, and potentially arcseconds. It may as well be given in decimals.

    The function would return the coordinates in decimal degrees easting and northing, so for
    the example given above that would be (14.283333, 48.300000).
    """
    lat, long = _separate_lat_and_long(coordinate_string)
    easting = _to_decimal_degree(*_split_coordinate(long))
    northing = _to_decimal_degree(*_split_coordinate(lat))
    return easting, northing


def transform_coordinates(x, y, from_epsg, to_epsg):
    """Tranforms coordinates from one coordinate reference system to the other."""
    point = shapely.geometry.Point(x, y)

    project = partial(
        pyproj.transform,
        pyproj.Proj(init=from_epsg),
        pyproj.Proj(init=to_epsg))

    transformed_point = shapely.ops.transform(project, point)
    return transformed_point.x, transformed_point.y


def transform_bounds(x_min, y_min, x_max, y_max, from_epsg, to_epsg):
    """Tranforms bounds from one coordinate reference system to the other.

    Returns bounds as tuple: (x_min, y_min, x_max, y_max).
    """
    points = [transform_coordinates(x, y, from_epsg, to_epsg)
              for x, y in product([x_min, x_max], [y_min, y_max])]

    return (min([p[0] for p in points]), min([p[1] for p in points]),
            max([p[0] for p in points]), max([p[1] for p in points]),)


def _separate_lat_and_long(coordinate_string):
    coordinate_string = coordinate_string.replace(",", "")
    assert "N" in coordinate_string, coordinate_string
    length_of_lat = coordinate_string.find("N") + 1
    lat = coordinate_string[:length_of_lat].strip()
    long = coordinate_string[length_of_lat:].strip()
    assert lat[-1] == "N", lat
    assert long[-1] in ["E", "O", "W"], long
    return lat, long


def _split_coordinate(coordinate_string):
    coordinate_string = coordinate_string.replace("′", "'")
    if "°" in coordinate_string: # coordinates are given in arcminutes
        degrees, residual = coordinate_string.split("°")
        if "'" in residual:
            arcminutes, residual = residual.split("'")
        else:
            arcminutes = 0.0
        if '"' in residual:
            arcseconds, residual = residual.split('"')
        else:
            arcseconds = 0.0
    else: # coordinates are given in decimals already
        degrees = coordinate_string[:-1]
        residual = coordinate_string[-1]
        arcminutes = 0.0
        arcseconds = 0.0
    degrees, arcminutes, arcseconds = [float(value) for value in (degrees, arcminutes, arcseconds)]
    if "W" in residual:
        # longitude given in westings instead of eastings
        degrees, arcminutes, arcseconds = [- value for value in (degrees, arcminutes, arcseconds)]
    return degrees, arcminutes, arcseconds


def _to_decimal_degree(degrees, arcminutes, arcseconds):
    return degrees + arcminutes / 60 + arcseconds / 3600