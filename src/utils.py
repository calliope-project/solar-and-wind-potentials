"""Module containing utilities."""
import math
from pathlib import Path

import click
import numpy as np
import rasterio
import yaml

PATH_TO_CONFIGS = Path(__file__).parent / '..' / 'config'


class Config(click.ParamType):
    """A configuration parameter on the command line.

    Configurations will always be read from the config directory.
    """
    name = "configuration"

    def convert(self, value, param, ctx):
        name_of_file = Path(value).name
        path_to_file = PATH_TO_CONFIGS / name_of_file
        return read_config(path_to_file)


def read_config(path_to_file):
    """Reads a configuration file."""
    path_to_file = Path(path_to_file)
    if not path_to_file.exists():
        raise ValueError("Config {} does not exist.".format(path_to_file))
    with path_to_file.open('r') as stream:
        try:
            return yaml.safe_load(stream)
        except yaml.YAMLError as exc:
            raise IOError(exc)


def determine_pixel_areas(crs, bounds, resolution):
    """Returns a raster in which the value corresponds to the area [km2] of the pixel.

    This assumes the data comprises square pixel in WGS84.

    Parameters:
        crs: the coordinate reference system of the data (must be WGS84)
        bounds: an object with attributes left/right/top/bottom given in degrees
        resolution: the scalar resolution (remember: square pixels) given in degrees
    """
    # the following is based on https://gis.stackexchange.com/a/288034/77760
    # and assumes the data to be in WGS84
    assert crs == rasterio.crs.CRS.from_epsg("4326") # WGS84
    width = int((bounds.right - bounds.left) / resolution)
    height = int((bounds.top - bounds.bottom) / resolution)
    latitudes = np.linspace(
        start=bounds.top,
        stop=bounds.bottom,
        num=height,
        endpoint=True,
        dtype=np.float64
    )
    varea_of_pixel = np.vectorize(lambda lat: _area_of_pixel(resolution, lat))
    pixel_area = varea_of_pixel(latitudes) # vector
    return pixel_area.repeat(width).reshape(height, width).astype(np.float64)


def _area_of_pixel(pixel_size, center_lat):
    """Calculate km^2 area of a wgs84 square pixel.

    Adapted from: https://gis.stackexchange.com/a/127327/2397

    Parameters:
        pixel_size (float): length of side of pixel in degrees.
        center_lat (float): latitude of the center of the pixel. Note this
            value +/- half the `pixel-size` must not exceed 90/-90 degrees
            latitude or an invalid area will be calculated.

    Returns:
        Area of square pixel of side length `pixel_size` centered at
        `center_lat` in km^2.

    """
    a = 6378137  # meters
    b = 6356752.3142  # meters
    e = math.sqrt(1 - (b / a)**2)
    area_list = []
    for f in [center_lat + pixel_size / 2, center_lat - pixel_size / 2]:
        zm = 1 - e * math.sin(math.radians(f))
        zp = 1 + e * math.sin(math.radians(f))
        area_list.append(
            math.pi * b**2 * (
                math.log(zp / zm) / (2 * e) +
                math.sin(math.radians(f)) / (zp * zm)))
    return pixel_size / 360. * (area_list[0] - area_list[1]) / 1e6


def buffer_if_necessary(shape):
    """Fix shapes which are invalid.

    Following the advice given here:
    https://github.com/Toblerity/Shapely/issues/344
    """
    if shape.is_valid:
        return shape

    new_shape = shape.buffer(0.0)
    assert new_shape.is_valid
    assert np.isclose(new_shape.area, new_shape.area)

    return new_shape