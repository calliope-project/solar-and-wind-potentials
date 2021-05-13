from enum import IntEnum, Enum
import numpy as np

DATATYPE = np.uint8

class Eligibility(IntEnum):
    """Categories defining land eligibility for renewable power."""
    NOT_ELIGIBLE = 0
    ROOFTOP_PV = 250
    ONSHORE_WIND_AND_PV = 180
    ONSHORE_WIND = 110
    OFFSHORE_WIND = 40

    @property
    def area_column_name(self):
        return "eligibility_{}_km2".format(self.name.lower())

    @staticmethod
    def onshore():
        """Returns all onshore eligibilities."""
        return [
            Eligibility.NOT_ELIGIBLE,
            Eligibility.ROOFTOP_PV,
            Eligibility.ONSHORE_WIND_AND_PV,
            Eligibility.ONSHORE_WIND,
        ]

    @staticmethod
    def offshore():
        """Returns all offshore eligibilities."""
        return [
            Eligibility.OFFSHORE_WIND
        ]


class GlobCover(IntEnum):
    """Original categories taken from GlobCover 2009 land cover."""
    POST_FLOODING = 11
    RAINFED_CROPLANDS = 14
    MOSAIC_CROPLAND = 20
    MOSAIC_VEGETATION = 30
    CLOSED_TO_OPEN_BROADLEAVED_FOREST = 40
    CLOSED_BROADLEAVED_FOREST = 50
    OPEN_BROADLEAVED_FOREST = 60
    CLOSED_NEEDLELEAVED_FOREST = 70
    OPEN_NEEDLELEAVED_FOREST = 90
    CLOSED_TO_OPEN_MIXED_FOREST = 100
    MOSAIC_FOREST = 110
    MOSAIC_GRASSLAND = 120
    CLOSED_TO_OPEN_SHRUBLAND = 130
    CLOSED_TO_OPEN_HERBS = 140
    SPARSE_VEGETATION = 150
    CLOSED_TO_OPEN_REGULARLY_FLOODED_FOREST = 160 # doesn't exist in Europe
    CLOSED_REGULARLY_FLOODED_FOREST = 170 # doesn't exist in Europe
    CLOSED_TO_OPEN_REGULARLY_FLOODED_GRASSLAND = 180 # roughly 2.3% of land in Europe
    ARTIFICAL_SURFACES_AND_URBAN_AREAS = 190
    BARE_AREAS = 200
    WATER_BODIES = 210
    PERMANENT_SNOW = 220
    NO_DATA = 230


FARM = [GlobCover.POST_FLOODING, GlobCover.RAINFED_CROPLANDS,
        GlobCover.MOSAIC_CROPLAND, GlobCover.MOSAIC_VEGETATION]
FOREST = [GlobCover.CLOSED_TO_OPEN_BROADLEAVED_FOREST, GlobCover.CLOSED_BROADLEAVED_FOREST,
          GlobCover.OPEN_BROADLEAVED_FOREST, GlobCover.CLOSED_NEEDLELEAVED_FOREST,
          GlobCover.OPEN_NEEDLELEAVED_FOREST, GlobCover.CLOSED_TO_OPEN_MIXED_FOREST,
          GlobCover.MOSAIC_FOREST, GlobCover.CLOSED_TO_OPEN_REGULARLY_FLOODED_FOREST,
          GlobCover.CLOSED_REGULARLY_FLOODED_FOREST]
VEGETATION = [GlobCover.MOSAIC_GRASSLAND, GlobCover.CLOSED_TO_OPEN_SHRUBLAND,
              GlobCover.CLOSED_TO_OPEN_HERBS, GlobCover.SPARSE_VEGETATION,
              GlobCover.CLOSED_TO_OPEN_REGULARLY_FLOODED_GRASSLAND]
BARE = [GlobCover.BARE_AREAS]
OTHER = VEGETATION + BARE
URBAN = [GlobCover.ARTIFICAL_SURFACES_AND_URBAN_AREAS]
WATER = [GlobCover.WATER_BODIES]


class ProtectedArea(IntEnum):
    """Derived from UNEP-WCMC data set."""
    PROTECTED = 255
    NOT_PROTECTED = 0


class Potential(Enum):
    """Classes of renewable electricity potentials."""
    ROOFTOP_PV = (1, [Eligibility.ROOFTOP_PV])
    OPEN_FIELD_PV = (2, [Eligibility.ONSHORE_WIND_AND_PV])
    ONSHORE_WIND = (3, [Eligibility.ONSHORE_WIND_AND_PV, Eligibility.ONSHORE_WIND])
    OFFSHORE_WIND = (4, [Eligibility.OFFSHORE_WIND])

    def __init__(self, int_id, corresponding_eligibilities):
        self.int_id = int_id
        self.eligible_on = corresponding_eligibilities

    @property
    def area_name(self):
        return "{}_km2".format(self.name.lower())

    @property
    def capacity_name(self):
        return "{}_mw".format(self.name.lower())

    @property
    def electricity_yield_name(self):
        return "{}_twh_per_year".format(self.name.lower())

    @staticmethod
    def onshore():
        """Returns all onshore potentials."""
        return [
            Potential.ROOFTOP_PV,
            Potential.OPEN_FIELD_PV,
            Potential.ONSHORE_WIND,
        ]

    @staticmethod
    def offshore():
        """Returns all offshore potentials."""
        return [
            Potential.OFFSHORE_WIND
        ]

    def __repr__(self):
        return self.electricity_yield_name

    def __str__(self):
        return self.__repr__()


def eligibility_land_mask(
    land_cover, slope_pv, slope_wind, bathymetry, building_share, urban_green_share,
    max_building_share, max_urban_green_share, max_depth_offshore, slope_threshold
):

    # parameters
    assert slope_pv <= slope_wind # wind can be built whereever pv can be built

    # prepare masks
    settlements = (building_share > max_building_share) | (urban_green_share > max_urban_green_share)
    farm = np.isin(land_cover, FARM)
    forest = np.isin(land_cover, FOREST)
    other = np.isin(land_cover, OTHER)
    water = np.isin(land_cover, WATER)
    pv = (slope_pv >= slope_threshold) & ~settlements & (farm | other)
    wind = (slope_wind >= slope_threshold) & ~settlements & (farm | forest | other)
    offshore = (bathymetry > max_depth_offshore) & ~settlements & water

    # allocate eligibility
    land = np.ones_like(land_cover, dtype=DATATYPE) * Eligibility.NOT_ELIGIBLE
    _add_eligibility(land, Eligibility.ROOFTOP_PV, settlements)
    _add_eligibility(land, Eligibility.ONSHORE_WIND_AND_PV, wind & pv)
    _add_eligibility(land, Eligibility.ONSHORE_WIND, wind & ~pv)
    _add_eligibility(land, Eligibility.OFFSHORE_WIND, offshore)
    return land


def _add_eligibility(land, eligibility, mask):
    assert all(land[mask] == Eligibility.NOT_ELIGIBLE), f"Overwriting other eligibility with {eligibility}."
    land[mask] = eligibility
