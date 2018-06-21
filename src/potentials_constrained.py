"""Quantify the constrained potentials for renewable power in regions.

Based on the unconstrained potentials and rules to constrain it.
"""
import click
import pandas as pd

from src.utils import Config
from src.eligible_land import Eligibility


@click.command()
@click.argument("path_to_unconstrained_potentials_prefer_pv")
@click.argument("path_to_unconstrained_potentials_prefer_wind")
@click.argument("path_to_result")
@click.argument("scenario")
@click.argument("config", type=Config())
def constrained_potentials(path_to_unconstrained_potentials_prefer_pv, path_to_unconstrained_potentials_prefer_wind,
                           path_to_result, scenario, config):
    unconstrained_prefer_pv = pd.read_csv(path_to_unconstrained_potentials_prefer_pv, index_col=0)
    unconstrained_prefer_wind = pd.read_csv(path_to_unconstrained_potentials_prefer_wind, index_col=0)

    constrained = _constrain_potential(unconstrained_prefer_pv, unconstrained_prefer_wind,
                                       config["scenarios"][scenario])
    constrained.to_csv(path_to_result, header=True)


def _constrain_potential(unconstrained_prefer_pv, unconstrained_prefer_wind, scenario_config):
    _assert_pv_has_higher_energy_density(unconstrained_prefer_pv, unconstrained_prefer_wind)
    factor_pv = pd.DataFrame(
        index=unconstrained_prefer_pv.index,
        data={eligibility.energy_column_name: _scaling_factor(eligibility, scenario_config)[0]
              for eligibility in Eligibility}
    )
    factor_wind = pd.DataFrame(
        index=unconstrained_prefer_wind.index,
        data={eligibility.energy_column_name: _scaling_factor(eligibility, scenario_config)[1]
              for eligibility in Eligibility}
    )
    return factor_pv * unconstrained_prefer_pv + factor_wind * unconstrained_prefer_wind


def _scaling_factor(eligibility, scenario_config):
    # Returns a tuple: first for prefer_pv, second for prefer_wind
    # Uses the assumption that pv's energy density is higher than the one from onshore wind.
    share_protected_areas_used = scenario_config["share-protected-areas-used"]
    share_rooftops_used = scenario_config["share-rooftops-used"]
    share_other_land_used = scenario_config["share-other-land-used"]
    share_farmland_used = scenario_config["share-farmland-used"]
    share_forest_used_for_wind = scenario_config["share-forest-used-for-wind"]
    share_pv_on_farmland = scenario_config["share-pv-on-farmland"]
    share_remain_on_farmland = share_farmland_used - share_pv_on_farmland
    share_offshore_used = scenario_config["share-offshore-used"]
    assert share_pv_on_farmland <= share_farmland_used
    return {
        Eligibility.NOT_ELIGIBLE: [1, 0],
        Eligibility.ROOFTOP_PV: [share_rooftops_used, 0],
        Eligibility.ONSHORE_WIND_AND_PV_OTHER: [share_other_land_used, 0],
        Eligibility.ONSHORE_WIND_OTHER: [share_other_land_used, 0],
        Eligibility.ONSHORE_WIND_FARMLAND: [share_farmland_used, 0],
        Eligibility.ONSHORE_WIND_FOREST: [share_forest_used_for_wind, 0],
        Eligibility.ONSHORE_WIND_AND_PV_FARMLAND: [share_pv_on_farmland, share_remain_on_farmland],
        Eligibility.OFFSHORE_WIND: [share_offshore_used, 0],
        Eligibility.ONSHORE_WIND_AND_PV_OTHER_PROTECTED: [share_other_land_used * share_protected_areas_used, 0],
        Eligibility.ONSHORE_WIND_OTHER_PROTECTED: [share_other_land_used * share_protected_areas_used, 0],
        Eligibility.ONSHORE_WIND_FARMLAND_PROTECTED: [share_farmland_used * share_protected_areas_used, 0],
        Eligibility.ONSHORE_WIND_FOREST_PROTECTED: [share_forest_used_for_wind * share_protected_areas_used, 0],
        Eligibility.ONSHORE_WIND_AND_PV_FARMLAND_PROTECTED: [share_pv_on_farmland * share_protected_areas_used,
                                                             share_remain_on_farmland * share_protected_areas_used],
        Eligibility.OFFSHORE_WIND_PROTECTED: [share_offshore_used * share_protected_areas_used, 0]
    }[eligibility]


def _assert_pv_has_higher_energy_density(prefer_pv, prefer_wind):
    # In the following I am assuming that the energy density of PV is always better than
    # wind onshore and hence PV is always preferred over wind.
    assert (prefer_pv >= prefer_wind).all().all()


if __name__ == "__main__":
    constrained_potentials()