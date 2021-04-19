"""Create PV simulation input for renewables.ninja."""
import pandas as pd
import geopandas as gpd

from renewablepotentialslib.shape_utils import point_raster_on_shapes


def pv_simulation_parameters(path_to_shapes_of_land_surface, path_to_roof_categories,
                             bounds, ninja, maximum_power_density, path_to_output):
    """Create PV simulation input for renewables.ninja."""
    points = point_raster_on_shapes(
        bounds_wgs84=bounds,
        shapes=gpd.read_file(path_to_shapes_of_land_surface),
        resolution_km2=ninja["resolution-grid"]
    )

    roof_categories = pd.read_csv(path_to_roof_categories, index_col=[0, 1])
    roof_categories = area_to_capacity(
        roof_categories,
        power_density_flat=maximum_power_density["pv-on-flat-areas"],
        power_density_tilted=maximum_power_density["pv-on-tilted-roofs"]
    ).reset_index()
    lat_long = pd.DataFrame(
        data={
            "lat": [point.y for point in points.geometry],
            "long": [point.x for point in points.geometry]
        }
    )

    index = pd.MultiIndex.from_product((points.index, roof_categories.index), names=["id", "roof_cat_id"])
    data = pd.DataFrame(index=index).reset_index()
    data = data.merge(roof_categories, left_on="roof_cat_id", right_index=True).rename(
        columns={"share_of_roof_areas": "weight"}
    )
    data = data.merge(lat_long, left_on="id", right_index=True)
    data["azim"] = data["orientation"].map(orientation_to_azimuth)
    data["site_id"] = data.id
    data["sim_id"] = data.apply(
        lambda row: "{}_{}_{}".format(row.id, row.orientation, round(row.average_tilt)),
        axis=1
    )
    flat_mask = data["orientation"] == "flat"
    data.loc[flat_mask, "average_tilt"] = data.loc[flat_mask, "lat"].map(optimal_tilt)
    data["pr"] = ninja["pv-performance-ratio"]
    data[
        ["sim_id", "weight", "site_id", "lat", "long", "average_tilt",
         "orientation", "azim", "pr"]
    ].sort_index().to_csv(
        path_to_output,
        header=True,
        index=False
    )


def orientation_to_azimuth(orientation):
    if orientation == "S":
        return 180
    elif orientation == "W":
        return -90
    elif orientation == "N":
        return 0
    elif orientation == "E":
        return 90
    elif orientation == "flat":
        return 180
    else:
        raise ValueError()


def optimal_tilt(latitude):
    # based on @Jacobson:2018
    optimal_tilt = 1.3793 + latitude * (1.2011 + latitude * (-0.014404 + latitude * 0.000080509))
    assert 90 > optimal_tilt >= 0
    return optimal_tilt


def area_to_capacity(statistical_roof_model_area_based, power_density_flat, power_density_tilted):
    """Maps area shares to capacity shares of statistical roof model.

    The statistical roof model defines roof categories (e.g. south-facing with tilt 10°) and their
    shares in a population of roofs. This function maps areas shares to shares of installable pv
    capacity. It discriminates between flat and tilted roofs, i.e. the power density of flat roofs
    can be different than the one from tilted roofs.

    Parameters:
        * statistical_roof_model_area_based: model as described above, values are shares of total roof area
        * power_density_flat: power density of flat pv installations, unit must be consistent with next
        * power_density_tilted: power density of tilted pv installations, unit must be consistent with previous
    Returns:
        * statistical_roof_model_cpacity_based: model as described above, values are shares of total
        installable capacity
    """
    cap_based = statistical_roof_model_area_based.copy()
    flat_roofs = cap_based.index.get_level_values(0) == "flat"
    tilted_roofs = cap_based.index.get_level_values(0) != "flat"
    cap_based[flat_roofs] = cap_based[flat_roofs] * power_density_flat
    cap_based[tilted_roofs] = cap_based[tilted_roofs] * power_density_tilted
    return cap_based / cap_based.sum()


if __name__ == "__main__":
    pv_simulation_parameters(
        path_to_shapes_of_land_surface=snakemake.input.units,
        path_to_roof_categories=snakemake.input.roof_categories,
        bounds=snakemake.params.bounds,
        ninja=snakemake.params.ninja,
        maximum_power_density=snakemake.params.maximum_power_density,
        path_to_output=snakemake.output.points
    )
