import click
import numpy as np
import pandas as pd
from SALib.sample import saltelli
from SALib.analyze import sobol
from SALib.util import compute_groups_matrix

from src.potentials_constrained import _constrain_potential


@click.command()
@click.argument("path_to_unconstrained_potentials_prefer_pv")
@click.argument("path_to_unconstrained_potentials_prefer_wind")
@click.argument("path_to_demand")
@click.argument("path_to_population")
@click.argument("path_to_result")
def sensititivity_analysis(path_to_unconstrained_potentials_prefer_pv, path_to_unconstrained_potentials_prefer_wind,
                           path_to_demand, path_to_population, path_to_result):
    demand = pd.read_csv(path_to_demand, index_col=0)["demand_twh_per_year"]
    population = pd.read_csv(path_to_population, index_col=0)["population_sum"].reindex(demand.index)
    unconstrained_prefer_pv = pd.read_csv(
        path_to_unconstrained_potentials_prefer_pv,
        index_col=0
    ).reindex(demand.index)
    unconstrained_prefer_wind = pd.read_csv(
        path_to_unconstrained_potentials_prefer_wind,
        index_col=0
    ).reindex(demand.index)

    problem = {
        'num_vars': 7,
        'names': ['share-protected-areas-used',
                  'share-pv-on-farmland',
                  'share-farmland-used',
                  'share-forest-used-for-wind',
                  'share-other-land-used',
                  'share-offshore-used',
                  'share-rooftops-used'],
        'bounds': [[0, 1],
                   [0, 1],
                   [0, 1],
                   [0, 1],
                   [0, 1],
                   [0, 1],
                   [0, 1]]
    }
    param_values = saltelli.sample(problem, 100)

    print("Start evaluation")
    Y = np.zeros([param_values.shape[0]])
    for i, X in enumerate(param_values):
        Y[i] = evaluate_model(unconstrained_prefer_pv, unconstrained_prefer_wind, demand, population, X)

    print("Start analysis")
    Si = sobol.analyze(problem, Y)

    with open(path_to_result, "w") as f_out:
        print_indices(Si, problem, True, f_out)


def evaluate_model(unconstrained_prefer_pv, unconstrained_prefer_wind, demand, population, x):
    config = {
        "share-protected-areas-used": x[0],
        "share-pv-on-farmland": min(x[1], x[2]),
        "share-farmland-used": x[2],
        "share-forest-used-for-wind": x[3],
        "share-other-land-used": x[4],
        "share-offshore-used": x[5],
        "share-rooftops-used": x[6]
    }
    constrained_potentials = _constrain_potential(
        unconstrained_prefer_pv,
        unconstrained_prefer_wind,
        config
    ).sum(axis=1)
    return population[constrained_potentials > demand].sum() / population.sum()


def print_indices(S, problem, calc_second_order, file):
    # taken from the SAlib source code and modified to print to file
    # https://github.com/SALib/SALib/blob/3bc2ddfb50f091e5e5dd1ed4e7fae05853f150e5/SALib/analyze/sobol.py#L240
    # Output to console
    if not problem.get('groups'):
        title = 'Parameter'
        names = problem['names']
        D = problem['num_vars']
    else:
        title = 'Group'
        _, names = compute_groups_matrix(problem['groups'])
        D = len(names)

    print('%s S1 S1_conf ST ST_conf' % title, file=file)

    for j in range(D):
        print('%s %f %f %f %f' % (names[j], S['S1'][
            j], S['S1_conf'][j], S['ST'][j], S['ST_conf'][j]), file=file)

    if calc_second_order:
        print('\n%s_1 %s_2 S2 S2_conf' % (title, title), file=file)

        for j in range(D):
            for k in range(j + 1, D):
                print("%s %s %f %f" % (names[j], names[k],
                      S['S2'][j, k], S['S2_conf'][j, k]), file=file)


if __name__ == "__main__":
    sensititivity_analysis()
