# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2020-2024 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT
"""
Interpolate costs using 2030, 2040 and 2050 as input and produce costs for 2035 and 2045.
"""

import numpy as np
import pandas as pd
import yaml
from pathlib import Path


def load_scenarios(scenarios_file):
    with open(scenarios_file, 'r') as f:
        scenarios = yaml.safe_load(f)
    return scenarios


def interpolate(costs, prev_year, next_year, y, interval):

    costs[y] = costs[prev_year].copy()
    costs[y]['value'] = (costs[prev_year]['value'] +
                         costs[next_year]['value'])/2
    costs[y]['source'] = f"Linear interpolation between {str(prev_year)} and {str(next_year)}"
    costs[y]['further description'] = ''

    FOM_index = costs[y].index.get_level_values(1).isin(['FOM'])

    FOM_next = costs[next_year].reset_index(
        'parameter').query("parameter.isin(['FOM'])").value
    FOM_prev = costs[prev_year].reset_index(
        'parameter').query("parameter.isin(['FOM'])").value
    inv_next = costs[next_year].reset_index('parameter').query(
        "parameter.isin(['investment'])").value
    inv_prev = costs[prev_year].reset_index('parameter').query(
        "parameter.isin(['investment'])").value

    FOM = ((FOM_prev*inv_prev + FOM_next*inv_next) /
           (inv_next+inv_prev)).to_frame().dropna()
    FOM['parameter'] = 'FOM'
    FOM.reset_index(inplace=True)
    FOM.set_index(['technology', 'parameter'], inplace=True)

    costs[y].loc[FOM_index, 'value'] = FOM

    return costs[y]


def __main__():
    costs = {}
    scenarios = load_scenarios(Path('config/scenarios.veka.yaml'))

    boundaries = [2030, 2040, 2050]
    to_interpolate = [2035, 2045]
    assert (len(boundaries) == (len(to_interpolate)+1)
            ), 'Dimension mismatch between years to interpolate and boundary years'
    interval = np.diff(boundaries)/2

    # Get costs for planning horizons
    for scen in list(scenarios.keys()):
        costs[scen] = {}
        for y in boundaries:
            costs[scen][y] = pd.read_csv(
                Path(f"data/costs/{scen}/costs_{y}.csv"), index_col=[0, 1])

        for count, y in enumerate(to_interpolate):
            prev_year = boundaries[count]
            next_year = boundaries[count + 1]
            interpolate(costs[scen], prev_year, next_year, y, interval[count])

            costs[scen][y].to_csv(Path(f"data/costs/{scen}/costs_{y}.csv"))

__main__()