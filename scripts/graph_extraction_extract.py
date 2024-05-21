# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2024 Climact for The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT

"""
Create data ready to present (extract)
"""
import logging
from pathlib import Path

import numpy as np
import pypsa
from scripts.graph_extraction_utils import bus_mapper
from scripts.make_summary import assign_carriers
from scripts.make_summary import assign_locations

logger = logging.getLogger(__name__)


def assign_countries(n):
    n.buses = (
        n.buses.merge(
            n.buses[(n.buses.location != "EU") & (n.buses.carrier == "AC")]["country"],
            how="left", left_on="location", right_index=True, suffixes=("_old", '')
        )
        .fillna({"country": "EU"})
        .drop(columns="country_old")
    )
    return


def assign_coordinates(n):
    n.buses.loc[n.buses.query('x == 0 or y == 0').index, ['x', 'y']] = ''
    n.buses = (
        n.buses.merge(
            n.buses[(n.buses.location != "EU") & (n.buses.carrier == "AC")][["x", "y"]],
            how="left", left_on="location", right_index=True, suffixes=("_old", '')
        )
        .fillna({"x": 0, "y": 0})
        .drop(columns=["x_old", "y_old"])
    )
    return


def searcher(x, carrier):
    if carrier in x.to_list():
        return str(x.to_list().index(carrier))
    else:
        return np.nan


def change_p_nom_opt_carrier(n, carriers=['AC'], temporal=True):
    """
    This function expresses for each asset p_nom_opt (whose carrier is not always the same)
    in function of the given carrier, considering the efficiency
    E.g.: having p_nom_opt = 333MW for a nuclear powerplant is equal to say that
            p_carrier_nom_opt = 100MW_e, if the carrier is 'AC', as the efficiency
            from p_nom_opt to AC is 0.333 for this technology

    Currently, carriers should be limited to ['AC'] as no strict testing was made with other carriers.
    In a future development, this function should be able to tackle several carriers

    Parameters
    ----------
    n : pypsa.network

    carriers : List, optional
        List of carriers to use. Must be limited to ['AC'] for the moment, as bugs might arise. The default is ['AC'].

    Returns
    -------
    None.

    """

    # Beware this also has extraneous locations for country (e.g. biomass) or continent-wide (e.g. fossil gas/oil) stuff
    li = n.links
    li_t = n.links_t
    li["efficiency0"] = 1
    li["p_carrier_nom_opt"] = li.p_nom_opt.copy()
    li_t["p_carrier_nom_opt"] = li_t.p0.copy()
    efficiency_map = li[[c for c in li.columns if "efficiency" in c]].rename(columns={"efficiency": "efficiency1"})
    buses_links = [c for c in li.columns if "bus" in c]
    carrier_map = li[buses_links].map(lambda x: bus_mapper(x, n, column="carrier"))

    for carrier in carriers:
        index_map = carrier_map.apply(lambda x: searcher(x, carrier), axis=1).dropna()

        efficiency_map = efficiency_map.loc[index_map.index]
        efficiency_map = efficiency_map.apply(lambda x: x / x[f"efficiency{index_map.loc[x.name]}"], axis=1)
        li.loc[efficiency_map.index, "p_carrier_nom_opt"] = li.loc[efficiency_map.index, "p_nom_opt"] / efficiency_map[
            "efficiency0"]
        if len(li_t.p0.columns):
            li_t.p_carrier_nom_opt.loc[:, efficiency_map.index] = li_t.p0.loc[:, efficiency_map.index] / efficiency_map[
                "efficiency0"]

    return


def extract_data(config):
    logger.info("Extracting data")

    n = {}
    for y in config["scenario"]["planning_horizons"]:
        run_name = Path(config["path"]["results_path"], "postnetworks", config["n_name"] + f"{y}.nc")
        n[y] = pypsa.Network(run_name)
        assign_carriers(n[y])
        assign_locations(n[y])
        assign_countries(n[y])
        assign_coordinates(n[y])
        change_p_nom_opt_carrier(n[y])

    # get historical capacities
    n_bf = pypsa.Network(Path(config["path"]["results_path"], "prenetworks-brownfield", config["n_name"] + f"{2030}.nc"))
    assign_countries(n_bf)
    assign_carriers(n_bf)
    assign_locations(n_bf)
    assign_coordinates(n_bf)

    year_hist = 2026
    n_bf.generators = n_bf.generators.query(f'build_year < {year_hist}')
    n_bf.links = n_bf.links.query(f'build_year < {year_hist}')
    n_bf.storage_units = n_bf.storage_units.query(f'build_year < {year_hist}')
    n_bf.generators.p_nom_opt = n_bf.generators.p_nom
    n_bf.links.p_nom_opt = n_bf.links.p_nom
    n_bf.storage_units.p_nom_opt = n_bf.storage_units.p_nom
    n_bf.lines.carrier = "AC"
    n_bf.lines.s_nom_opt = n_bf.lines.s_nom_min
    n_bf.links.loc[n_bf.links.carrier.isin(["DC"]), 'p_nom_opt'] = n_bf.links.loc[
        n_bf.links.carrier.isin(["DC"]), 'p_nom_min']
    change_p_nom_opt_carrier(n_bf)

    n_ext = n.copy()
    n_ext['hist'] = n_bf.copy()

    return n, n_ext
