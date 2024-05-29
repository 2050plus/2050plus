# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2020-2024 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT
"""
Build industrial energy demand per model region taking into account data from VLAIO. This is a custom rule
designed for VEKA 2050+ project.
"""
import logging

import numpy as np
import pandas as pd

from _helpers import set_scenario_config


def build_flanders_energy_demand(demand, energy_fl_ets, energy_fl_non_ets, map):
    """
    Build Flanders energy demand using VLAIO data and energy carrier mapping. Steel is not included
    as it is defined by PyPSA.
    """

    vlaio_industries = [
        "ceramics",
        "chemicals-ammonia",
        "chemicals-chlorine",
        "chemicals-HVC",
        "chemicals-others",
        "food",
        "glass",
        "non-ferrous",
        "paper",
        # "steel",  # Will always be configured through PyPSA
    ]

    vlaio_carriers = [
        "chemicals-ethanol",
        "chemicals-ethanol-import",
        "chemicals-methanol",
        "electricity",
        "electricity-for-CC",
        "ethane-import",
        "gaseous-fuel-bio",
        "gaseous-fuel-fossil",
        "gaseous-fuel-syn",
        "gaseous-fuel-syn-Cbased",
        "gaseous-fuel-syn-nonCbased",
        "heat-medium",
        "liquid-fuel-bio",
        "liquid-fuel-fossil",
        "liquid-fuel-syn",
        "liquid-fuel-syn-Cbased",
        "liquid-fuel-syn-nonCbased",
        "refineries-LPG",
        "refineries-naphtha",
        "refineries-naphtha-import",
        "solid-fuel-bio",
        "solid-fuel-fossil",
    ]

    demand_fl_ets = (
        energy_fl_ets.loc[int(snakemake.wildcards.planning_horizons)].copy()
        .query("industry.isin(@vlaio_industries) and `energy carrier`.isin(@vlaio_carriers)")
        .merge(map, left_on="energy carrier", right_index=True, how="left")
        .assign(data_pypsa=lambda x: x["data"] * x["map_factor"])
        .pivot_table(columns="energy carrier pypsa", values="data_pypsa", index="industry", aggfunc="sum")
        .replace(0, np.nan).dropna(how="all", axis=1).dropna(how="all", axis=0)
        .fillna(0)
    )

    demand_fl_non_ets = (
        energy_fl_non_ets.loc[int(snakemake.wildcards.planning_horizons)]
        .to_frame("non-ETS")
        .T
        .replace(0, np.nan).dropna(how="all", axis=1).dropna(how="all", axis=0)
    )

    # Add steel demand from PyPSA
    steel_pypsa = [
        "DRI CH4 + Electric arc",
        "DRI H2 + Electric arc",
        "Electric arc",
        "Integrated steelworks",
    ]
    demand_fl_steel = (
        demand.query("`TWh/a (MtCO2/a)` == 'BE1 0' and industry.isin(@steel_pypsa)")
        .loc["BE1 0"]
        .replace(0, np.nan).dropna(how="all", axis=1).dropna(how="all", axis=0)
    )

    # Add current electricity reference
    demand_fl_current_el = (
        demand.loc["BE1 0", ["current electricity"]]
        .replace(0, np.nan).dropna(how="all", axis=1).dropna(how="all", axis=0)
    )

    demand_fl = (
        pd.concat([
            demand_fl_ets,
            demand_fl_non_ets,
            demand_fl_steel,
            demand_fl_current_el,
        ])
        .fillna(0)
        .reset_index()
        .assign(node="BE1 0")
        .rename(columns={"node": "TWh/a (MtCO2/a)", "index": "industry"})
        .set_index(["TWh/a (MtCO2/a)", "industry"])
    )

    return demand_fl


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake(
            "build_industrial_energy_demand_per_node_for_vlaio",
            simpl="181",
            clusters="39m",
            planning_horizons=2030,
            run="electrification",
            configfiles="config/config.veka.yaml",
        )
    set_scenario_config(snakemake)

    # industrial energy demand per node
    fn = snakemake.input.industrial_energy_demand_per_node_ind
    demand = pd.read_csv(fn, header=0, index_col=[0, 1])

    scenario = snakemake.config["industry"].get("vlaio_scenario", "MIX CENTRAL")
    demand_vlaio = demand.copy()
    if scenario:
        logging.info("Overriding industrial energy demand based on VLAIO study")

        # import vlaio's Flanders ETS energy demand
        fn = snakemake.input.vlaio_energy_demand_fl_ets
        energy_fl_ets = (
            pd.read_csv(fn, header=0, index_col=0)
            .query("scenario == @scenario")
        )

        # import vlaio's Flanders non-ETS energy demand
        fn = snakemake.input.vlaio_energy_demand_fl_non_ets
        energy_fl_non_ets = pd.read_csv(fn, header=0, index_col=0)

        # import vlaio's vector mapping
        fn = snakemake.input.vlaio_vector_mapping
        map = (
            pd.read_csv(fn, header=0, index_col=0)
            .replace(0, np.nan)
            .melt(ignore_index=False, var_name="energy carrier pypsa", value_name="map_factor")
            .dropna()
        )

        demand_fl = build_flanders_energy_demand(demand, energy_fl_ets, energy_fl_non_ets, map)

        demand_vlaio = pd.concat([demand_vlaio.drop("BE1 0"), demand_fl])

    fn = snakemake.output.industrial_energy_demand_per_node_for_vlaio
    demand_vlaio.groupby(by="TWh/a (MtCO2/a)").sum().to_csv(fn, float_format="%.2f")
    fn = snakemake.output.industrial_energy_demand_per_node_ind_for_vlaio
    demand_vlaio.to_csv(fn, float_format="%.2f")
