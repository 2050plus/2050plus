# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2023 Climact for The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT

"""
Add tomorrow electrical load to network.
"""

import logging

import geopandas as gpd
import numpy as np
import pandas as pd
import pypsa
from vresutils import transfer as vtransfer

from _helpers import configure_logging

idx = pd.IndexSlice

logger = logging.getLogger(__name__)


def normed(s):
    return s / s.sum()


def attach_load_tomorrow(n, n_clustered, regions, load, nuts3_shapes, countries, bus_map_s, bus_map_s_c, scaling=1.0):
    """
    This function is very similar to `scripts.add_electricity.attach_load()`.
    See differentiation point.
    """
    substation_lv_i = n.buses.index[n.buses["substation_lv"]]
    regions = gpd.read_file(regions).set_index("name").reindex(substation_lv_i)
    opsd_load = pd.read_csv(load, index_col=0, parse_dates=True).filter(items=countries)

    logger.info(f"Load data scaled with scalling factor {scaling}.")
    opsd_load *= scaling

    nuts3 = gpd.read_file(nuts3_shapes).set_index("index")

    def upsample(cntry, group):
        l = opsd_load[cntry]
        if len(group) == 1:
            return pd.DataFrame({group.index[0]: l})
        else:
            nuts3_cntry = nuts3.loc[nuts3.country == cntry]
            transfer = vtransfer.Shapes2Shapes(
                group, nuts3_cntry.geometry, normed=False
            ).T.tocsr()
            gdp_n = pd.Series(
                transfer.dot(nuts3_cntry["gdp"].fillna(1.0).values), index=group.index
            )
            pop_n = pd.Series(
                transfer.dot(nuts3_cntry["pop"].fillna(1.0).values), index=group.index
            )

            # relative factors 0.6 and 0.4 have been determined from a linear
            # regression on the country to continent load data
            factors = normed(0.6 * normed(gdp_n) + 0.4 * normed(pop_n))
            return pd.DataFrame(
                factors.values * l.values[:, np.newaxis],
                index=l.index,
                columns=factors.index,
            )

    load = pd.concat(
        [
            upsample(cntry, group)
            for cntry, group in regions.geometry.groupby(regions.country)
        ],
        axis=1,
    )

    # Differentiated here
    load_nodes = substation_lv_i.astype(int).values
    nodes_s = (
        pd.read_csv(bus_map_s, index_col=0)
        .loc[load_nodes]
        .reset_index()
        .rename(columns={"0": "nodes_s", "index": "nodes"})
    )
    nodes_s_c = (
        pd.read_csv(bus_map_s_c, index_col=0)
        .loc[nodes_s["nodes_s"]]
        .reset_index()
        .astype(str)
    )
    nodes_map = (
        nodes_s_c
        .merge(nodes_s.astype(str), how="outer", left_on="Bus", right_on="nodes_s")
        .set_index("nodes")["busmap"]
        .to_dict()
    )

    load = load.rename(columns=nodes_map).groupby(by='Bus', axis=1).sum()
    load.index.names = ["snapshots"]

    n_clustered.loads_t.p_set = load


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake("add_electricity_tomorrow", simpl="",
                                   planning_horizons="2030", clusters="37", ll="v1.1", opts="")
    configure_logging(snakemake)

    n = pypsa.Network(snakemake.input.base_network)
    n_clustered = pypsa.Network(snakemake.input.network)

    if snakemake.config["enable"].get("modify_residual_load", False):
        attach_load_tomorrow(
            n,
            n_clustered,
            snakemake.input.regions,
            snakemake.input.load_future,
            snakemake.input.nuts3_shapes,
            snakemake.config["countries"],
            snakemake.input.bus_map_s,
            snakemake.input.bus_map_s_c,
            snakemake.config["load"]["scaling_factor"]
        )

    n_clustered.export_to_netcdf(snakemake.output[0])
