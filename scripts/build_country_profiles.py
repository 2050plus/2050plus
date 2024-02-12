# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2023 Climact for The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT

"""
Creates sub sectors profiles based on desegregated annual values and hourly reference profiles.
"""

import logging
import re
import datetime
from itertools import product
from os.path import exists

import numpy as np
import pandas as pd
import xarray as xr

from _helpers import configure_logging
from _helpers import generate_periodic_profiles


def build_heat_profiles(snapshots):
    """
    Build heat demand profiles for all countries at once
    This duplicates part of the code used in `build_heat_demand` defined in `prepare_sector_network`
    """
    heat_map = {"residential water": "HE_res_wat",
                "residential space": "HE_res_spa",
                "services water": "HE_ter_wat",
                "services space": "HE_ter_spa"}

    nyears = len(snapshots) / 8760

    # ToDo Avoid Duplicated code from `prepare_sector_network`
    # copy forward the daily average heat demand into each hour, so it can be multiplied by the intraday profile
    daily_space_heat_demand = (
        xr.open_dataarray(snakemake.input.heat_demand_total)
        .to_pandas()
        .reindex(index=snapshots, method="ffill")
    )

    intraday_profiles = pd.read_csv(snakemake.input.heat_profile, index_col=0)

    sectors = ["residential", "services"]
    uses = ["water", "space"]

    heat_profiles = {}
    for sector, use in product(sectors, uses):
        weekday = list(intraday_profiles[f"{sector} {use} weekday"])
        weekend = list(intraday_profiles[f"{sector} {use} weekend"])
        weekly_profile = weekday * 5 + weekend * 2
        intraday_year_profile = generate_periodic_profiles(
            daily_space_heat_demand.index.tz_localize("UTC"),
            nodes=daily_space_heat_demand.columns,
            weekly_profile=weekly_profile,
        )

        if use == "space":
            heat_demand_shape = daily_space_heat_demand * intraday_year_profile
        else:
            heat_demand_shape = intraday_year_profile

        heat_profiles[heat_map[f"{sector} {use}"]] = heat_demand_shape / heat_demand_shape.sum() * nyears

    heat_profiles = pd.concat(heat_profiles, axis=0)
    heat_profiles.index.names = ["uses", "utc timestamp"]
    heat_profiles = heat_profiles.pivot_table(index="utc timestamp", columns="uses")

    return heat_profiles, heat_map


def build_transport_profiles(snapshots, nodes=[]):
    """
    Build transport demand profiles for all countries at once
    """
    nyears = len(snapshots) / 8760
    mobility_fn = snakemake.input.mobility_fn
    transport_map = {"TR_cars": mobility_fn + "/Pkw__count",
                     "TR_bus": mobility_fn + "/Bus__count",
                     "TR_ldv": mobility_fn + "/Lfw__count",
                     "TR_hdv": mobility_fn + "/HDV__count",
                     "TR_rail": mobility_fn + "/KFZ__count"}  # Take total road number

    transport_profiles = {}
    for sector, traffic_fn in transport_map.items():
        if (sector == "TR_hdv") & (not exists(traffic_fn)):  # produce HDV profile
            loa = pd.read_csv(mobility_fn + "/LoA__count", skiprows=2)
            lzg = pd.read_csv(mobility_fn + "/Lzg__count", skiprows=2)
            sat = pd.read_csv(mobility_fn + "/Sat__count", skiprows=2)
            HDV_count = lzg + sat + loa
            with open(traffic_fn, 'w') as fd:
                fd.write(
                    f"File generated for type: HDV_\n"
                    f"Time of generation:{datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
                HDV_count.to_csv(fd, index=False)
            traffic = HDV_count["count"]
        else:
            traffic = pd.read_csv(traffic_fn, skiprows=2, usecols=["count"]).squeeze("columns")

        transport_shape = generate_periodic_profiles(dt_index=snapshots.tz_localize("UTC"),
                                                     nodes=nodes,
                                                     weekly_profile=traffic.values)

        transport_profiles[sector] = transport_shape / transport_shape.sum() * nyears

    transport_profiles = pd.concat(transport_profiles, axis=0)
    transport_profiles.index.names = ["uses", "utc timestamp"]
    transport_profiles = transport_profiles.pivot_table(index="utc timestamp", columns="uses")

    return transport_profiles, transport_map


def build_country_profiles(heat_profiles, transport_profiles, snapshots):
    """
    Build sectoral profiles for each country
    """
    sectors = set([re.sub(r"^[A-Z]{2}_(.*)$", r"\1", i) for i in load_annual.columns]) - set(["tot", "TR_tot"])
    clustered_pop = pd.read_csv(snakemake.input.clustered_pop_layout).set_index("name")

    profiles = []
    for country in snakemake.config["countries"]:
        # Use population fraction to weight profiles in each country
        nodes = list(set([i for i, j in heat_profiles.columns if country in i]))
        nodes_fraction = clustered_pop.loc[nodes, "fraction"].to_dict()

        heat_profiles_country = (heat_profiles[nodes]
                                 .apply(lambda x: nodes_fraction[x.name[0]] * x)
                                 .groupby(level=["uses"], axis=1)
                                 .sum()
                                 )
        transport_profiles_country = (transport_profiles[nodes]
                                      .apply(lambda x: nodes_fraction[x.name[0]] * x)
                                      .groupby(level=["uses"], axis=1)
                                      .sum()
                                      )
        industry_supply_profile = np.ones(len(snapshots)) / 8760

        for sector in sectors:
            if "IN" in sector or "SU" in sector:
                profiles.append(
                    pd.DataFrame(industry_supply_profile, columns=[country + '_' + sector], index=snapshots))
                logging.debug(f"{'Industry' if 'IND' in sector else 'Supply'} ({sector}) profile sum for {country}"
                              f": {industry_supply_profile.sum():.2f}")
            elif "HE_" in sector:
                heat_profile = (
                    heat_profiles_country[[sector]]
                    .rename(columns={sector: country + '_' + sector})
                )
                profiles.append(heat_profile)
                logging.debug(f"Heat ({sector}) profile sum for {country}: {heat_profile.sum()[0]:.2f}")

            elif "TR_" in sector:
                transport_profile = (
                    transport_profiles_country[[sector]]
                    .rename(columns={sector: country + '_' + sector})
                )
                profiles.append(transport_profile)
                logging.debug(f"Transport ({sector}) profile sum for {country}: {transport_profile.sum()[0]:.2f}")

    return pd.concat(profiles, axis=1)


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake("build_country_profiles", simpl="",
                                   planning_horizons="2030", clusters="37")

    configure_logging(snakemake)

    load_annual = pd.read_csv(snakemake.input.load_annual, parse_dates=True, index_col=0)
    sectors = set([re.sub(r"^[A-Z]{2}_(.*)$", r"\1", i) for i in load_annual.columns]) - set(["tot", "TR_tot"])
    load_hourly = pd.read_csv(snakemake.input.load_hourly, parse_dates=True, index_col=0)
    snapshots = load_hourly.index.rename("utc timestamp")

    heat_profiles, heat_map = build_heat_profiles(snapshots)
    transport_profiles, transport_map = build_transport_profiles(snapshots,
                                                                 nodes=heat_profiles.columns.get_level_values(
                                                                     0).drop_duplicates())

    profiles = build_country_profiles(heat_profiles, transport_profiles, snapshots)
    profiles.to_csv(snakemake.output.profiles)

    pd.DataFrame.from_dict([heat_map]).to_csv(snakemake.output.heat_map, index=False)
    pd.DataFrame.from_dict([transport_map]).to_csv(snakemake.output.transport_map, index=False)
