# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2023 Climact for The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT

"""
Creates futur load profiles.
"""

import logging

import pandas as pd

logger = logging.getLogger(__name__)
from _helpers import configure_logging


def apply_profiles_tomorrow(load_annual, countries, profiles, heat_map, transport_map, snapshots):
    """
    Apply defined profiles to futur sub sector annual load.
    1. Heat: Compute load profile for each country and sub sector as fraction of annual demand
    2. Transport: Compute load profile for each country and sub sector as fraction of annual demand
    3. Industry: Compute total profile for each country
    4. Supply: We use a rule of thumb on total to go from annual to hourly values. Then, assume transmission losses as load.
    5. Compute residual load for each country
    """
    load_annual = load_annual.reindex(index=snapshots, method="ffill")

    load_futur = pd.DataFrame([])
    for c in countries:
        # Heat
        heat = {}
        for he in heat_map.values():
            heat[c + '_' + he] = load_annual[c + '_' + he].values * profiles[c + '_' + he]
        heat = pd.concat(heat, axis=1).sum(axis=1)

        # Transport
        transport = {}
        for tr in transport_map.keys():
            transport[c + '_' + tr] = load_annual[c + '_' + tr].values * profiles[c + '_' + tr]
        transport = pd.concat(transport, axis=1).sum(axis=1)

        # Industry
        industry = load_annual[c + "_IN_tot"].values * profiles[c + "_IN_tot"]

        # Residual load
        # ToDo Switch utc_timestamp to snapshots
        residual_profile = pd.read_csv(snakemake.input.res_load_profile, index_col=['utc_timestamp'])
        residual = (
                (load_annual[c + "_tot"] -
                 (
                     load_annual
                     .filter(regex=f"^{c}_.*")
                     .filter(regex=f"^(?!{c}_tot|{c}_TR_tot).*$")
                     .sum(axis=1)
                 )
                 ).mean()
                * residual_profile[c].values
        )
        # ToDo Set tr_losses as option in yaml
        tr_losses = 0.05
        supply = (industry + heat + transport + residual) * tr_losses

        # Futur load
        load_futur[c] = industry + heat + transport + supply + residual
        logger.info(f"Build total load for {c} is {load_futur[c].sum() / 1e6:.2f} TWh")

    return load_futur.set_index(snapshots)


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake("build_future_load", simpl="",
                                   planning_horizons="2030", clusters="37")

    configure_logging(snakemake)

    horizon = snakemake.wildcards.planning_horizons
    logger.info(f"Building futur load for year {horizon}")

    # ToDo What if leap year (e.g.: 2040)
    # ToDo Adjust weekly pattern to new year
    snapshots = pd.date_range(freq="h", **snakemake.config["snapshots"])
    snapshots = pd.DatetimeIndex([i.replace(year=int(horizon)) for i in snapshots.to_list()])

    # FixMe Rule retrieve_futur_load should define a name for year index column
    load_annual_futur = pd.read_csv(snakemake.input.load_annual, delimiter=',', index_col=0, parse_dates=True)
    load_annual_futur.index.names = ["year"]

    profiles = pd.read_csv(snakemake.input.profiles)
    heat_map = pd.read_csv(snakemake.input.heat_map).to_dict(orient="index")[0]
    transport_map = pd.read_csv(snakemake.input.transport_map).to_dict(orient="index")[0]

    load_future = apply_profiles_tomorrow(load_annual_futur, snakemake.config["countries"], profiles, heat_map,
                                          transport_map, snapshots)

    # ToDo check wildcards pour la création de wild card en output
    load_future.to_csv(snakemake.output.load_future)
