# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2023 Climact for The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT

"""
Creates custom residual demand profile based on desegregated annual values and hourly reference profiles.
"""

import logging

import pandas as pd

logger = logging.getLogger(__name__)
from _helpers import configure_logging


def apply_hourly_load_profiles(load_hourly, load_annual, profiles, heat_map, transport_map):
    """
    Apply reference hourly load profiles on sub sectors of annual load.
    We apply detailed profiles on sub sectors of annual load.
    We remove those profiles from annual hourly load to get residual load.
    1. Heat: Compute load profile for each country and sub sector as fraction of annual demand
    2. Transport: Compute load profile for each country and sub sector as fraction of annual demand
    3. Industry: Compute total profile for each country
    4. Supply: We use a rule of thumb on total to go from annual to hourly values. Then, assume transmission losses as load.
    5. Compute residual and non-residual load for each country
    """
    countries = load_hourly.columns
    snapshots = load_hourly.index.rename("utc timestamp")
    nyears = len(snapshots) / 8760
    if nyears != 1:
        logging.warning(f"Snapshots used to compute profiles are not on one year, which can be hazardous.")

    load_annual = load_annual.reindex(index=snapshots, method="ffill")

    load_residual = load_hourly.copy()
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

        # Supply
        tr_losses = 0.05
        load_annual_hourly = load_hourly[c] / load_hourly[c].sum() * load_annual[c + "_tot"].mean() * nyears
        supply = load_annual_hourly.values * tr_losses

        # Non-residual load
        non_residual = supply + industry + heat + transport
        non_residual = non_residual.fillna(0)
        non_residual.index = snapshots

        # Residual load
        logger.info(
            "Total non-residual consumption for {} is {:.2f}TWh (for reference hourly load used of {:.2f}TWh)"
            .format(c, non_residual.sum() / 1e6, load_hourly[c].sum() / 1e6))
        load_residual[c] -= non_residual

        logging.debug(f"Reference hourly load from {c}: {load_hourly[c].sum().sum() / 1e6:.2f}TWh")
        logging.debug(f"Residual load from {c}: {load_residual[c].sum().sum() / 1e6:.2f}TWh")
        logging.debug(f"Supply from {c}: {supply.sum() / 1e6:.2f}TWh")
        logging.debug(f"Industry from {c}: {industry.sum() / 1e6:.2f}TWh")
        logging.debug(f"Heat from {c}: {heat.sum().sum() / 1e6:.2f}TWh")
        logging.debug(f"Transport from {c}: {transport.sum().sum() / 1e6:.2f}TWh")
        logging.debug(f"Total from {c}: {non_residual.sum().sum() / 1e6:.2f}TWh")

    return load_residual


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake("build_residual_load_profile", simpl="", planning_horizons="2050",
                                   countries=["BE"], )

    configure_logging(snakemake)

    countries = snakemake.config["countries"]
    load_annual = (
            pd.read_csv(snakemake.input.load_annual, parse_dates=True, index_col=0)
            .filter(regex='|'.join([f"^{c}.*" for c in countries]))
            * 1e6
    )
    load_hourly = (
        pd.read_csv(snakemake.input.load_hourly, parse_dates=True, index_col="utc_timestamp")
        .filter(items=countries)
    )
    profiles = pd.read_csv(snakemake.input.profiles)
    heat_map = pd.read_csv(snakemake.input.heat_map).to_dict(orient="index")[0]
    transport_map = pd.read_csv(snakemake.input.transport_map).to_dict(orient="index")[0]

    # ToDo What if snapshots are considered through multiple years
    logger.info(
        f"Building residual load from {len(snakemake.config['countries'])} countries using year "
        f"{load_annual.index[0].year} ({', '.join(snakemake.config['countries'])})")

    load_residual = apply_hourly_load_profiles(load_hourly, load_annual, profiles, heat_map, transport_map)

    residual_profile = load_residual / load_residual.sum() * len(load_residual) / 8760
    residual_profile.to_csv(snakemake.output.res_load_profile)
