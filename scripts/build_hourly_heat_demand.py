# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2020-2023 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT
"""
Build hourly heat demand time series from daily ones.
"""

import pandas as pd
import xarray as xr
from _helpers import generate_periodic_profiles, update_config_with_sector_opts
from itertools import product



if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake(
            "build_heat_demands",
            simpl="",
            clusters=48,
        )

    snapshots = pd.date_range(freq="h", **snakemake.params.snapshots)

    daily_space_heat_demand = (
        xr.open_dataarray(snakemake.input.heat_demand)
        .to_pandas()
        .reindex(index=snapshots, method="ffill")
    )

    intraday_profiles = pd.read_csv(snakemake.input.heat_profile, index_col=0)

    sectors = ["residential", "services"]
    uses = ["water", "space"]

    heat_demand = {}
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
            heat_demand[f"{sector} {use}"] = daily_space_heat_demand * intraday_year_profile
        else:
            heat_demand[f"{sector} {use}"] = intraday_year_profile

    heat_demand = pd.concat(heat_demand,
                            axis=1,
                            names = ["sector use", "node"])

    heat_demand.index.name="snapshots"

    print(heat_demand)

    print(heat_demand.stack())

    ds = heat_demand.stack().to_xarray()#xr.Dataset.from_dataframe(heat_demand)

    print(ds)

    ds.to_netcdf(snakemake.output.heat_demand)
