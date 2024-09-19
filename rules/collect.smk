# SPDX-FileCopyrightText: : 2023-2024 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT


localrules:
    all,
    cluster_networks,
    extra_components_networks,
    prepare_elec_networks,
    prepare_sector_networks,
    solve_elec_networks,
    solve_sector_networks,


rule cluster_networks:
    input: lambda w:
        expand(
            resources("networks/elec_s{simpl}_{clusters}.nc"),
            **config_provider("scenario")(w),
            run=config["run"]["name"],
        ),


rule extra_components_networks:
    input: lambda w:
        expand(
            resources("networks/elec_s{simpl}_{clusters}_ec.nc"),
            **config_provider("scenario")(w),
            run=config["run"]["name"],
        ),


rule prepare_elec_networks:
    input: lambda w:
        expand(
            resources("networks/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}.nc"),
            **config_provider("scenario")(w),
            run=config["run"]["name"],
        ),


rule prepare_sector_networks:
    input: lambda w:
        expand(
            RESULTS
            + "prenetworks/elec_s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}.nc",
            **config_provider("scenario")(w),
            run=config["run"]["name"],
        ),


rule solve_elec_networks:
    input: lambda w:
        expand(
            RESULTS + "networks/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}.nc",
            **config_provider("scenario")(w),
            run=config["run"]["name"],
        ),


rule solve_sector_networks:
    input: lambda w:
        expand(
            RESULTS
            + "postnetworks/elec_s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}.nc",
            **config_provider("scenario")(w),
            run=config["run"]["name"],
        ),


rule solve_sector_networks_perfect:
    input: lambda w:
        expand(
            RESULTS
            + "maps/elec_s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}-costs-all_{planning_horizons}.pdf",
            **config_provider("scenario")(w),
            run=config["run"]["name"],
        ),


rule validate_elec_networks:
    input: lambda w:
        expand(
            RESULTS
            + "figures/.statistics_plots_elec_s{simpl}_{clusters}_ec_l{ll}_{opts}",
            **config_provider("scenario")(w),
            run=config["run"]["name"],
        ),
        lambda w:
        expand(
            RESULTS
            + "figures/.validation_{kind}_plots_elec_s{simpl}_{clusters}_ec_l{ll}_{opts}",
            **config_provider("scenario")(w),
            run=config["run"]["name"],
            kind=["production", "prices", "cross_border"],
        ),
