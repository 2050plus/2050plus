# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2024 Climact for The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT

"""
Create data ready to present (load to streamlit)
"""

import logging
import re
from pathlib import Path

import pandas as pd

from scripts.graph_extraction_utils import CLIP_VALUE_TWH
from scripts.graph_extraction_utils import NICE_RENAMER
from scripts.graph_extraction_utils import _load_supply_energy
from scripts.graph_extraction_utils import query_imp_exp

COST_SEGMENTS = {'prod': 'Energy production', 'sto': 'Storage', 'tran': 'Transmission', 'distr': 'Distribution',
                 'net_imp': "Net_Imports"}

logger = logging.getLogger(__name__)


def load_supply_energy_df(config, load=True):
    """
    Allow to split _load_supply_energy into its carriers for it
    to be written tab by tab

    """
    dfl = []
    supply_energy_carrier = (
        _load_supply_energy(config, load=load)
        .carrier
        .replace({"AC": "electricity", "low voltage": "electricity"})
        .unique()
    )

    for ca in supply_energy_carrier:

        # Todo : should we move the HV/LV/imports/exports to the calling function to keep this function read only (no modifications) ?
        if ca == "electricity":
            df_ac = _load_supply_energy(config, load=load, carriers="AC", aggregate=False)
            df_low = _load_supply_energy(config, load=load, carriers="low voltage", aggregate=False)
            df = pd.concat([df_ac, df_low])
            del df["carrier"]
            df = df.groupby(by=["sector", "node"]).sum().reset_index()
            df.drop(df.query('sector in ["V2G", "Battery charging", "Hydroelectricity"]').index, inplace=True)
            df["carrier"] = ca
            dfl.append(df)
        else:
            dfl.append(_load_supply_energy(config, load=load, carriers=ca, aggregate=False))

    df = pd.concat(dfl).replace(NICE_RENAMER)
    return df


def load_res_potentials(config):
    df = pd.read_csv(Path(config["path"]["csvs"], "res_potentials.csv"), header=0)
    return (
        df.loc[:, ['carrier', 'region', config['years_str'][-1]]]
        .rename(columns={'region': 'country', config['years_str'][-1]: 'Potential [GW]'})
        .pivot(columns='carrier', index='country', values='Potential [GW]')
        .fillna(0)
        .reindex(columns=['hydro', 'ror', 'offwind', 'onwind', 'solar'])
        .reset_index()
        .replace(NICE_RENAMER)
    )


def load_imports_exports(config):
    """
    This function loads the imports and exports for all countries, carriers and years
    considered during the runs. The table loaded is imports_exports, as it is only filtering that is done
    at streamlit level.

    """

    return pd.read_csv(Path(config["path"]["csvs"], "imports_exports.csv")).replace(NICE_RENAMER)


def _load_load_temporal(config, countries=None):
    load_raw = _load_supply_energy(config, load=True, aggregate=True, temporal=True, countries=countries)
    load = {}
    for y in [c for c in load_raw.columns if re.match(r"[0-9]{4}", c)]:
        load_i = load_raw.pivot_table(values=y, index=["carrier", "sector"], columns="snapshot")
        load_i.columns = pd.to_datetime(pd.DatetimeIndex(load_i.columns, name='snapshots').strftime(f'{y}-%m-%d-%H'))
        load_i = load_i.loc[load_i.sum(axis=1) / 1e3 > CLIP_VALUE_TWH, :]
        load_i = load_i.loc[(load_i.std(axis=1) / load_i.mean(axis=1)).sort_values().index]
        load[y] = load_i.reset_index().replace(NICE_RENAMER).T.reset_index()
    return load


def load_load_temporal(config):
    return _load_load_temporal(config)


def load_load_temporal_be(config):
    return _load_load_temporal(config, "BE")


def load_load_temporal_fl(config):
    return _load_load_temporal(config, "FL")


def _load_supply_temporal(config, countries=None):
    supply_raw = _load_supply_energy(config, load=False, aggregate=True, temporal=True, countries=countries)
    supply = {}
    for y in [c for c in supply_raw.columns if re.match(r"[0-9]{4}", c)]:
        supply_i = supply_raw.pivot_table(values=y, index=["carrier", "sector"], columns="snapshot")
        supply_i.columns = pd.to_datetime(
            pd.DatetimeIndex(supply_i.columns, name='snapshots').strftime(f'{y}-%m-%d-%H'))
        supply_i = supply_i.loc[supply_i.sum(axis=1) / 1e3 > CLIP_VALUE_TWH, :]
        supply_i = supply_i.loc[(supply_i.std(axis=1) / supply_i.mean(axis=1)).sort_values().index]
        supply[y] = supply_i.reset_index().replace(NICE_RENAMER).T.reset_index()
    return supply


def load_supply_temporal(config):
    return _load_supply_temporal(config)


def load_supply_temporal_be(config):
    return _load_supply_temporal(config, "BE")


def load_supply_temporal_fl(config):
    return _load_supply_temporal(config, "FL")


def load_res_temporal(config):
    res_raw = pd.read_csv(Path(config["path"]["csvs"], "temporal_res_supply.csv"), header=0).replace(NICE_RENAMER)
    res = {}
    for y in res_raw["year"].unique():
        res_i = res_raw.query("year==@y").drop("year", axis=1)
        res_i = res_i.set_index(['country', 'carrier'])
        res_i = res_i.astype(float).abs()
        res_i = res_i.loc[res_i.sum(axis=1) * 8760 / res_i.shape[0] > CLIP_VALUE_TWH].reset_index()
        res_i.index = pd.to_datetime(pd.DatetimeIndex(res_i.index, name='snapshots').strftime(f'{y}-%m-%d-%H'))
        res[str(y)] = res_i
    return res


def load_power_capacities(config):
    return pd.read_csv(Path(config["path"]["csvs"], "power_production_countries.csv"))


def load_balancing_capacities(config):
    return pd.read_csv(Path(config["path"]["csvs"], "balancing_capacities_countries.csv"))


def load_balancing_supply(config):
    return pd.read_csv(Path(config["path"]["csvs"], "balancing_supply_countries.csv"))


# generic function for calling costs
def _load_costs_year_segment(config, year=None, _countries=None, cost_segment=None):
    """
    Return the costs per segment for a given year or per year for a given segment,
    considering a subset of countries to consider
    Parameters
    ----------
    year : TYPE, optional
        DESCRIPTION. The default is None.
    _countries : TYPE, optional
        DESCRIPTION. The default is None.
    cost_segment : TYPE, optional
        DESCRIPTION. The default is None.

    Returns
    -------
    df : TYPE
        DESCRIPTION.

    """
    df = pd.read_csv(Path(config["path"]["csvs"], "costs_countries.csv"), header=0)
    prices = pd.read_csv(Path(config["path"]["csvs"], 'marginal_prices_countries.csv'), header=0)

    if _countries:
        df = df.query("country in @_countries")
        countries = list(set(_countries).intersection(set(df.country.unique())))
        if "BE" in countries:
            countries = [c for c in countries if c != "BE"] + ["BX", "FL", "WL"]
        else:
            countries = countries
    else:
        countries = None

    df = (
        df
        .groupby(["cost_segment", "cost"]).sum(numeric_only=True)
        .reset_index()
    )

    if cost_segment:
        net_cost = pd.DataFrame([], columns=config["imp_exp_carriers"], index=config["scenario"]["planning_horizons"])
        if cost_segment != "Net_Imports":
            df = df.query('cost_segment in @cost_segment')

        if cost_segment == "Energy production" or cost_segment == "Net_Imports":
            for y in config["scenario"]["planning_horizons"]:
                for ca in config["imp_exp_carriers"]:
                    imp = _load_imp_exp(config, export=False, countries=countries, carriers=ca, years=[y]).set_index(
                        'countries') * 1e6  # MWh
                    exp = _load_imp_exp(config, export=True, countries=countries, carriers=ca, years=[y]).set_index(
                        'countries') * 1e6  # MWh
                    price_ca = prices.query("carrier == @ca").set_index('countries').loc[:, str(y)]  # €/MWh
                    net_cost.loc[y, ca] = 0
                    if len(imp) > 0:
                        net_cost.loc[y, ca] += price_ca.loc[imp.index].dot(imp).sum()
                    if len(exp) > 0:
                        net_cost.loc[y, ca] -= (price_ca.loc[countries].mean() * exp).values.sum()
            df.loc[df.query("'fuel' in cost").index, config["years_str"]] += net_cost.sum(axis=1).astype(float).values
        if cost_segment == "Net_Imports":
            df = net_cost.reset_index()
            df_melted = pd.melt(df, id_vars=['index'], var_name='cost', value_name='value')
            df_pivoted = df_melted.pivot(index='cost', columns='index', values='value')
            df_pivoted.reset_index(inplace=True)
            df_pivoted.columns=df_pivoted.columns.astype(str)
            df = df_pivoted.copy()
            df["cost_segment"] = cost_segment

        df_tot = (df.groupby("cost_segment")
                  .sum()
                  .reset_index()
                  .assign(cost_segment=cost_segment,cost="total"))
        df = pd.concat([df,df_tot])
    else:
        df = (
            df.pivot(columns="cost", values=year, index="cost_segment")
            .fillna(0)
            .reset_index()
        )
    return df


def _load_costs(config, per_segment=False, per_year=False):
    dico = {}
    for co_name, subset in config["countries"].items():
        if per_segment:
            for seg_name, seg in COST_SEGMENTS.items():
                dico[f"{seg_name}_{co_name}"] = _load_costs_year_segment(config, _countries=subset, cost_segment=seg).rename(columns={"cost": "cost/carrier"})              
        elif per_year:
            for y in config["years_str"]:
                dico[f"{y}_{co_name}"] = _load_costs_year_segment(config, _countries=subset, year=y)
        else:
            logging.warning("Unkown configuration to load costs.")
    dico = pd.concat(dico.values(), keys=dico.keys()).droplevel(1).reset_index().rename(columns={"index": "config"})
    return dico


def load_costs_segments(config):
    return _load_costs(config, per_segment=True)


def load_costs_years(config):
    return _load_costs(config, per_year=True)


def _load_imp_exp(config, export=True, countries=None, carriers=None, years=None):
    """
    Return the imports or export of a subset of countries per country external to the subset
    for a given carrier. Since the network imports/exports are zero-sum, the exports can be obtained
    from the imports matrix

    """
    imp_exp = []
    df = pd.read_csv(Path(config["path"]["csvs"], "imports_exports.csv"), header=0)
    if countries is None:
        countries = df.columns[df.columns.str.match('[A-Z]{2}')]
    for y in years:
        imports_exports = 'exports' if export else 'imports'
        df_carrier = query_imp_exp(df.copy(), carriers, countries, y, imports_exports)
        imp_exp.append(df_carrier.rename(y))
    imp_exp = pd.concat(imp_exp, axis=1)
    imp_exp.index.rename('countries',inplace=True)
    return (
        imp_exp.loc[~(imp_exp == 0).all(axis=1)].reset_index()
    )


def load_marginal_prices(config):
    return pd.read_csv(Path(config["path"]["csvs"], "marginal_prices_countries.csv"))


def load_marginal_prices_t(config):
    return pd.read_csv(Path(config["path"]["csvs"], "marginal_prices_t_countries.csv"))


def load_buses(config):
    return pd.read_csv(Path(config["path"]["csvs"], "buses.csv"))


# %% Load main
def load_data_st(config):
    logger.info(f"Exporting data to streamlit")

    outputs = [
        "load_temporal",
        "load_temporal_be",
        "load_temporal_fl",
        "supply_temporal",
        "supply_temporal_be",
        "supply_temporal_fl",
        "supply_energy_df",
        "imports_exports",
        "res_temporal",
        "res_potentials",
        "power_capacities",
        "balancing_capacities",
        "balancing_supply",

        # Costs
        "costs_segments",
        "costs_years",
        "marginal_prices",
        "marginal_prices_t",

        # Geo
        "buses",
    ]

    dir = config["path"]["streamlit"]
    dir.mkdir(parents=True, exist_ok=True)
    for output in outputs:
        o = globals()["load_" + output](config)
        if isinstance(o, pd.DataFrame):
            o.to_csv(Path(dir, output + ".csv"), index=False)
        elif isinstance(o, dict):
            for k, v in o.items():
                # Determine sheet name
                sheet_name = output + "_" + k
                max_sheet_name_length = 31  # Limit sheet name to 31 char to enhance compatibility support (Excel)
                overflow_char = ".."
                sheet_name = (sheet_name[:max_sheet_name_length - len(overflow_char)] + overflow_char) \
                    if len(sheet_name) > max_sheet_name_length - len(overflow_char) else sheet_name
                v.to_csv(Path(dir, sheet_name + ".csv"), index=False)
        else:
            logging.warning(f"Given output for {output} is not mapped out to output file.")
