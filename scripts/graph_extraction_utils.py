# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2024 Climact for The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT

"""
Create data ready to present (utils)
"""
import logging
import re
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

CLIP_VALUE_TWH = 1e-1  # TWh
CLIP_VALUE_ANNUAL_MWH = 1e3  # MWh
CLIP_VALUE_GW = 1e-3  # GW
RES = ["solar", "solar rooftop", "offwind", "offwind-ac", "offwind-dc", "onwind"]
HYDRO = ['ror', 'hydro']
HEAT_RENAMER = {"residential rural heat": "dec_heat", "services rural heat": "dec_heat",
                "residential urban decentral heat": "dec_heat", "services urban decentral heat": "dec_heat",
                "urban central heat": "cent_heat", "rural heat": 'dec_heat', "urban decentral heat": "dec_heat"}
NICE_RENAMER = {'dec_heat': "Decentralized heat", 'cent_heat': "Centralized heat",
                'elec': 'Electricity', 'rural heat': "Decentralized heat", 'electricity': 'Electricity',
                'urban decentral heat': "Decentralized heat", 'urban central heat': "Centralized heat",
                'methanol': 'Methanol', 'oil': 'Oil', 'co2': 'CO2 emissions',
                "coal": "Coal", 'gas': 'Methane', 'solid biomass': 'Solid biomass',
                "hydro": "Hydro Dams", "offwind": "Offshore wind", "onwind": "Onshore wind",
                "solar": "Solar PV", "solar thermal": 'Solar thermal', 'ror': "Run-of-the-river",
                "solid biomass CHP": "Solid biomass CHP", "solid biomass CHP CC": "Solid biomass CHP CC",
                "co2 stored": "CO2 captured", "H2": "Hydrogen", "NH3": "Ammonia"}
ELEC_RENAMER = {'AC': 'elec', 'DC': 'elec', 'low voltage': 'elec'}
TRANSMISSION_RENAMER = {"AC": "elec", "DC": "elec", "H2 pipeline": "H2",
                        "H2 pipeline retrofitted": "H2", "gas pipeline": "gas", "gas pipeline new": "gas"}
PREFIX_TO_REMOVE = [
    "residential ",
    "services ",
    "urban ",
    "rural ",
    "central ",
    "decentral ",
]

LONG_LIST_LINKS = ["coal/lignite", "oil", "CCGT", "OCGT", "H2 Electrolysis", "H2 Fuel Cell", "battery charger",
                   "home battery charger", "Haber-Bosch", "Sabatier", "ammonia cracker", "helmeth", "SMR", "SMR CC",
                   "hydro"]
LONG_LIST_GENS = ["solar", "solar rooftop", "onwind", "offwind", "offwind-ac", "offwind-dc", "ror", "nuclear",
                  "urban central solid biomass CHP", "home battery", "battery", "H2", "NH3"]
LONG_TERM_STORAGE = ["H2", "NH3"]
SHORT_TERM_STORAGE = ["battery", "home battery", "EV batteries"]
FF_ELEC = ["OCGT", "CCGT", "coal/lignite"]
FF_HEAT = ["residential / services oil boiler", "residential / services gas boiler"]
PRODUCTION = FF_ELEC + ["PHS", "hydro", "nuclear", "urban central biomass CHP", "solid biomass"] + RES
H2 = ["H2 Electrolysis", "H2 Fuel Cell"]
BALANCE = H2 + ["battery discharger", "home battery discharger", "V2G", "ammonia cracker"]
COST_SEGMENTS = {'prod': 'Production', 'bal': 'Balancing', 'tran': 'Transport', 'net_imp': "Net_Imports"}

logger = logging.getLogger(__name__)


def load_config(config_file, analysis_path, dir_export, scenario=''):
    logger.info("Loading configuration")

    with open(config_file, 'r') as f:
        config = yaml.safe_load(f)

    run = config["run"]["name"] if not scenario else scenario
    config["path"] = {
        "analysis_path": analysis_path,
        "resources_path": Path(analysis_path, "resources", run),
        "results_path": Path(analysis_path, "results", run),
    }

    years = config["scenario"]["planning_horizons"]
    config["years_str"] = list(map(str, years))
    simpl = config["scenario"]["simpl"][0]
    cluster = config["scenario"]["clusters"][0]
    opts = config["scenario"]["opts"][0]
    sector_opts = config["scenario"]["sector_opts"][0]
    ll = config["scenario"]["ll"][0]
    config["label"] = (cluster, ll, sector_opts, years)

    config["n_name"] = f"elec_s{simpl}_{cluster}_l{ll}_{opts}_{sector_opts}_"
    config["path"]["csvs"] = Path(analysis_path, dir_export, f"{scenario}_{config['n_name']}")
    config["path"]["streamlit"] = Path(config["path"]["analysis_path"], "graph_extraction_st", scenario)
    config["path"]["excel"] = Path(config["path"]["analysis_path"], "graph_extraction_xl", scenario)

    # Todo : type cast all references to year into str or int
    config["excel_columns"] = {"all_years": ["carrier", "hist"] + config["years_str"],
                               "all_years_units": ["carrier", "hist"] + config["years_str"] + ["units"],
                               "future_years": ["carrier"] + config["years_str"],
                               "future_years_sector": ["carrier", "sector"] + config["years_str"],
                               "first_year_units": ["carrier"] + [config["years_str"][0], "units"],
                               "last_hist_units": ["carrier", "hist"] + [config["years_str"][-1], "units"],
                               "last_units": ["carrier"] + [config["years_str"][-1], "units"],
                               "first_hist_units": ["carrier", "hist"] + [config["years_str"][0], "units"]}

    return config


def remove_prefixes(label):
    for ptr in PREFIX_TO_REMOVE:
        if label[: len(ptr)] == ptr:
            label = label[len(ptr):]
    return label


def query_imp_exp(df, carriers, countries, year, imports_exports):
    df_imp_exp = (
        df.query("carriers == @carriers")
        .query("year == @year")
        .query("imports_exports == @imports_exports")
        .drop(["carriers", "year", "imports_exports"], axis=1)
        .set_index('countries')
    )
    
    if countries is not None:
        country_list = df_imp_exp.columns.intersection(countries)
    else:
        country_list = df_imp_exp.columns.intersection(df.countries.unique())
                   
    if countries is not None:
        df_imp_exp = (df_imp_exp
        .loc[df_imp_exp.index.symmetric_difference(country_list), country_list]
        )
    df_imp_exp = df_imp_exp.sum(axis=1)
    return df_imp_exp


bus_map = {"BE1 0": None, "BE1 1": "BX", "BE1 2": "WL", "FL1 0": "FL",
           "BE1 0 H2": None, "BE1 1 H2": "BX", "BE1 2 H2": "WL", "FL1 0 H2": "FL",
           "BE1 0 gas": None, "BE1 1 gas": "BX", "BE1 2 gas": "WL", "FL1 0 gas": "FL",}
rx = re.compile("BE1 [0-2].*")


def bus_mapper(x, n, column=None):
    if x in n.buses.index:
        return bus_map[x] if (column == "country") and rx.fullmatch(x) else n.buses.loc[x, column]
    else:
        return np.nan


def renamer_to_country(x):
    rx_c = re.compile(r"([A-Z]{2})[0-9]\s[0-9]")
    rx_eu = re.compile(r"(EU)\s")
    if rx_c.match(x):
        return rx_c.match(x).group(1)
    elif rx_eu.match(x):
        return rx_eu.match(x).group(1)
    else:
        return np.nan


def _load_supply_energy(config, load=True, carriers=None, countries=None, aggregate=True, temporal=False):
    """
    Load nodal supply energy data and aggregate on carrier and sector, given some conditions.
    :param load: If True, keep only load data (negatives values)
    :param carriers: If specified, keep only a given carrier
    :param countries: If specified, keep a specific list of countries
    :param aggregate: If specified, determine if data are aggregated
    :return:
    """
    file = "supply_energy_sectors.csv"
    if temporal:
        file = 'temporal_' + file
        if countries:
            file = file.replace(".csv", f"_{countries}.csv")
    df = (
        pd.read_csv(Path(config["path"]["csvs"], file), header=0)
    )

    def get_load_supply(x):
        if load:
            return x.where(x <= 0, np.nan) * -1
        else:
            return x.where(x > 0, np.nan)

    df[config["years_str"]] = df[config["years_str"]].apply(get_load_supply)
    df = df.dropna(subset=config["years_str"], how="all")

    if carriers:
        df = df.query("carrier in @carriers")
    if countries and not temporal:
        df = df.query("node in @countries")
    if aggregate:
        idx = ["carrier", "sector"]
        if temporal:
            idx.append('snapshot')
            df = (
                df.groupby(by=idx).sum().reset_index()
                .reindex(columns=config["excel_columns"]["future_years_sector"] + ['snapshot'])
            )
        else:
            df = (
                df.query("node != 'FL'")  # Flanders already removed from temporal
                .groupby(by=idx).sum().reset_index()
                .reindex(columns=config["excel_columns"]["future_years_sector"])
            )

    else:
        idx = ["carrier", "sector", "node"]
        if temporal:
            idx.append('snapshot')
            idx.remove('node')
        df = (
            df.groupby(by=idx).sum().reset_index()
            .reindex(columns=config["excel_columns"]["future_years_sector"] + ['node'])
        )

    if temporal:
        idx.remove('snapshot')
        index = (df
                 .groupby(idx)
                 .sum(numeric_only=True)
                 .sum(axis=1) >= CLIP_VALUE_TWH
                 )
        df = df.set_index(idx).loc[index].reset_index()
    else:
        df = df.set_index(idx)
        df = df[df.sum(axis=1) >= CLIP_VALUE_TWH]
        df = df.reset_index()

    return df


def groupby_bus(n, c, nice_names=True):
    if c in n.one_port_components:
        return [n.df(c).bus, n.df(c).carrier]
    else:
        return [n.df(c).bus1, n.df(c).carrier]


def groupby_buses(n, c, nice_names=True):
    if c in n.one_port_components:
        return [n.df(c).bus, n.df(c).carrier]
    else:
        return [n.df(c).bus0, n.df(c).bus1, n.df(c).carrier]
