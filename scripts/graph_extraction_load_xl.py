# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2024 Climact for The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT

"""
Create data ready to present (load to excel)
"""
import logging
from pathlib import Path

import pandas as pd
from scripts.graph_extraction_utils import CLIP_VALUE_GW
from scripts.graph_extraction_utils import HEAT_RENAMER
from scripts.graph_extraction_utils import RES
from scripts.graph_extraction_utils import _load_supply_energy
from scripts.graph_extraction_utils import query_imp_exp
from scripts.graph_extraction_utils import LONG_TERM_STORAGE
from scripts.graph_extraction_utils import SHORT_TERM_STORAGE
from scripts.graph_extraction_utils import FF_ELEC
from scripts.graph_extraction_utils import FF_HEAT
from scripts.graph_extraction_utils import PRODUCTION
from scripts.graph_extraction_utils import H2
from scripts.graph_extraction_utils import BALANCE
from scripts.graph_extraction_utils import COST_SEGMENTS



logger = logging.getLogger(__name__)


# %% Unit countries capacities load
def _load_capacities(config, techs, historical="Historical (installed capacity by 2025)", countries=None):
    """

    Parameters
    ----------
    techs : list
        List of techs to filter on.
    historical : str, optional
        String to replace hist with. The default is "Historical (installed capacity by 2025)".
    countries : list, optional
        List of countries to filter on. The default is None.

    Returns
    -------
    DataFrame
        Generic load function, for which are specified the technologies to filter from. If a different
        name from historical is to be set, or if a subset of countries is needed, it can be provided.
.

    """
    techs
    df = (
        pd.read_csv(Path(config["path"]["csvs"], "units_capacities_countries.csv"), header=0)
        .drop(columns=["units"])
        .query("carrier in @techs")
    )
    countries
    if countries:
        df = df.query("node in @countries")

    idx = ["carrier"]
    df = (
        df.groupby(by=idx).sum().reset_index()
        .reindex(columns=config["excel_columns"]["all_years"])
        .rename(columns={"hist": historical})
    )

    df = df.set_index(idx)
    df = df[df.sum(axis=1) >= CLIP_VALUE_GW * (len(config["scenario"]["planning_horizons"]) + 1)]
    df = df.reset_index()

    return df


def load_capacities(config, tech_list, historical):
    """
    Generic function loading for each countries subset the type given in input and the name associated,
    as well as the historical memo.

    """
    dico = {}
    for co_name, subset in config["countries"].items():
        dico[f"{co_name}"] = _load_capacities(config, tech_list, countries=subset, historical=historical)
    return dico


def load_res_capacities(config):
    return load_capacities(config, RES, historical="Historical (planned by 2022)")


def load_production(config):
    return load_capacities(config, PRODUCTION, historical="Historical (installed capacity by 2025)")


def load_balance(config):
    return load_capacities(config, BALANCE, historical="Historical")


def load_long_term_storage(config):
    return load_capacities(config, LONG_TERM_STORAGE, historical="Historical (installed capacity by 2025)")


def load_short_term_storage(config):
    return load_capacities(config, SHORT_TERM_STORAGE, historical="Historical (installed capacity by 2025)")


def load_fossil_fuels(config):
    return load_capacities(config, FF_ELEC + FF_HEAT, historical="Historical (installed capacity by 2025)")


def load_h2_capacities(config):
    return load_capacities(config, H2, historical="Historical (installed capacity by 2025)")


# %% Costs load

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

    else:
        countries = prices.countries.unique()

    cost_mapping = pd.read_csv(
        Path(config["path"]["analysis_path"].resolve().parents[1], "data", "cost_mapping.csv"), index_col=[0, 1],
        header=0).dropna()
    df = (
        df.merge(cost_mapping, left_on=["carrier", "type"], right_index=True, how="left")
        .groupby(["cost_segment", "cost"]).sum(numeric_only=True)
        .reset_index()
    )

    if cost_segment:
        net_cost = pd.DataFrame([], columns=config["imp_exp_carriers"], index=config["scenario"]["planning_horizons"])
        if cost_segment != "Net_Imports":
            df = df.query('cost_segment in @cost_segment')

        if cost_segment == "Production" or cost_segment == "Net_Imports":
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
            df.loc[df.query("'fuel' in cost").index, config["years_str"]] += net_cost.sum(axis=1).values
        if cost_segment == "Net_Imports":
            df = net_cost.reset_index()

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
                dico[f"{seg_name}_{co_name}"] = _load_costs_year_segment(config, _countries=subset, cost_segment=seg)
                # print(seg_name,seg)
        elif per_year:
            for y in config["years_str"]:
                dico[f"{y}_{co_name}"] = _load_costs_year_segment(config, _countries=subset, year=y)
        else:
            logging.warning("Unkown configuration to load costs.")
    return dico


def load_costs_years(config):
    return _load_costs(config, per_year=True)


def load_costs_segments(config):
    return _load_costs(config, per_segment=True)


# %%Supply energy
def _load_supply_energy_dico(config, load=True, countries=None):
    """
    Allow to split _load_supply_energy into its carriers for it
    to be written tab by tab

    """
    dico = {}
    supply_energy_carrier = (
        _load_supply_energy(config, load=load, countries=countries)
        .carrier
        .replace({"AC": "electricity", "low voltage": "electricity"})
        .unique()
    )
    # TODO add h2 and gas imports and exports
    for ca in supply_energy_carrier:
        # TODO : should we move the HV/LV/imports/exports to the calling function to keep this function read only (no modifications) ?
        if ca == "electricity":
            df_ac = _load_supply_energy(config, load=load, countries=countries, carriers="AC")
            df_low = _load_supply_energy(config, load=load, countries=countries, carriers="low voltage")
            df = pd.concat([df_ac, df_low])
            del df["carrier"]
            df = df.groupby(by="sector").sum().reset_index()

            if not (load) and countries:
                df_imp = _load_imp_exp(config, export=False, countries=countries, carriers='elec',
                                       years=config["scenario"]["planning_horizons"])
                df_exp = _load_imp_exp(config, export=True, countries=countries, carriers='elec',
                                       years=config["scenario"]["planning_horizons"])
                df_net_imp = (df_imp[config["scenario"]["planning_horizons"]] - df_exp[
                    config["scenario"]["planning_horizons"]]).sum()
                df = pd.concat([df, pd.DataFrame(['Net Imports'] + df_net_imp.values.tolist(), index=df.columns).T])
            df = pd.concat([pd.DataFrame(df.sum(axis=0).values.tolist(), index=df.columns).T, df])
            df.iloc[0, 0] = 'Total'
            df.drop(df.query('sector == "V2G"').index, inplace=True)
            dico[ca] = df
        else:
            dico[ca] = _load_supply_energy(config, load=load, countries=countries, carriers=ca)

    for heat in ['dec_heat', 'cent_heat']:
        carriers = [k for k in list(dico.keys()) if HEAT_RENAMER.get(k) == heat]
        df = pd.DataFrame()
        for d in carriers:
            df = pd.concat([df, dico[d]])
            dico.pop(d)
        df = (df.groupby('sector')
              .agg({c: 'sum' if c in config['years_str'] else 'first' for c in df.columns})
              .drop('carrier', axis=1))
        dico[heat] = df
    return dico


def _load_imp_exp(config, export=True, countries=None, carriers=None, years=None):
    """
    Return the imports or export of a subset of countries per country external to the subset
    for a given carrier. Since the network imports/exports are zero-sum, the exports can be obtained
    from the imports matrix

    """
    imp_exp = []
    df = pd.read_csv(Path(config["path"]["csvs"], "imports_exports.csv"), header=0)
    for y in years:
        imports_exports = 'exports' if export else 'imports'
        df_carrier = query_imp_exp(df.copy(), carriers, countries, y, imports_exports)
        imp_exp.append(df_carrier.rename(y))
    imp_exp = pd.concat(imp_exp, axis=1)
    # TODO treat the display of exports/imports that does not specifically have the same countries displayed
    return (
        imp_exp.loc[~(imp_exp == 0).all(axis=1)].reset_index()
    )


def load_imports_exports(config):
    """
    This function loads the imports and exports per countries susbet with the
    rest of the system, for all countries subset and per carrier.
    If the countries subset is None, imports and exports are given per country
    with the rest of the system:
        e.g. :
            * countries = {'tot' : None} will call for
                - a DataFrame with the total imports per country and for each carrier
                - a DataFrame with the total exports per country and for each carrier
            * countries = {'EU27' : EU27_COUNTRIES} will call for
                - a DataFrame with the total imports of EU27 countries from external countries
                    per country and for each carrier
                - a DataFrame with the total exports of EU27 countries to external countries
                    per country and for each carrier

    """
    dico = {}

    for imp_exp, exp_value in {'imp': False, 'exp': True}.items():
        for ca in config["imp_exp_carriers"]:
            for name, subset in config["countries"].items():
                dico[f"{imp_exp}_{ca}_{name}"] = _load_imp_exp(config, export=exp_value,
                                                               countries=subset, carriers=ca,
                                                               years=config["scenario"]["planning_horizons"])
    return dico


def load_load_sectors(config):
    dico = {}
    for co_name, subset in config["countries"].items():
        to_rename = _load_supply_energy_dico(config, load=True, countries=subset)
        for k, v in to_rename.items():
            dico[f"{co_name}_{k}"] = v
    return dico


def load_supply_sectors(config):
    dico = {}
    for co_name, subset in config["countries"].items():
        to_rename = _load_supply_energy_dico(config, load=False, countries=subset)
        for k, v in to_rename.items():
            dico[f"{co_name}_{k}"] = v
    return dico


def load_supply_heat_be(config):
    dico = _load_supply_energy_dico(config, load=False, countries=["BE"])
    data = pd.DataFrame()
    df = {}

    for k, v in dico.items():
        to_add = v.copy()
        to_add['carrier'] = HEAT_RENAMER.get(k)
        data = pd.concat([data, to_add], ignore_index=True)

    for heat in ["dec_heat", "cent_heat"]:
        df[heat] = (
            data.loc[data.carrier.isin([heat])]
            .drop(columns="carrier")
            .groupby('sector').sum()
            .reset_index()
        )
    return df


# %% Non standard loads


def load_gas_phase_out(config):
    return (
        pd.read_csv(Path(config["path"]["csvs"], "gas_phase_out.csv"), header=0)
        .reindex(columns=config["excel_columns"]["first_hist_units"])
        .rename(columns={"hist": "Historical (planned by 2025)"})
    )


def load_grid_capacities(config):
    return (
        pd.read_csv(Path(config["path"]["csvs"], "grid_capacities.csv"), header=0)
        .reindex(columns=config["excel_columns"]["all_years_units"])
        .rename(columns={"hist": "Historical (planned by 2025)"})
    )


def load_grid_capacities_countries(config):
    return (
        pd.read_csv(Path(config["path"]["csvs"], "grid_capacities_countries.csv"), header=0)
        .groupby('Year').sum(numeric_only=True)
        .loc[:, ['LU', 'GB', 'NL', 'DE', 'FR']]
        .reset_index()
        .rename(columns={"hist": "Historical (planned by 2025)"})
    )


def load_h2_network_capacities(config):
    return (
        pd.read_csv(Path(config["path"]["csvs"], "H2_network_capacities.csv"), header=0)
        .reindex(columns=config["excel_columns"]["all_years_units"])
        .rename(columns={"hist": "Historical (planned by 2025)"})
    )


def load_h2_network_capacities_countries(config):
    return (
        pd.read_csv(Path(config["path"]["csvs"], "H2_network_capacities_countries.csv"), header=0)
        .groupby('Year').sum(numeric_only=True)
        .loc[:, ['LU', 'GB', 'NL', 'DE', 'FR']]
        # .rename(columns={"hist": "Historical (planned by 2025)"})
        .reset_index()
    )


def load_res_potentials(config):
    return (
        pd.read_csv(Path(config["path"]["csvs"], "res_potentials.csv"), header=0)
        .drop(columns=config["years_str"][:-1])
        .groupby(by="carrier").agg({config["years_str"][-1]: "sum", "units": "first"}).reset_index()
        .reindex(columns=config["excel_columns"]["last_units"])
    )


def load_res_potentials_be(config):
    return (
        pd.read_csv(Path(config["path"]["csvs"], "res_potentials.csv"), header=0)
        .query("region == 'BE'")
        .drop(columns=config["years_str"][:-1])
        .reindex(columns=config["excel_columns"]["last_units"])
    )


# def load_production_profile(config):
#     return (
#         pd.read_csv(Path(config["path"]["csvs"], "generation_profiles.csv"), header=0)
#         [["Year", "Country", "Carrier", "Annual sum [TWh]"]]
#         .query("Carrier in ['ammonia cracker', 'battery charger', 'H2 Electrolysis', 'Haber-Bosch', 'helmet', "
#                "'home battery charger', 'Sabatier']")
#         .groupby(by=["Year", "Carrier"]).sum(numeric_only=True).reset_index()
#         .rename(columns={"Carrier": "carrier"})
#         .replace(RENAMER)
#         .pivot_table(index="carrier", columns="Year", values="Annual sum [TWh]")
#         .reset_index()
#     )

def load_data_xl(config):
    logger.info(f"Exporting data to excel")

    outputs = [
        # Standard load
        "res_capacities",
        "production",
        "balance",
        "long_term_storage",
        "short_term_storage",
        "fossil_fuels",
        # "h2_capacities",
        "load_sectors",
        "supply_sectors",
        "supply_heat_be",

        # Costs
        "costs_years",
        "costs_segments",

        # Non standard
        "gas_phase_out",
        "grid_capacities",
        "grid_capacities_countries",
        "h2_network_capacities",
        "h2_network_capacities_countries",
        "res_potentials",
        "res_potentials_be",
        "imports_exports",
        # "production_profile",
    ]

    dir = config["path"]["excel"]
    dir.mkdir(parents=True, exist_ok=True)
    with pd.ExcelWriter(Path(dir, "graph_extraction_raw.xlsx")) as xl:
        for output in outputs:
            o = globals()["load_" + output](config)
            if isinstance(o, pd.DataFrame):
                o.to_excel(xl, sheet_name=output, index=False)
            elif isinstance(o, dict):
                for k, v in o.items():
                    # Determine sheet name
                    sheet_name = output + "_" + k
                    max_sheet_name_length = 31  # Limit sheet name to 31 char to enhance compatibility support (Excel)
                    overflow_char = ".."
                    sheet_name = (sheet_name[:max_sheet_name_length - len(overflow_char)] + overflow_char) \
                        if len(sheet_name) > max_sheet_name_length - len(overflow_char) else sheet_name
                    v.to_excel(xl, sheet_name=sheet_name, index=False)
            else:
                logging.warning(f"Given output for {output} is not mapped out to output file.")
    return
