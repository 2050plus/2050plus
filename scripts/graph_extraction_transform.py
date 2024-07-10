# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2024 Climact for The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT

"""
Create data ready to present (transform)
"""
import logging
import re
from pathlib import Path

import matplotlib as mpl
import pandas as pd
from matplotlib import pyplot as plt
from scripts.graph_extraction_utils import CLIP_VALUE_TWH
from scripts.graph_extraction_utils import ELEC_RENAMER
from scripts.graph_extraction_utils import HEAT_RENAMER
from scripts.graph_extraction_utils import HYDRO
from scripts.graph_extraction_utils import RES
from scripts.graph_extraction_utils import TRANSMISSION_RENAMER
from scripts.graph_extraction_utils import bus_mapper
from scripts.graph_extraction_utils import groupby_bus
from scripts.graph_extraction_utils import remove_prefixes
from scripts.graph_extraction_utils import renamer_to_country
from scripts.make_summary import calculate_nodal_supply_energy

from prepare_sector_network import prepare_costs

logger = logging.getLogger(__name__)

RES
RENAMER = {
    # Carriers
    "offwind-dc": "offwind",
    "offwind-ac": "offwind",
    "solar rooftop": "solar",
    "coal": "coal/lignite",
    "lignite": "coal/lignite",
    "battery storage": "EV batteries",

    # Boilers
    "residential rural biomass boiler": "residential / services biomass boiler",
    "residential urban decentral biomass boiler": "residential / services biomass boiler",
    "services rural biomass boiler": "residential / services biomass boiler",
    "services urban decentral biomass boiler": "residential / services biomass boiler",
    "residential rural gas boiler": "residential / services gas boiler",
    "residential urban decentral gas boiler": "residential / services gas boiler",
    "services rural gas boiler": "residential / services gas boiler",
    "services urban decentral gas boiler": "residential / services gas boiler",
    "urban central gas boiler": "residential / services gas boiler",
    "residential rural oil boiler": "residential / services oil boiler",
    "residential urban decentral oil boiler": "residential / services oil boiler",
    "services rural oil boiler": "residential / services oil boiler",
    "services urban decentral oil boiler": "residential / services oil boiler",
    "urban central oil boiler": "residential / services oil boiler",

    # Water tanks
    "residential rural water tanks charger": "residential / services water tanks charger",
    "residential urban decentral water tanks charger": "residential / services water tanks charger",
    "services rural water tanks charger": "residential / services water tanks charger",
    "services urban decentral water tanks charger": "residential / services water tanks charger",
    "urban central water tanks charger": "residential / services water tanks charger",
    "residential rural water tanks discharger": "residential / services water tanks discharger",
    "residential urban decentral water tanks discharger": "residential / services water tanks discharger",
    "services rural water tanks discharger": "residential / services water tanks discharger",
    "services urban decentral water tanks discharger": "residential / services water tanks discharger",
    "urban central water tanks discharger": "residential / services water tanks discharger",

    # Heat pumps
    "residential rural ground heat pump": "residential / services rural ground heat pump",
    "services rural ground heat pump": "residential / services rural ground heat pump",
    "residential urban decentral air heat pump": "residential / services air heat pump",
    "services urban decentral air heat pump": "residential / services air heat pump",
    "urban central air heat pump": "residential / services air heat pump",

    # Resistive heaters
    "residential rural resistive heater": "residential / services resistive heater",
    "residential urban decentral resistive heater": "residential / services resistive heater",
    "services rural resistive heater": "residential / services resistive heater",
    "services urban decentral resistive heater": "residential / services resistive heater",
    "urban central resistive heater": "residential / services resistive heater",

    # Solar thermals
    "residential rural solar thermal collector": "solar thermal",
    "residential urban decentral solar thermal collector": "solar thermal",
    "services rural solar thermal collector": "solar thermal",
    "services urban decentral solar thermal collector": "solar thermal",
    "urban central solar thermal collector": "solar thermal",

    # Solid biomass CHP
    "urban central solid biomass CHP": "solid biomass CHP",
    "urban central solid biomass CHP CC": "solid biomass CHP CC",

}


# %%  Extract loads

def extract_res_potential(n):
    """
    Extract renewable potentials in GW.
    :param n: Network
    :return: Potentials of renewables in GW
    """
    dfx = []
    rx = re.compile("([A-z]+)[0-9]+\s[0-9]+\s([A-z\-\s]+)-*([0-9]*)")

    for y, ni in n.items():
        df = pd.concat([ni.df("Generator"), ni.df("Link"), ni.df("StorageUnit")])[
            ["p_nom_max", "p_nom_opt"]].reset_index()
        df[["region", "carrier", "build_year"]] = df["index"].str.extract(rx)
        df["carrier"] = df["carrier"].str.rstrip("-").replace(RENAMER)
        df["planning horizon"] = y
        df = df.query("carrier in ['onwind', 'offwind', 'solar', 'solar thermal', 'ror', 'hydro', 'PHS']")
        dfx.append(
            df.groupby(["planning horizon", "carrier", "build_year", "region"]).sum(numeric_only=True) / 1e3
        )  # GW      

    dfx = pd.concat(dfx)
    # TODO : dfx["p_nom_opt"].index.get_level_values("build_year") ==  dfx.index.get_level_values("build_year")
    df_potential = pd.concat([
        (
            dfx.loc[
                dfx["p_nom_opt"].index.get_level_values("build_year") <
                dfx["p_nom_opt"].index.get_level_values("planning horizon").astype(
                    str)  # selects all rows whose build_year index has a smaller value than the planning_horizon index
                , "p_nom_opt"]  # if the year of construction is before planning horizon, this means that the generator has already been installed. However, what was installed = p_nom_opt
            .groupby(["planning horizon", "carrier", "region"]).sum()
            # all the lines that have the sames indexes ("planning horizon", "carrier", "region") but a different "build_year" will be summed based on p_nom_opt
        ),
        (
            dfx.loc[
                dfx["p_nom_max"].index.get_level_values("build_year") ==
                dfx["p_nom_max"].index.get_level_values("planning horizon").astype(str)
                , "p_nom_max"]
            .groupby(["planning horizon", "carrier", "region"]).sum()
        )
    ], axis=1)

    df_potential["potential"] = df_potential["p_nom_max"].fillna(0) + df_potential["p_nom_opt"].fillna(
        0)  # potential = what could have been installed (p_nom_opt) + what has been installed (p_nom_opt)
    df_potential = (
        df_potential.reset_index()
        .pivot(index=["carrier", "region"], columns="planning horizon", values="potential")
    )
    df_potential["units"] = "GW_e"
    return df_potential


# %%


def extract_inst_capa_elec_node(config, n, carriers_renamer):
    inst_capa_elec_node = []

    sector_mapping = pd.read_csv(
        Path(config["path"]["analysis_path"].resolve().parents[1], "data", "sector_mapping.csv"), index_col=[0, 1, 2],
        header=0).dropna()
    sector_mapping = sector_mapping.reset_index()
    sector_mapping["carrier"] = sector_mapping["carrier"].map(lambda x: carriers_renamer.get(x, x))
    sector_mapping.loc[sector_mapping.query("component == 'links'").index, "item"] = sector_mapping.loc[
        sector_mapping.query("component == 'links'").index, "item"].apply(lambda x: x[:-1])
    sector_mapping.component = sector_mapping.component.apply(
        {"generators": "Generator", "storage_units": "StorageUnit", "links": "Link"}.get)
    sector_mapping = sector_mapping.query("carrier=='elec'").groupby(["carrier", "component", "item"]).first()

    for y, ni in n.items():
        stats = ni.statistics
        inst_capa_elec_node_i = stats.optimal_capacity(bus_carrier=["AC", "low voltage"], groupby=groupby_bus)

        inst_capa_elec_node_i = (inst_capa_elec_node_i
                                 .rename(y)
                                 .reset_index()
                                 .rename(columns={"carrier": "item", "bus": "country", "bus1": "country"})
                                 .merge(sector_mapping.droplevel("carrier"), left_on=["component", "item"],
                                        right_index=True, how="left")
                                 .apply({"country": lambda x: x[:2], "sector": lambda x: x, y: lambda x: x})
                                 .groupby(["country", "sector"])
                                 .sum(numeric_only=True)
                                 / 1e3  # GW
                                 )

        inst_capa_elec_node.append(inst_capa_elec_node_i)

    return pd.concat(inst_capa_elec_node, axis=1)


def extract_balancing_data(method, n):
    balancing_data = []
    techno_to_keep = {"PHS", "hydro", "H2 Fuel Cell", "battery charger",
                      "home battery charger", "V2G", "BEV charger",
                      "ground heat pump", "air heat pump", "water tanks charger"}

    for y, ni in n.items():
        stats = ni.statistics

        data_i = (
                getattr(stats, method)(groupby=groupby_bus)
                .reset_index()
                .drop('component', axis=1)
                .rename(columns={"bus": "country", "bus1": "country", 0: y})
                .pipe(
                    lambda df: df.assign(carrier=df.carrier.apply(remove_prefixes).apply(lambda x: RENAMER.get(x, x))))
                .query("carrier.isin(@techno_to_keep)")
                .assign(country=lambda df: df["country"].map(ni.buses.country))
                .groupby(["country", "carrier"]).sum(numeric_only=True)
                / 1e3  # GW
        )

        balancing_data.append(data_i)

    balancing_data = pd.concat(balancing_data, axis=1)

    return balancing_data


def extract_electricity_network(n):
    # DC
    links_DC = []

    for y, ni in n.items():
        links_DC_i = ni.links[ni.links['carrier'] == 'DC']
        links_DC_i = links_DC_i[['bus0', 'bus1', 'length', 'p_nom', 'p_nom_max', 'carrier', 'p_nom_opt']]
        links_DC_i = links_DC_i.assign(year=str(y))
        links_DC_i = links_DC_i.reset_index()
        links_DC_i.rename(columns={'Link': 'Cable'}, inplace=True)
        links_DC_i.set_index(['Cable', 'year'], inplace=True)

        links_DC.append(links_DC_i)

    links_DC = pd.concat(links_DC)

    # AC
    lines_AC = []

    for y, ni in n.items():
        lines_AC_i = ni.lines[['bus0', 'bus1', 'length', 's_nom', 's_nom_max', 'carrier', 's_nom_opt']].copy()
        lines_AC_i.rename(columns={'s_nom': 'p_nom', 's_nom_max': 'p_nom_max', 's_nom_opt': 'p_nom_opt'}, inplace=True)
        lines_AC_i = lines_AC_i.assign(year=str(y))
        lines_AC_i = lines_AC_i.reset_index()
        lines_AC_i.rename(columns={'Line': 'Cable'}, inplace=True)
        lines_AC_i.set_index(['Cable', 'year'], inplace=True)

        lines_AC.append(lines_AC_i)

    lines_AC = pd.concat(lines_AC)

    elec_grid = pd.concat([links_DC, lines_AC])

    elec_grid[['p_nom', 'p_nom_max', 'p_nom_opt']] = elec_grid[['p_nom', 'p_nom_max', 'p_nom_opt']].div(1e3)

    return elec_grid


# %%

def extract_res_temporal_energy(config, n):
    HYDRO
    df = []
    for y, ni in n.items():
        units = pd.concat([ni.generators, ni.storage_units, ni.links])
        units_t = pd.concat([ni.generators_t.p, ni.storage_units_t.p, ni.links_t.p1], axis=1)
        res = units.copy().query("carrier in @RES or carrier in @HYDRO or " \
                                 "carrier.str.contains('solar thermal') or carrier.str.contains('solid biomass CHP')")
        res_t = units_t[res.index]
        res.loc[res.bus.isna(), "bus"] = res.loc[res.bus.isna(), "bus1"]
        res["country"] = res.bus.map(renamer_to_country)
        res["carrier"] = res.carrier.apply(remove_prefixes).apply(lambda x: RENAMER.get(x, x))

        res_t.columns = res_t.columns.map(lambda x: (y, res.loc[x].carrier, res.loc[x].country))
        res_t.rename_axis(["year", "carrier", "country"], axis=1, inplace=True)
        res_t = res_t.T.groupby(["year", "carrier", "country"]).sum().T.abs().astype(float)
        df.append(res_t)
    df = pd.concat(df, axis=1) / 1e3  # GW
    return df.T


def calculate_imp_exp(country_map, transmission_t, y):
    countries = sorted(country_map.stack().unique())

    table_li_co = pd.DataFrame([], index=country_map.index)
    other_bus = pd.DataFrame([], index=transmission_t.columns)
    mat_imp = pd.DataFrame([], columns=countries, index=countries).rename_axis(index="countries")
    mat_exp = pd.DataFrame([], columns=countries, index=countries).rename_axis(index="countries")
    mat_imp["year"] = y
    mat_exp["year"] = y

    for co in countries:
        if "hist" != y:
            table_li_co[co] = country_map.apply(lambda x: -1 if x.bus0 == co else 0, axis=1)
            table_li_co[co] += country_map.apply(lambda x: 1 if x.bus1 == co else 0, axis=1)

            other_bus[co] = country_map.apply(lambda x: x.bus1 if x.bus0 == co else "", axis=1)
            other_bus[co] += country_map.apply(lambda x: x.bus0 if x.bus1 == co else "", axis=1)

            ie_raw = transmission_t.mul(table_li_co[co]) / 1e6  # TWh
            imp = ie_raw.where(ie_raw > 0, 0).sum(axis=0)
            exp = ie_raw.mask(ie_raw > 0, 0).sum(axis=0)

            idx1 = other_bus[co].loc[imp[imp > CLIP_VALUE_TWH].index].replace('', None).dropna()
            mat_imp.loc[idx1, co] = imp[idx1.index].values
            idx2 = other_bus[co].loc[exp[exp < -CLIP_VALUE_TWH].index].replace('', None).dropna()
            mat_exp.loc[idx2, co] = -exp[idx2.index].values

    return mat_imp.fillna(0), mat_exp.fillna(0), table_li_co, other_bus


def extract_transmission(n, carriers=["AC", "DC"],
                         units={"AC": "GW_e", "DC": "GW_e",
                                "gas pipeline": "GW_lhv,ch4", "gas pipeline new": "GW_lhv,ch4",
                                "H2 pipeline": "GW_lhv,h2", "H2 pipeline retrofitted": "GW_lhv,h2"}):
    capacities = []
    capacities_countries = []
    imports, exports = [], []
    # Add projected values
    for y, ni in n.items():

        transmission = []
        if "hist" != y:
            transmission_t = []
        for ca in carriers:
            if ca == "AC":
                transmission.append(ni.lines.rename(columns={"s_nom_opt": "p_nom_opt"}))
                if "hist" != y:
                    transmission_t.append(ni.lines_t.p0 * 8760 / len(ni.snapshots))
            else:
                transmission.append(ni.links.query("carrier == @ca"))
                if "hist" != y:
                    transmission_t.append(
                        ni.links_t.p0[ni.links.query("carrier == @ca").index] * 8760 / len(ni.snapshots))

        transmission = pd.concat(transmission)
        if "hist" != y:
            transmission_t = pd.concat(transmission_t, axis=1)

        buses_links = [c for c in transmission.columns if "bus" in c]
        country_map = transmission[buses_links].map(lambda x: bus_mapper(x, ni, column="country"))
        transmission_co = {}
        mono_co = {}
        for co in ni.buses.country.unique():
            transmission_co[co] = (transmission
                                   .query("@co == @country_map.bus0 or @co == @country_map.bus1")
                                   .groupby("carrier")
                                   .p_nom_opt.sum()
                                   )

            mono_co[co] = (
                transmission.loc[(transmission.index.str.contains("->")) & (transmission.index.str.contains("<"))]
                .query("@co == @country_map.bus0")
                .groupby("carrier")
                .p_nom_opt.sum()
            )

            if len(mono_co[co]):
                transmission_co[co].loc[mono_co[co].index] -= mono_co[co]

        transmission_co = pd.DataFrame.from_dict(transmission_co, orient="columns").fillna(0) / 1e3
        capacities_countries.append(pd.concat({y: transmission_co}, names=["Year"]))

        transmission_total = pd.DataFrame(transmission.groupby("carrier").p_nom_opt.sum()) / 1e3
        capacities.append(transmission_total.rename(columns={"p_nom_opt": y}))

        # imports/exports
        mat_imp, mat_exp, table_li_co, _ = calculate_imp_exp(country_map, transmission_t, y)
        if "hist" != y:
            imports.append(mat_imp.reset_index().set_index(["countries", "year"]))
            exports.append(mat_exp.reset_index().set_index(["countries", "year"]))

    df = pd.concat(capacities, axis=1)
    df_co = pd.concat(capacities_countries, axis=0)
    df_imp = pd.concat(imports, axis=0).fillna(0)
    df_imp["imports_exports"] = "imports"
    df_exp = pd.concat(exports, axis=0).fillna(0)
    df_exp["imports_exports"] = "exports"

    df["units"] = df.index.map(units)
    df_co["units"] = df_co.index.get_level_values(level=1).map(units)
    df_imp_exp = pd.concat([df_imp, df_exp])
    df_imp_exp["carriers"] = TRANSMISSION_RENAMER.get(carriers[0])

    # FixMe Remove double counting in df (FLanders)
    logger.warning("Function output double count Flanders transmission capacities in df")
    return df, df_co, df_imp_exp


def extract_nodal_costs(config, n):
    # Todo : add handling of multiple runs
    def costs_groupby(n, c, nice_names=True):
        if c in n.one_port_components:
            return [n.df(c).location, n.df(c).carrier]
        else:
            return [n.df(c).location, n.df(c).carrier, n.df(c).build_year.astype(int), n.df(c).lifetime]

    # Prepare costs for
    costs = prepare_costs(Path(config["path"]["resources_path"], "costs_" + config["years_str"][0] + ".csv"),
                          config["costs"], 1 / 3)

    gas_distribution_old = costs.query('index.str.contains("gas distribution")').eval("investment*FOM/100*1/3").mean()
    gas_distribution = costs.query('index.str.contains("electricity distribution")').fixed.mean() * config["sector"][
        "gas_distribution_grid_cost_factor"]

    df_comp = {}
    for y, ni in n.items():
        df_agg = {}
        for c in ni.iterate_components(["Generator", "Link", "Line", "StorageUnit", "Store"]):

            capex = ni.statistics.capex(comps=c.name, groupby=costs_groupby).rename("capital")
            opex = ni.statistics.opex(comps=c.name, groupby=costs_groupby).rename("marginal")
            capacities = ni.statistics.optimal_capacity(comps=c.name, groupby=costs_groupby).rename("capacities")
            df = pd.concat([capex, opex, capacities], axis=1).fillna(0)

            if c.list_name == "links":
                df["distribution_cost"] = 0.0
                ind_gas_distri = (df.query('(carrier.str.contains("gas boiler") or carrier.str.contains("micro CHP"))' \
                                           ' and not(carrier.str.contains("urban central"))')).index
                ind_old = df.query("build_year<2030").index
                df.loc[ind_gas_distri.intersection(ind_old), "distribution_cost"] = gas_distribution_old
                df.loc[ind_gas_distri.difference(ind_old), "distribution_cost"] = gas_distribution
                df["capital"] -= df["distribution_cost"] * df.capacities
                df_gas_distri = df.loc[ind_gas_distri]
                df_gas_distri.index = df_gas_distri.index.set_levels(
                    df_gas_distri.index.levels[1] + " gas distribution", level=1)
                df_gas_distri["capital"] = df_gas_distri["distribution_cost"] * df_gas_distri.capacities
                df_gas_distri["marginal"] = 0
                df = pd.concat([df, df_gas_distri])
                df.index = df.index.droplevel(3)
            df = pd.melt(df, value_vars=["capital", "marginal"], var_name="cost", value_name=y, ignore_index=False)
            df = df.groupby(["cost", "location", "carrier"]).sum()
            df_agg[c.list_name] = df
        df_comp[y] = pd.concat(df_agg.values(), keys=df_agg.keys())
        df_comp[y].index = df_comp[y].index.set_names("type", level=0)
    df_comp = pd.concat(df_comp.values(), axis=1)

    def renamer_to_country_costs(x):
        rx_eu = re.compile(r"(EU)\s")
        rx_c = re.compile(r"([A-Z]{2})")
        if rx_eu.match(x):
            return rx_eu.match(x).group(1)
        elif rx_c.match(x):
            return rx_c.match(x).group(1)
        else:
            return "EU"

    df_comp.reset_index(inplace=True)
    df_comp.rename(columns={"location": "country"}, inplace=True)
    df_comp.country = df_comp.country.map(renamer_to_country_costs)
    fuels = df_comp.query(
        "carrier in ['gas','oil','coal','lignite','uranium','H2'] and cost == 'marginal' and type == 'generators'").index
    biomass = df_comp.query("(carrier.str.contains('biomass') or carrier.str.contains('biogas')) and"
                            " cost == 'marginal' and type == 'stores'").index
    # df_comp.loc[fuels.union(biomass), "cost"] = "fuel"  # deactivated for now
    df_comp = df_comp.set_index(["type", "cost", "country", "carrier"])
    df_comp = df_comp.fillna(0).groupby(["type", "cost", "country", "carrier"]).sum()
    df_comp = df_comp.loc[~df_comp.apply(lambda x: x < 1e3).all(axis=1)]
    df_comp.insert(0, column="units", value="Euro")

    cost_mapping = pd.read_csv(
        Path(config["path"]["analysis_path"].resolve().parents[1], "data", "cost_mapping.csv"), index_col=[0, 1],
        header=0).dropna()
    df_comp = (
        df_comp.merge(cost_mapping, left_on=["carrier", "type"], right_index=True, how="left")
    )
    return df_comp


def extract_marginal_prices(n, carrier_list=["AC"]):
    df = []
    df_t = []
    for ca in carrier_list:
        prices = pd.DataFrame([]).rename_axis(index="countries")
        prices_t = pd.DataFrame([])

        for y, ni in n.items():
            if "hist" != y:
                price_y = (
                    ni.buses_t.marginal_price[ni.buses.query("carrier == @ca ").index]
                )
                prices[y] = price_y.mean().groupby(lambda x: bus_mapper(x, ni, column="country")).mean()
                prices[f"{y}_std"] = price_y.std().groupby(lambda x: bus_mapper(x, ni, column="country")).mean()
                marginal_t = price_y.T.groupby(lambda x: bus_mapper(x, ni, column="country")).mean().rename_axis(
                    index=["countries"])
                marginal_t["year"] = y

                prices_t = pd.concat([prices_t, marginal_t.reset_index().set_index(["year", "countries"])])

        prices["carrier"] = ca.replace("AC", "elec")
        prices_t["carrier"] = ca.replace("AC", "elec")
        df.append(prices.reset_index().set_index(["countries", "carrier"]))
        df_t.append(prices_t.reset_index().set_index(["countries", "carrier", "year"]))
    df = pd.concat(df, axis=0)
    df_t = pd.concat(df_t)
    return df, df_t


def extract_nodal_supply_energy(config, n):
    labels = {y: config["label"][:-1] + (y,) for y in n.keys()}
    columns = pd.MultiIndex.from_tuples(labels.values(), names=["cluster", "ll", "opt", "planning_horizon"])
    df = pd.DataFrame(columns=columns, dtype=float)
    for y, ni in n.items():
        df = calculate_nodal_supply_energy(ni, label=labels[y], nodal_supply_energy=df, country_aggregate=False)
    idx = ["node", "carrier", "component", "item"]
    df.index.names = idx
    df.columns = df.columns.get_level_values(3)

    df = df.reset_index()
    df["node"] = df["node"].map(renamer_to_country)
    df = df.dropna(how="all").set_index(idx)

    df = df * 1e-6  # TWh
    df["units"] = "TWh"

    sector_mapping = pd.read_csv(
        Path(config["path"]["analysis_path"].resolve().parents[1], "data", "sector_mapping.csv"), index_col=[0, 1, 2],
        header=0).dropna()
    df = df.merge(sector_mapping, left_on=["carrier", "component", "item"], right_index=True, how="left")
    return df


def extract_temporal_supply_energy(config, n, carriers=None, carriers_renamer=None,
                                   time_aggregate=False, country_aggregate=True):
    labels = {y: config["label"][:-1] + (y,) for y in n.keys()}
    columns = pd.MultiIndex.from_tuples(labels.values(), names=["cluster", "ll", "opt", "planning_horizon"])
    df = pd.DataFrame(columns=columns, dtype=float)

    sector_mapping = pd.read_csv(
        Path(config["path"]["analysis_path"].resolve().parents[1], "data", "sector_mapping.csv"), index_col=[0, 1, 2],
        header=0).dropna()
    sector_mapping = sector_mapping.reset_index()
    sector_mapping["carrier"] = sector_mapping["carrier"].map(lambda x: carriers_renamer.get(x, x))
    sector_mapping["item"] = sector_mapping["item"].map(remove_prefixes)
    if carriers is None:
        carriers = sector_mapping["carrier"].unique()
    sector_mapping = sector_mapping.groupby(["carrier", "component", "item"]).first()

    # Snapshot values are weigthed by the timestep length (hence if timestep length = 3h, snapshot values are equal to 3h consumption)
    # Hence snapshot sum is equal to the annual sum
    for y, ni in n.items():
        df = calculate_nodal_supply_energy(ni, label=labels[y], nodal_supply_energy=df,
                                           carriers=carriers, carriers_renamer=carriers_renamer,
                                           time_aggregate=time_aggregate,
                                           country_aggregate=country_aggregate, fun=remove_prefixes)
        logger.info(f"Calculated nodal supply for year {y}")
    idx = ["carrier", "component", "item", "snapshot"]
    if not (country_aggregate):
        idx = ["node"] + idx
    df.index.names = idx
    df.columns = df.columns.get_level_values(3)

    df = df * 1e-3  # GWh

    if not (country_aggregate):
        df = df.reset_index()
        df["node"] = df["node"].map(renamer_to_country)
        # Filter on annual total
        df = df.dropna().loc[(df.groupby(["node", "carrier", "component", "item"])
                              [config["scenario"]["planning_horizons"]]
                              .filter(lambda x: (x.sum(axis=0) > 1e2).any())
                              ).index].set_index(idx)

    # Convert hourly values back to 1h-consumption for visualisation
    df = df * len(ni.snapshots) / 8760
    df = df.merge(sector_mapping, left_on=["carrier", "component", "item"], right_index=True, how="left").dropna(axis=0,
                                                                                                                 subset="sector")
    return df


def get_state_of_charge_t(n, carrier):
    df = n.storage_units_t.state_of_charge.T.reset_index()
    df = df.merge(n.storage_units.reset_index()[["carrier", "StorageUnit"]], on="StorageUnit")
    df = df.groupby(by="carrier").sum()
    df.drop(columns=["StorageUnit"], inplace=True)
    return df.T[[carrier]] / 1e6  # TWh


def get_e_t(n, carrier):
    df = n.stores_t.e.T.reset_index()
    df = df.merge(n.stores.reset_index()[["carrier", "Store"]], on="Store")
    df = df.groupby(by="carrier").sum()
    df.drop(columns=["Store"], inplace=True)
    return df.T[[carrier]] / 1e6  # TWh


def get_p_carrier_nom_t(n, carrier):
    df = n.links_t.p_carrier_nom_opt.T.reset_index()
    df = df.merge(n.links.reset_index()[["carrier", "Link"]], on="Link")
    df = df.groupby(by="carrier").sum()
    df.drop(columns=["Link"], inplace=True)
    return df.T[[carrier]] / 1e3  # GW


def _extract_graphs(config, n, storage_function, storage_horizon, both=False, units={}, color_shift=None):
    carrier = list(storage_horizon.keys())
    if color_shift:
        pass
    else:
        color_shift = dict(zip(config["scenario"]["planning_horizons"],
                               ["C" + str(i) for i in range(len(config["scenario"]["planning_horizons"]))]))
    fig = plt.figure(figsize=(14, 8))

    def plotting(ax, title, data, y, unit):
        data.index = pd.to_datetime(pd.DatetimeIndex(data.index.values).strftime("2030-%m-%d-%H"))
        ax.plot(data, label=y, color=color_shift.get(y))
        ax.set_title(title)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.set_ylabel(unit, ha="left", y=1.1, rotation=0, labelpad=0.2)
        plt.xticks(rotation=30)
        plt.tight_layout()
        return

    for y, ni in n.items():
        lt = {}
        st = {}
        for car in carrier:
            storage = globals()[storage_function.get(car, "get_e_t")](ni, car)
            if both:
                lt[car] = storage
                st[car] = storage.iloc[:int(8 * 31)]
            elif "L" in storage_horizon.get(car):
                lt[car] = storage
            else:
                st[car] = storage.iloc[:int(8 * 31)]
        for i, (car, s) in enumerate(st.items()):
            ax = plt.subplot(3, 2, 2 * i + 1)
            plotting(ax, car, s, y, unit=r"${}$".format(units.get(car, "[TWh]")))
        for i, (car, l) in enumerate((lt).items()):
            ax = plt.subplot(3, 2, 2 * (i + 1))
            plotting(ax, car, l, y, unit=r"${}$".format(units.get(car, "[TWh]")))

    ax.legend()
    return {"": fig}


def extract_graphs(config, n, color_shift=None):
    ## Figures to extract
    plt.close("all")
    # Storage
    mpl.rcParams.update(mpl.rcParamsDefault)
    if color_shift:
        pass
    else:
        color_shift = dict(zip(config["scenario"]["planning_horizons"], ["C0", "C2", "C1"]))

    storage_function = {"hydro": "get_state_of_charge_t", "PHS": "get_state_of_charge_t"}
    storage_horizon = {"hydro": "LT", "PHS": "ST", "H2": "LT",
                       "battery": "ST", "home battery": "ST",
                       "NH3": "LT"}
    n_sto = _extract_graphs(config, n, storage_function, storage_horizon, color_shift=color_shift)
    # h2
    storage_function = {"H2 Fuel Cell": "get_p_carrier_nom_t", "H2 Electrolysis": "get_p_carrier_nom_t"}
    storage_horizon = {"H2 Fuel Cell": "LT", "H2 Electrolysis": "LT", "H2": "LT"}
    n_h2 = _extract_graphs(config, n, storage_function, storage_horizon, color_shift=color_shift,
                           both=True, units={"H2 Fuel Cell": "[GW_e]", "H2 Electrolysis": "[GW_e]",
                                             "H2": "[TWh_{lhv,h2}]"})
    return n_sto, n_h2


def export_csvs_figures(csvs, outputs, figures):
    csvs.mkdir(parents=True, exist_ok=True)

    for f_name, f in figures.items():
        for y, plot in f.items():
            if "px" in f_name:
                plot.write_html(Path(csvs, f"{f_name}_{y}.html"))
            else:
                plot.savefig(Path(csvs, f"{f_name}_{y}.png"), transparent=True)

    for o_name, o in outputs.items():
        o.to_csv(Path(csvs, f"{o_name}.csv"))

    logger.info(f"Exported files and figures to folder : {csvs}")

    return


def extract_geo_buses(n):
    buses = (
        next(iter(n.values())).buses
        .query("carrier=='AC'")
        [["x", "y", "country"]]
        .rename(columns={"x": "lon", "y": "lat"})
        .drop(["BE1 0", "BE1 2"])  # Keep only one coordinate set for Belgium
    )
    buses = buses.set_index("country")
    return buses


def transform_data(config, n, n_ext, color_shift=None):
    logger.info("Transforming data")

    carriers_renamer = {}
    carriers_renamer.update(HEAT_RENAMER)
    carriers_renamer.update(ELEC_RENAMER)

    elec_grid = extract_electricity_network(n)

    n_balancing_capa = extract_balancing_data("optimal_capacity", n)
    n_balancing_supply = extract_balancing_data("supply", n)
    n_power_capa = extract_inst_capa_elec_node(config, n, carriers_renamer)

    # # DataFrames to extract
    temporal_res_supply = extract_res_temporal_energy(config, n)
    n_res_pot = extract_res_potential(n)
    ACDC_grid, ACDC_countries, el_imp_exp = extract_transmission(n_ext)
    H2_grid, H2_countries, H2_imp_exp = extract_transmission(n_ext, carriers=["H2 pipeline", "H2 pipeline retrofitted"])
    gas_grid, gas_countries, gas_imp_exp = extract_transmission(n_ext, carriers=["gas pipeline", "gas pipeline new"])
    n_costs = extract_nodal_costs(config, n)
    marginal_prices, marginal_prices_t = extract_marginal_prices(n, carrier_list=["gas", "AC", "H2"])
    nodal_supply_energy = extract_nodal_supply_energy(config, n)
    temporal_supply_energy = extract_temporal_supply_energy(config, n, carriers_renamer=carriers_renamer)
    temporal_supply_energy_BE = extract_temporal_supply_energy(config, n, carriers_renamer=carriers_renamer,
                                                               country_aggregate="BE")
    temporal_supply_energy_FL = extract_temporal_supply_energy(config, n, carriers_renamer=carriers_renamer,
                                                               country_aggregate="FL")

    imp_exp = pd.concat([y.reset_index()
                        .set_index(["imports_exports", "countries", "year", "carriers"])
                         for y in [el_imp_exp, H2_imp_exp, gas_imp_exp]])

    buses = extract_geo_buses(n)

    # Figures to extract
    n_sto, n_h2 = extract_graphs(config, n, color_shift)

    # Define outputs and export them
    outputs = {
        # assets
        "res_potentials": n_res_pot,
        "power_production_countries": n_power_capa,
        "balancing_capacities_countries": n_balancing_capa,
        "balancing_supply_countries": n_balancing_supply,

        # networks
        "grid_capacities_countries": ACDC_countries,
        "H2_network_capacities_countries": H2_countries,
        "gas_network_capacities_countries": gas_countries,
        "elec_grid": elec_grid,

        # energy balance
        "imports_exports": imp_exp,
        "supply_energy_sectors": nodal_supply_energy,
        "temporal_supply_energy_sectors": temporal_supply_energy,
        "temporal_supply_energy_sectors_BE": temporal_supply_energy_BE,
        "temporal_supply_energy_sectors_FL": temporal_supply_energy_FL,

        # insights
        "costs_countries": n_costs,
        "marginal_prices_countries": marginal_prices,
        "marginal_prices_t_countries": marginal_prices_t,
        "temporal_res_supply": temporal_res_supply,

        # geo
        "buses": buses,
    }

    figures = {
        "storage_unit": n_sto,
        "h2_production": n_h2,
    }

    export_csvs_figures(config["path"]["csvs"], outputs, figures)

    return outputs, figures
