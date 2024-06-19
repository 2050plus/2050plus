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
from scripts.make_summary import assign_carriers
from scripts.make_summary import assign_locations
from scripts.make_summary import calculate_nodal_capacities
from scripts.make_summary import calculate_nodal_supply_energy
from scripts.plot_power_network import rename_techs_tyndp

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
def extract_loads(n):
    profiles = {}
    for y, ni in n.items():
        loads_t = ni.loads_t.p.T
        loads_t.index.names = ["Load"]
        loads_t["country"] = ni.buses.loc[ni.loads.loc[loads_t.index].bus].country.values
        loads_t.reset_index(inplace=True)
        loads_t["Load"].mask(loads_t["Load"].str.contains("NH3"), "NH3 for sectors", inplace=True)
        loads_t["Load"].mask(loads_t["Load"].str.contains("H2"), "H2 for sectors", inplace=True)
        loads_t["Load"].where(loads_t["Load"].str.contains("sectors"), "Electricity demand for sectors", inplace=True)

        loads_t = loads_t.groupby(["country", "Load"]).sum()
        loads_t.insert(0, column="Annual sum [TWh]", value=loads_t.sum(axis=1) / 1e6 * 8760 / len(ni.snapshots))
        profiles[y] = loads_t

    df = pd.concat(profiles, names=["Years"])
    df.insert(0, column="units", value="MW_e")
    df.loc[(slice(None), slice(None), "H2 for sectors"), "units"] = "MW_lhv,h2"
    df.loc[(slice(None), slice(None), "NH3 for sectors"), "units"] = "MW_lhv,nh3"
    return df


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


# %% Plot functions

def plot_series(network, carrier="AC", name="test", load_only=True, supply_only=False, path=None,
                _stop="-02-01", _start="-01-01", year="2013", colors=None,
                save=True, return_data=False, regionalized=False):
    n = network.copy()
    assign_locations(n)
    assign_carriers(n)
    assert not (load_only and supply_only), "Cannot have both supply-only and load-only modes"

    if carrier == "electricity":
        buses = n.buses.query("'AC' in carrier or index.str.contains('low voltage')").index
    else:
        buses = n.buses.index[n.buses.carrier.str.contains(carrier)]

    supply = pd.DataFrame(index=n.snapshots)
    for c in n.iterate_components(n.branch_components):
        n_port = 3 if c.name == "Link" else 2
        for i in range(n_port):
            supply = pd.concat(
                (
                    supply,
                    (-1)
                    * c.pnl["p" + str(i)]
                      .loc[:, c.df.index[c.df["bus" + str(i)].isin(buses)]]
                      .T.groupby([c.df.carrier, c.df["bus" + str(i)]])
                      .sum().T,
                ),
                axis=1,
            )

    for c in n.iterate_components(n.one_port_components):
        comps = c.df.index[c.df.bus.isin(buses)]
        supply = pd.concat(
            (
                supply,
                ((c.pnl["p"].loc[:, comps]).multiply(c.df.loc[comps, "sign"]))
                .T.groupby([c.df.carrier, c.df.bus])
                .sum().T,
            ),
            axis=1,
        )

    # Map buses to countries in MultiIndex
    supply.index = pd.to_datetime(pd.DatetimeIndex(supply.index.values, name='snapshots').strftime(f'{year}-%m-%d-%H'))
    supply.columns = pd.MultiIndex.from_tuples(
        [(x, n.buses.loc[y, "country"]) for (x, y) in supply.columns.values],
        names=["carrier", "country"])
    supply = (supply.T.groupby(level=[0, 1])
              .sum().T
              .stack(level=1, future_stack=True)
              .sort_index(level=[1, 0]))

    # Group by carriers if not regionalized
    if not regionalized:
        supply = supply.groupby(level=0).sum()

    supply = supply.T.groupby(rename_techs_tyndp).sum().T

    both = supply.columns[(supply < 0.0).any() & (supply > 0.0).any()]

    positive_supply = supply[both]
    negative_supply = supply[both]

    positive_supply[positive_supply < 0.0] = 0.0
    negative_supply[negative_supply > 0.0] = 0.0

    supply[both] = positive_supply

    suffix = " charging"

    negative_supply.columns = negative_supply.columns + suffix

    supply = pd.concat((supply, negative_supply), axis=1)

    if supply_only:
        supply = supply.loc[:, ~(supply <= 0).all(axis=0)]

    if load_only:
        supply = - supply.loc[:, (supply <= 0).all(axis=0)]

    # 14-21.2 for flaute
    # 19-26.1 for flaute

    start = year + _start
    stop = year + _stop

    threshold = 10e3

    to_drop = supply.columns[(abs(supply) < threshold).all()]

    if len(to_drop) != 0:
        logger.info(f"Dropping {to_drop.tolist()} from supply")
        supply.drop(columns=to_drop, inplace=True)

    if not (return_data):
        supply.index.name = None

    supply = supply / 1e3

    supply.rename(
        columns={"electricity": "electric demand", "heat": "heat demand"}, inplace=True
    )
    supply.columns = supply.columns.str.replace("residential ", "")
    supply.columns = supply.columns.str.replace("services ", "")
    supply.columns = supply.columns.str.replace("urban decentral ", "decentral ")

    supply = supply.T.groupby(supply.columns).sum().T
    new_columns = (supply.std() / supply.mean()).sort_values().index

    if return_data:
        return supply[new_columns]

    fig, ax = plt.subplots()
    fig.set_size_inches((8, 5))

    if colors:
        (
            supply.loc[start:stop, new_columns].plot(
                ax=ax,
                kind="area",
                stacked=True,
                linewidth=0.0,
                color=[
                    colors[i.replace(suffix, "")]
                    for i in new_columns
                ],
            )
        )
        loads = n.loads_t.p[n.loads.query('bus in @buses').index].sum(axis=1) / 1e3
        loads.index = pd.to_datetime(pd.DatetimeIndex(loads.index.values).strftime(f'{year}-%m-%d-%H'))
        loads[start:stop].plot(ax=ax, color="k", linestyle='-')

    handles, labels = ax.get_legend_handles_labels()

    handles.reverse()
    labels.reverse()

    new_handles = []
    new_labels = []

    for i, item in enumerate(labels):
        if "charging" not in item:
            new_handles.append(handles[i])
            new_labels.append(labels[i])

    ax.legend(new_handles, new_labels, ncol=3, loc="upper left", frameon=False)
    ax.set_xlim([start, stop])
    if load_only or supply_only:
        ax.set_ylim([supply.sum().min() * 2 / 1e3, supply.sum().max() * 2 / 1e3])
    else:
        ax.set_ylim([-1000, 1900])
    ax.grid(True)
    ax.set_ylabel("Power [GW]")
    fig.tight_layout()

    if save:
        if path:
            fig.savefig(path, transparent=True)
    else:
        return fig


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
def extract_res_statistics(n):
    df = []
    for y, ni in n.items():
        res = ni.generators.copy().query("carrier in @RES")
        res_t = ni.generators_t.p[res.index]
        res["year"] = y
        cf = (res_t / (res["p_nom_opt"])).mean()
        p_tot = res_t.sum() * ni.snapshot_weightings.generators.mean()
        opex = p_tot * ni.generators.marginal_cost
        capex = res.capital_cost * res.p_nom_opt

        res.loc[cf.index, "cf"] = cf
        res.loc[cf.index, "p_tot"] = p_tot
        res.loc[cf.index, "opex"] = opex
        res.loc[cf.index, "capex"] = capex
        res.loc[cf.index, "totex"] = res.loc[cf.index,
        "capex"] + res.loc[cf.index, "opex"]

        LCOE = res.groupby(by=["carrier"]).totex.sum() / \
               res.groupby(by=["carrier"]).p_tot.sum()
        LCOE = LCOE.rename("carrier").to_frame().rename(
            columns={"carrier": "LCOE"})
        res = res.reset_index().merge(LCOE, on="carrier").set_index(["Generator", "year"])
        res = res.loc[:,
              ["carrier", "bus", "capital_cost", "marginal_cost", "p_nom_opt", "build_year", "p_nom", "cf", "p_tot",
               "p_nom_max", "LCOE", "opex", "capex", "totex"]]
        res = res.sort_values(by="carrier")
        df.append(res)

    return pd.concat(df)


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


def extract_country_capacities(config, n):
    # TODO : dictionnary is useless here
    df = {}
    df["nodal_capacities"] = pd.DataFrame(columns=config["scenario"]["planning_horizons"], dtype=float)

    # Duplicates scripts.make_summary.calculate_nodal_capacites
    for y, ni in n.items():
        df["nodal_capacities"] = calculate_nodal_capacities(ni, y, df["nodal_capacities"],
                                                            _opt_name={"Store": "e", "Line": "s", "Transformer": "s",
                                                                       "Link": "p_carrier"})

    df_capa = (df["nodal_capacities"]
               .rename(RENAMER)
               .reset_index()
               .rename(columns={"level_0": "unit_type",
                                "level_1": "node",
                                "level_2": "carrier"}))

    df_capa.node = df_capa.node.apply(lambda x: x[:2])

    # add extraction and storage suffixes
    dico = {"generators": "_extraction", "stores": "_stores"}
    for d, suffix in dico.items():
        to_modify = df_capa.query(
            "unit_type in [@d] and carrier in ['gas', 'oil', 'coal/lignite', 'uranium', 'solid biomass']").index
        df_capa.loc[to_modify, ["carrier"]] += suffix

    df_capa = df_capa.groupby(["unit_type", "node", "carrier"]).sum().reset_index(["carrier", "unit_type"])

    df_capa = df_capa.drop(columns="unit_type").groupby(["node", "carrier"]).sum() / 1e3

    df_capa.loc[(slice(None), "Haber-Bosch"), :] *= 4.415385
    df_capa["units"] = "GW_e"
    df_capa.loc[(slice(None), ["Haber-Bosch", "ammonia cracker"]), "units"] = "GW_lhv,nh3"
    df_capa.loc[(slice(None), ["Sabatier"]), "units"] = "GW_lhv,h2"
    df_capa.loc[(slice(None), ["H2"]), "units"] = "GWh_lhv,h2"
    df_capa.loc[(slice(None), ["battery", "home battery"]), "units"] = "GWh_e"
    df_capa.loc[(slice(None), ["gas_extraction", "oil_extraction",
                               "coal/lignite_extraction", "uranium_extraction"]), "units"] = "GW_lhv"
    return df_capa


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

    return df, df_co, df_imp_exp


def extract_nodal_costs(config):
    # Todo : add handling of multiple runs
    df = (pd.read_csv(Path(config["path"]["results_path"], "csvs", "nodal_costs.csv"),
                      index_col=[0, 1, 2, 3],
                      skiprows=3,
                      header=0)
          .reset_index()
          .rename(columns={"planning_horizon": "type",
                           "level_1": "cost",
                           "level_2": "country",
                           "level_3": "carrier"})
          )
    df_capa = (pd.read_csv(Path(config["path"]["results_path"], "csvs", "nodal_capacities.csv"),
                           index_col=[0, 1, 2],
                           skiprows=3,
                           header=0)
               .reset_index()
               .rename(columns={"planning_horizon": "type",
                                "level_1": "country",
                                "level_2": "carrier"})
               )

    index_elec_capa = df_capa.query('carrier.str.contains("distribution")').index
    index_elec_cost = df.query('carrier.str.contains("distribution") and cost=="capital"').index
    capital_cost = (
            (df.loc[index_elec_cost, ['2030', '2040']].values /
             df_capa.loc[index_elec_capa, ['2030', '2040']].values)
            .mean() * config["sector"]["gas_distribution_grid_cost_factor"]
    )
    # gas boilers
    cond_str = '((carrier.str.contains("gas boiler") and not(carrier.str.contains("urban central"))) or carrier.str.contains("micro gas"))'
    df = df.set_index(['type', 'country', 'carrier', 'cost'])
    df_capa = df_capa.set_index(['type', 'country', 'carrier'])
    index_gas_cost = df.query(cond_str + 'and cost=="capital"').index
    index_gas_capa = df_capa.query(cond_str).index
    distri_gas = capital_cost * df_capa.loc[index_gas_capa, ['2030', '2040']]
    distri_gas['cost'] = 'capital'
    distri_gas = distri_gas.reset_index().set_index(['type', 'country', 'carrier', 'cost'])

    df.loc[index_gas_cost, ['2030', '2040']] = df.loc[index_gas_cost, ['2030', '2040']] - distri_gas[
        ['2030', '2040']].values

    assets_distri = df.loc[index_gas_cost].reset_index()
    assets_distri.carrier = 'gas distribution grid'
    assets_distri[['2030', '2040']] = distri_gas.values
    assets_distri = assets_distri.groupby(["type", "country", "carrier", "cost", ]).sum()
    df = pd.concat([df, assets_distri], axis=0).reset_index()

    df["country"] = df["country"].str[:2].fillna("EU")
    fuels = df.query(
        "carrier in ['gas','oil','coal','lignite','uranium'] and cost == 'marginal' and type == 'generators'").index
    biomass = df.query(
        "(carrier.str.contains('biomass') or carrier.str.contains('biogas')) and cost == 'marginal' and type == 'stores'").index
    df.loc[fuels.union(biomass), "cost"] = "fuel"
    df = df.set_index(["type", "cost", "country", "carrier"])
    df = df.fillna(0).groupby(["type", "cost", "country", "carrier"]).sum()
    df = df.loc[~df.apply(lambda x: x < 1e3).all(axis=1)]
    df.insert(0, column="units", value="Euro")
    return df


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
                prices[y] = price_y.mean().groupby(lambda x: x[:2]).mean()
                prices[f"{y}_std"] = price_y.std().groupby(lambda x: x[:2]).mean()
                marginal_t = price_y.T.groupby(lambda x: x[:2]).mean().rename_axis(index=["countries"])
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
    df = df.dropna().set_index(idx)

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


def extract_gas_phase_out(n, year):
    dimensions = ["country", "build_year"]
    n_cgt = (
        n[year].links[n[year].links.carrier.str.contains("CGT")]
        .merge(
            n[year].buses["country"].reset_index(),
            left_on="bus1",
            right_on="Bus",
            how="left"
        )
        .groupby(by=dimensions)
        ["p_carrier_nom_opt"]
        .sum(numeric_only=True)
        .reset_index()
    )
    n_cgt.loc[n_cgt["build_year"] < year, "build_year"] = "hist"
    n_cgt = (
                n_cgt
                .groupby(by=dimensions)
                .sum()
                .reset_index()
                .pivot(index="country", columns="build_year", values="p_carrier_nom_opt")
            ) / 1e3  # GW
    n_cgt["units"] = "GW_e"

    if year in n_cgt.columns:
        sorting = year
    else:
        sorting = "hist"

    n_cgt = n_cgt.sort_values(by=sorting, ascending=False)
    return n_cgt[n_cgt[sorting] >= 5].fillna(0)


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


def extract_series(config, n, carriers=["electricity"], load=False, supply=False):
    with plt.style.context(["ggplot"]):
        df = config["plotting"]
        plots = {}
        # TODO : add multiple carriers
        for carrier in carriers:
            for y, ni in n.items():
                with pd.option_context("mode.chained_assignment", None):
                    plots[y] = plot_series(ni, carrier=carrier, name=carrier, year=str(y),
                                           load_only=load, supply_only=supply, colors=df["tech_colors"],
                                           path=Path(config["path"]["csvs"], f"series_AC_{y}.png"), save=False)
    return plots


def extract_profiles(config, n, regionalized=False, supply=False, load=False):
    production_profiles = []
    assert (supply or load), "No data to return"

    for carrier in config["carriers_to_plot"]:
        df = []
        for y, ni in n.items():
            df.append(plot_series(ni, carrier=carrier, name=carrier, year=str(y),
                                  return_data=True, supply_only=supply,
                                  load_only=load, regionalized=regionalized))
        df = pd.concat(df)
        df["carrier"] = carrier
        production_profiles.append(df)
    production_profiles = (pd.concat(production_profiles)
                           .reset_index()
                           .set_index(["carrier", "snapshots", "country"] if regionalized else ["carrier", "snapshots"])
                           )
    return production_profiles


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


def transform_data(config, n, n_ext, color_shift=None):
    logger.info("Transforming data")

    carriers_renamer = {}
    carriers_renamer.update(HEAT_RENAMER)
    carriers_renamer.update(ELEC_RENAMER)

    elec_grid = extract_electricity_network(n)

    capa_country = extract_country_capacities(config, n_ext)
    n_balancing_capa = extract_balancing_data("optimal_capacity", n)
    n_balancing_supply = extract_balancing_data("supply", n)
    n_power_capa = extract_inst_capa_elec_node(config, n, carriers_renamer)

    # # DataFrames to extract
    temporal_res_supply = extract_res_temporal_energy(config, n)
    prod_profiles = extract_profiles(config, n, supply=True)
    load_profiles = extract_profiles(config, n, load=True)
    n_res_pot = extract_res_potential(n)
    res_stats = extract_res_statistics(n)
    capa_country = extract_country_capacities(config, n_ext)
    ACDC_grid, ACDC_countries, el_imp_exp = extract_transmission(n_ext)
    H2_grid, H2_countries, H2_imp_exp = extract_transmission(n_ext, carriers=["H2 pipeline", "H2 pipeline retrofitted"])
    gas_grid, gas_countries, gas_imp_exp = extract_transmission(n_ext, carriers=["gas pipeline", "gas pipeline new"])
    n_costs = extract_nodal_costs(config)
    marginal_prices, marginal_prices_t = extract_marginal_prices(n, carrier_list=["gas", "AC", "H2"])
    nodal_supply_energy = extract_nodal_supply_energy(config, n)
    temporal_supply_energy = extract_temporal_supply_energy(config, n, carriers_renamer=carriers_renamer)
    temporal_supply_energy_BE = extract_temporal_supply_energy(config, n, carriers_renamer=carriers_renamer,
                                                               country_aggregate="BE")
    temporal_supply_energy_FL = extract_temporal_supply_energy(config, n, carriers_renamer=carriers_renamer,
                                                               country_aggregate="FL")

    n_gas_out = extract_gas_phase_out(n, config["scenario"]["planning_horizons"][0])

    imp_exp = pd.concat([y.reset_index()
                        .set_index(["imports_exports", "countries", "year", "carriers"])
                         for y in [el_imp_exp, H2_imp_exp, gas_imp_exp]])

    # Figures to extract
    n_sto, n_h2 = extract_graphs(config, n, color_shift)
    series_production = extract_series(config, n, supply=True)
    series_consumption = extract_series(config, n, load=True)

    # Define outputs and export them
    outputs = {
        # assets
        "units_capacities_countries": capa_country,
        "gas_phase_out": n_gas_out,
        "res_potentials": n_res_pot,
        "power_production_countries": n_power_capa,
        "balancing_capacities_countries": n_balancing_capa,
        "balancing_supply_countries": n_balancing_supply,

        # networks
        "grid_capacities_countries": ACDC_countries,
        "H2_network_capacities_countries": H2_countries,
        "gas_network_capacities_countries": gas_countries,
        "grid_capacities": ACDC_grid,
        "H2_network_capacities": H2_grid,
        "gas_network_capacities": gas_grid,
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
        "res_statistics": res_stats,
        # "loads_profiles": n_loads,
        "generation_profiles": prod_profiles,
        "load_profiles": load_profiles,
        "temporal_res_supply": temporal_res_supply
    }

    figures = {
        "storage_unit": n_sto,
        "h2_production": n_h2,
        "series_production": series_production,
        "series_consumption": series_consumption,
    }

    export_csvs_figures(config["path"]["csvs"], outputs, figures)

    return outputs, figures
