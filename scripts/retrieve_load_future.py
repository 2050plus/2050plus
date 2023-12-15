# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2023 Climact for The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT

"""
This rule downloads yearly electric load data for each european country from the `2050 Pathways Explorer <https://pathwaysexplorer.climact.com>`_.
This rule downloads appropriate years using configuration file. Yearly electric load data are defined by per-country
scenarios.

**Releveant Settings**

.. code:: yaml

    planning_horizons:

    snapshots:
        start:

**Inputs**

- ``scenario_builder_tool_input.xlsx`` Projection scenario definitions for each region of Europe.
- ``resources/load.csv`` Hourly per-country load profiles.

**Outputs**

- ``data/patex`` Yearly electric load per-country.
"""

import json
import logging
import os
import re
from pathlib import Path

import pandas as pd
import requests

from _helpers import configure_logging

# ToDo Move this to yaml configuration
URL_TO_USE = "prod url"
# URL_TO_USE = "test url"

METRIC_MAP = pd.DataFrame([
    ["tra_energy-demand_domestic_electricity_BEV_freight_HDV[TWh]", "TR_hdv"],
    ["tra_energy-demand_domestic_electricity_PHEV_freight_HDV[TWh]", "TR_hdv"],
    ["tra_energy-demand_domestic_electricity_PHEVCE_freight_HDV[TWh]", "TR_hdv"],
    ["tra_energy-demand_domestic_electricity_BEV_freight_LDV[TWh]", "TR_ldv"],
    ["tra_energy-demand_domestic_electricity_PHEV_freight_LDV[TWh]", "TR_ldv"],
    ["tra_energy-demand_domestic_electricity_BEV_passenger_LDV[TWh]", "TR_cars"],
    ["tra_energy-demand_domestic_electricity_PHEV_passenger_LDV[TWh]", "TR_cars"],
    ["tra_energy-demand_domestic_electricity_BEV_passenger_bus[TWh]", "TR_bus"],
    ["tra_energy-demand_domestic_electricity_PHEV_passenger_bus[TWh]", "TR_bus"],
    ["tra_energy-demand_domestic_electricity_CEV_freight_rail[TWh]", "TR_rail"],
    ["bld_energy-demand_non-residential_heating_electricity[TWh]", "HE_ter_spa"],
    ["bld_energy-demand_non-residential_hotwater_electricity[TWh]", "HE_ter_wat"],
    ["bld_energy-demand_residential_heating_electricity[TWh]", "HE_res_spa"],
    ["bld_energy-demand_residential_hotwater_electricity[TWh]", 'HE_res_wat'],
    ["ind_energy-demand-by-carrier_electricity[TWh]", "IN_tot"],
    ["elc_elec-demand-by-energy-carrier-and-sector_electricity_refineries[TWh]", "SU_tot"],
    ["elc_elec-demand-by-energy-carrier-and-sector_electricity_losses[TWh]", "SU_tot"],
    ["elc_elec-demand-by-energy-carrier-and-sector_electricity_refineries[TWh]", "tot"],
    ["elc_elec-demand-by-energy-carrier-and-sector_electricity_losses[TWh]", "tot"],
    ["elc_elec-demand-by-energy-carrier-and-sector_electricity_agr[TWh]", "tot"],
    ["elc_elec-demand-by-energy-carrier-and-sector_electricity_bld[TWh]", "tot"],
    ["elc_elec-demand-by-energy-carrier-and-sector_electricity_ind[TWh]", "tot"],
    ["elc_elec-demand-by-energy-carrier-and-sector_electricity_tra[TWh]", "tot"],
    ["elc_elec-demand-by-energy-carrier-and-sector_electricity_hydrogen-for-sector[TWh]", "tot"],
    ["elc_elec-demand-by-energy-carrier-and-sector_electricity_hydrogen-for-power-prod[TWh]", "tot"],
    ["elc_elec-demand-by-energy-carrier-and-sector_electricity_efuels[TWh]", "tot"],
    ["elc_elec-demand-by-energy-carrier-and-sector_electricity_heat-CHP[TWh]",  "tot"],
    ["elc_elec-demand-by-energy-carrier-and-sector_electricity_heat-only[TWh]", "tot"]
], columns=["metric_id", "sector"])

TEMPLATE_REQUEST = """
{
  "levers":
  [LEVERS],
  "outputs": {"0": [VARIABLES]},
  "dimension": {"0": null},
  "aggregate": true
}
"""
API = {
    "prod url": {
        "db": "https://pathwaysexplorer.climact.com/api/v1.0/model_results",
        "model": "https://pathwaysexplorer.climact.com/model/v1.0/model_results"
    },
    "test url": {
        "db": "https://test.pathwaysexplorer.climact.com/api/v1.0/model_results",
        "model": "https://test.pathwaysexplorer.climact.com/model/v1.0/model_results"
    }
}

RX_SCENARIO = re.compile("^\[([A-z]+)\]\s(.*)$")


def __safely_to_int(f):
    """
    Cast a float value to int if possible. Otherwise, to float.

    :return : float or int
    """
    if isinstance(f, int):
        return f
    else:
        return int(f) if f.is_integer() else f


def load_scenario_builder():
    """
    Read scenarios from scenario builder in the same folder.
    This scenario builder should only contain the data to model EU27 as sum.

    :return df_scenarios: Dataframe of scenarios
    :return metrics: List of levers
    """
    return pd.read_excel(snakemake.input.scenario_builder, sheet_name="Levers", header=[0, 1]).iloc[:, 4:]


def get_eu27_countries():
    """
    Get the list of countries in EU 27
    """
    return ["AT", "BE", "BG", "HR", "CY", "CZ", "DK", "EE", "FI", "FR", "DE", "EL", "HU", "IE", "IT", "LV", "LT",
            "LU", "MT", "NL", "PL", "PT", "RO", "SK", "SI", "ES", "SE"]


def fill_scenario_list(df_scenarios):
    """
    Based on country list and scenario builder,
    create a complete list of scenario.

    :param df_scenarios: Defined scenarios and given exception
    :return scenarios_list: Complete list of scenarios as
        {
        [(master region 1, scenario name 1)] : {[country 1] : [levers 1],[country 2] : [levers 2]},
        [(master region 2, scenario name 2] : {[country 1] : [levers 3],[country 2] : [levers 4]}
        }
    """
    scenarios_short_list = [c for c in df_scenarios.columns if not RX_SCENARIO.match(c[1])]

    scenarios_dict = {}
    for master_region, name in scenarios_short_list:
        scenarios_dict_by_scenario = {}
        if master_region == "EU27 as sum":
            for eu27_c in get_eu27_countries():
                exception_name = f"[{eu27_c}] {name}"
                if exception_name in df_scenarios.columns.levels[1]:
                    levers = df_scenarios.loc[:, (master_region, exception_name)]
                else:
                    levers = df_scenarios.loc[:, (master_region, name)]
                scenarios_dict_by_scenario[eu27_c] = levers
                scenarios_dict[(master_region, exception_name)] = {eu27_c: levers}
        else:
            levers = df_scenarios.loc[:, (master_region, name)]
            scenarios_dict_by_scenario[master_region] = levers
        scenarios_dict[(master_region, name)] = scenarios_dict_by_scenario
    return scenarios_dict


def get_results(scenarios_dict, variables_list):
    """
    Get results for each scenario through API request

    :param scenarios_dict: dictionary with scenarios as key
    :param template_request: template of a request to API
    :param variables_list: list of variables to get
    :return results: Dict (key : scenario name, value : json result)
    """
    # Load the necessary variables and format them
    variables = "[\n"
    for var in set(variables_list):
        variables += f'"{var}",\n'
    variables = variables.rsplit(",\n", 1)[0] + "\n  ]"

    # Get result for each scenario (through API request)
    results = {}

    # Add variable in request payload
    my_template_request = TEMPLATE_REQUEST.replace("[VARIABLES]", variables)

    # Loop on countries and their scenario and format them
    for (master_region, scenario_name), values in scenarios_dict.items():
        logging.debug(f"          Getting '{scenario_name}' for {master_region}")

        levers = "{\n"
        for c, l in values.items():
            logging.debug(f"                Getting '{c}'")
            levers += f'        "{c}":\n        {{\n            "PyClimact":\n            {{\n                ' \
                      f'"static_levers":\n                {str([__safely_to_int(li) for li in l])},\n                ' \
                      f'"dynamic_levers":\n                {{}}\n            }}\n        }},\n'
        levers = levers.rsplit(",\n", 1)[0] + "\n    }"

        # Add scenarios in request payload
        my_request = my_template_request.replace("[LEVERS]", levers)

        # Get result
        logging.debug(f"                Requesting data to API")
        result = requests.post(API[URL_TO_USE]['db'], json=json.loads(my_request)).text
        result = json.loads(result)
        if not result:
            result = requests.post(API[URL_TO_USE]['model'], json=json.loads(my_request)).text
            result = json.loads(result)
            if not result:
                raise ValueError("API sent an empty response")
            if not result["outputs"]["0"]:
                raise ValueError("API sent an empty response")
        results[(master_region, scenario_name)] = result
        logging.debug(f"                Data received")

    return results


def format_results(results):
    """
    Format received JSON into dataframe
    Change unit from TWh to MWh (* 1e6)
    Add non-MS data using proportions

    :param results: data in JSON format
    :return df_results: data in dataframe
    """
    df_results = []
    for (master_region, scenario), values in results.items():
        values = values['outputs']['0']
        timeAxis = values[0]['timeAxis']
        data = []
        names = []
        for i in values:
            data.append(i['data'])
            if RX_SCENARIO.match(scenario):
                region = RX_SCENARIO.match(scenario).group(1)
            else:
                region = master_region
            names.append((scenario, region, i['id'], i['title']))
        df_results.append(
            pd.DataFrame(data, index=pd.MultiIndex.from_tuples(names, names=['scenario', 'region', 'metric_id',
                                                                             'metric_name']),
                         columns=timeAxis))

    df_results = pd.concat(df_results).reset_index(drop=False)

    df_results = (
        df_results.set_index("metric_id")
        .join(METRIC_MAP.set_index("metric_id"))
        .set_index("sector").reset_index()
    )
    # Hypothesis : One unique scenario for each country
    df_results = df_results.groupby(by=["region", "sector"]).sum().reset_index()
    df_results["key"] = df_results["region"] + "_" + df_results["sector"]
    df_results = df_results.set_index("key")

    horizons = [pd.Timestamp(snakemake.config["snapshots"]["start"]).year] + \
               snakemake.config["scenario"]["planning_horizons"]

    df_results = df_results[horizons] * 1e6
    df_results.index = df_results.index.str.replace("EL", "GR")

    # Manual scaling for non-MS if missing
    dict_non_MS = {"AL": "GR",
                   "MK": "GR",
                   "NO": "SK",
                   "CH": "FR",
                   "GB": "FR",
                   "BA": "HR",
                   "KV": "HU",
                   "RS": "HU",
                   "ME": "GR"}
    dict_non_MS = {k: v for k, v in dict_non_MS.items() if not any(df_results.index.str.startswith(k))}
    historical_load_h = pd.read_csv(snakemake.input.load_hourly, index_col="utc_timestamp", parse_dates=True)
    missing_load = []
    for non_MS, MS in dict_non_MS.items():
        df = df_results[df_results.index.str.startswith(MS)]
        df *= historical_load_h[non_MS].sum() / historical_load_h[MS].sum()
        df.index = df.index.str.replace(MS, non_MS)
        missing_load.append(df)
    df_results = pd.concat([df_results] + missing_load)

    return df_results


def write_files(df_results):
    """
    Export of data
    """
    df_results = df_results.T
    df_results.index = df_results.index.set_names("year")

    path = Path(snakemake.output[0]).parent
    for y in df_results.index:
        try:
            df_results.loc[[y]].to_csv(Path(path, f"patex_{y}.csv"))
        except OSError:
            os.makedirs(path)
            df_results.loc[[y]].to_csv(Path(path, f"patex_{y}.csv"))


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake("retrieve_load_future")

    configure_logging(snakemake)

    # Load configuration
    logging.info("Loading configuration")
    df_scenarios = load_scenario_builder()
    scenarios_dict = fill_scenario_list(df_scenarios)

    # Getting data from API
    logging.info("Getting data from API")
    results = get_results(scenarios_dict, METRIC_MAP["metric_id"])

    # Formatting data
    logging.info("Formatting data")
    df_results = format_results(results)

    # Writing data
    logging.info("Writing data")
    write_files(df_results)
