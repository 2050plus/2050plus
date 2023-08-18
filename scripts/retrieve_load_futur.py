# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2023 Climact for The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT

"""
Retrieve futur loads from 2050 Pathways Explorer (Climact) model.
Parse scenario builder and compute request to API.
"""

import json
import logging
import re
from pathlib import Path

import pandas as pd
import requests

from _helpers import configure_logging

# ToDo Move this to yaml configuration
URL_TO_USE = "prod url"
# URL_TO_USE = "test url"

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

# GLOBAL VARS
RX_SCENARIO_EU27_SUM = re.compile("^\[([A-z]+)\]\s(.*)$")


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
    return pd.read_excel(snakemake.input.scenario_builder, sheet_name="Levers", header=[0, 1]).iloc[:, 5:], \
        list(pd.read_excel(snakemake.input.scenario_builder, sheet_name="Variables")["metrics"])


def get_eu27_countries():
    """
    Get the list of countries in EU 27
    """
    return ["AT", "BE", "BG", "HR", "CY", "CZ", "DK", "EE", "FI", "FR", "DE", "EL", "HU", "IE", "IT", "LV", "LT",
            "LU", "MT", "NL", "PL", "PT", "RO", "SK", "SI", "ES", "SE"]


def fill_scenario_list(df_scenarios, eu27_countries):
    """
    Based on country list and scenario builder,
    create a complete list of scenario.

    :param df_scenarios: Defined scenarios and given exception
    :param eu27_countries: List of countries in EU27
    :return scenarios_list: Complete list of scenarios as
        {
        [(master region 1, scenario name 1)] : {[country 1] : [levers 1],[country 2] : [levers 2]},
        [(master region 2, scenario name 2] : {[country 1] : [levers 3],[country 2] : [levers 4]}
        }
    """
    scenarios_short_list = [c for c in df_scenarios.columns if not RX_SCENARIO_EU27_SUM.match(c[1])]

    scenarios_dict = {}
    for master_region, name in scenarios_short_list:
        scenarios_dict_by_scenario = {}
        if master_region == "EU27 as sum":
            for eu27_c in eu27_countries:
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
    for var in variables_list:
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
            levers += f'        "{c}":\n        {{\n            "PyClimact":\n            {{\n                "static_levers":\n                {str([__safely_to_int(li) for li in l])},\n                "dynamic_levers":\n                {{}}\n            }}\n        }},\n'
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
            if RX_SCENARIO_EU27_SUM.match(scenario):
                region = RX_SCENARIO_EU27_SUM.match(scenario).group(1)
            else:
                region = master_region
            names.append((scenario, region, i['id'], i['title']))
        df_results.append(
            pd.DataFrame(data, index=pd.MultiIndex.from_tuples(names, names=['scenario', 'region', 'metric_id',
                                                                             'metric_name']),
                         columns=timeAxis))
    return pd.concat(df_results).reset_index(drop=False)


def write_files(df_results):
    """
    Temporary export for PyPSA
    """
    horizons = [2013, 2030, 2040, 2050]

    metric_map = {
        "tra_energy-demand_domestic_electricity_BEV_freight_HDV[TWh]": "TR_hdv",
        "tra_energy-demand_domestic_electricity_PHEV_freight_HDV[TWh]": "TR_hdv",
        "tra_energy-demand_domestic_electricity_PHEVCE_freight_HDV[TWh]": "TR_hdv",
        "tra_energy-demand_domestic_electricity_BEV_freight_LDV[TWh]": "TR_ldv",
        "tra_energy-demand_domestic_electricity_PHEV_freight_LDV[TWh]": "TR_ldv",
        "tra_energy-demand_domestic_electricity_BEV_passenger_LDV[TWh]": "TR_cars",
        "tra_energy-demand_domestic_electricity_PHEV_passenger_LDV[TWh]": "TR_cars",
        "tra_energy-demand_domestic_electricity_BEV_passenger_bus[TWh]": "TR_bus",
        "tra_energy-demand_domestic_electricity_PHEV_passenger_bus[TWh]": "TR_bus",
        "tra_energy-demand_domestic_electricity_CEV_freight_rail[TWh]": "TR_rail",
        "tra_energy-demand_domestic_electricity_CEV_passenger_rail[TWh]": "TR_rail"
    }

    df_results = df_results[df_results["metric_id"].isin(metric_map.keys())]
    df_results["metric_id"] = df_results["metric_id"].replace(metric_map)
    # Hypothesis : One unique scenario for each country
    df_results = df_results.groupby(by=["region", "metric_id"]).sum().reset_index()
    df_results["key"] = df_results["region"] + "_" + df_results["metric_id"]

    for y in horizons:
        df_results_y = df_results.set_index("key")[[y]] * 1e6
        df_results_y.T.to_csv(Path(snakemake.output[0], f"patex_{y}.csv"))


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake("retrieve_load_futur")

    configure_logging(snakemake)

    # Load configuration
    logging.info("Loading configuration")
    df_scenarios, metrics = load_scenario_builder()
    eu27_countries = get_eu27_countries()
    scenarios_dict = fill_scenario_list(df_scenarios, eu27_countries)

    # Getting data from API
    logging.info("Getting data from API")
    results = get_results(scenarios_dict, metrics)

    # Formatting data
    logging.info("Formatting data")
    df_results = format_results(results)

    # Writing data
    logging.info("Writing data")
    write_files(df_results)
