# -*- coding: utf-8 -*-
from pathlib import Path

import pandas as pd
import pycountry as pyc
import streamlit as st
from st_common import st_page_config
from st_common import st_side_bar
from st_common import scenario_dict
from st_common import network_path

st_page_config(layout="wide")
# TODO : check if scenario is relevant
scenario = st_side_bar()

st.title("FAQ")

with st.expander("**Acronyms**"):
    st.markdown("Here is the list of used acronyms:")
    acronyms = {
        "ror": "Run-of-the-river",
        "V2G": "Vehicle to grid",
        "BEV": "Battery Electric Vehicle",
        "EV": "Electric Vehicle",
        "PHS": "Pumped hydro storage",
        "TES": "Thermal energy storage",
        "H2": "Hydrogen",
        "NH3": "Ammonia",
        "CH4": "Methane",
        "CHP": "Combined Heat and Power",
        "SMR": "Steam Methane Reforming",

    }
    acronyms = {k.upper(): v.capitalize() for k, v in acronyms.items()}
    for key, value in sorted(acronyms.items()):
        st.write(f"**{key}**: {value}")

with st.expander("**Description of carriers**"):
    carriers = {
        "centralized heat": "Heat supplied by large-scale district heating networks in urban areas with dense heat"
                            "population",
        "decentralized heat": "Heat supplied to buildings not using district heating",
        "electricity": "Electricity",
        "service rural heat": "Heat supplied to services buildings in rural areas with low population density."
                              " Heat demand from agriculture sector is also included here.",
        "oil": "Liquid fuel (oil) of fossil or synthetic origin (Fischer-Tropsch)",
        "H2": "Hydrogen (Steam Methane Reforming or Electrolyser)",
        "solid biomass": "Solid biomass",
        "gas": "Methane of fossil or synthetic origin (Sabatier)",
        "residential rural heat": "Heat supplied to residential buildings in rural areas with low population density",
        "residential urban decentral heat": "Heat supplied to residential buildings in urban areas not using district "
                                            "heating",
        "services urban decentral heat": "Heat supplied to services buildings in urban areas not using district "
                                         "heating",
        "urban central heat": "Heat supplied by large-scale district heating networks in urban areas with dense heat"
                              "population",
        "methanol": "Methanol produced from hydrogen and carbon dioxide using electricity (Methanolisation). Historical"
                    " methanol demand is accounted as electricity (0.167 MWh/t) and methane (10.25 MWh/t) demand.",
        "ammonia": "Ammonia (Haber-Bosch)",

    }
    carriers = {k.capitalize(): v for k, v in carriers.items()}
    for key, value in sorted(carriers.items()):
        st.write(f"**{key}** {value}")

with st.expander("**Description of energy assets cost segments**"):
    st.markdown(
        "To simplify the cost structure of the graphs, outputs are aggregated using a mapping table. This table is presented here bellow. The *cost_segment* column is the one used to display the graphs.")
    st.dataframe(
        pd.read_csv(Path(network_path, Path(scenario_dict["Mix"]["path"]).parent, "cost_mapping.csv")),
        hide_index=True,
        use_container_width=True
    )
    st.markdown(
        f"**Net Imports** category focuses on net import costs per country and per carrier. Already taken into account in the Energy production \'fuel\' category.")

with st.expander("**Description of technologies**"):
    st.markdown("Here is a short description of the technologies used:")
    technos = {
        "BEV charger": "Electrical consumption of the EV batteries from the grid.",
        "H2 Fuel cell": "Technology that converts hydrogen gas and oxygen into electricity.",
        "PHS": "Method of storing energy by pumping water to an elevated reservoir when energy is abundant and "
               "releasing it through turbines to generate electricity when demand is high.",
        "V2G": "(Vehicule to grid) Electricity production from the EV batteries to the grid.",
        "Air heat pump": "Device that uses electricity to transfer heat from the outside air to provide heating and "
                         "cooling for buildings.",
        "Battery charger": "Device or system designed to charge stored electrical energy in batteries from the "
                           "grid. The device can be reversed to return some of the electricity to the grid.",
        "Ground heat pump": "Device that uses electricity to transfer heat from the ground to provide heating and "
                            "cooling for buildings.",
        "Home battery charger": "Device or system used to charge stored electrical energy into home-based "
                                "batteries. The device can be reversed and return some of the electricity to the "
                                "grid.",
        "Hydro": "(Hydropower - Hydroelectricity) System that converts the potential energy of water stored in a "
                 "dam into electricity using a turbine. Different from *run-of-the-river*.",
        "Water tanks charger": ": Device used to charge TES (Thermal Energy Storage), i.e. long-term storage of "
                               "heat with water.",
        "Battery": "Device designed to store electrical energy from the grid. The device can be reversed to return "
                   "some of the electricity to the grid.",
        "Coal": "Thermal power station which burns coal to generate electricity.",
        "Lignite": "Thermal power station which burns lignite to generate electricity.",
        "Gas": "Thermal power station which burns methane (gas) to generate electricity.",
        "Nuclear": "Thermal power station using a nuclear reactor as heat source to generate electricity via a steam "
                   "turbine.",
        "Nuclear (SMR)": "(Small modular reactor) Specific class of thermal power station using a nuclear reactor as heat source to generate "
                         "electricity via a steam turbine. They are designed to be built in a factory, then shipped to"
                         "operational sites.",
        "Offwind": "Form of renewable energy that harnesses the power of the wind to generate electricity. This "
                   "particular class operates at sea.",
        "Oil": "Thermal power station which burns oil to generate electricity.",
        "Onwind": "Form of renewable energy that harnesses the power of the wind to generate electricity. This "
                  "particular class works on land.",
        "Solar": "Form of renewable energy that harnesses the power of the sun radiation to generate electricity.",
        "Solid biomass": "Thermal power station which burns solid biomass to generate electricity.",
        "Solid biomass CHP": "Thermal power station which burns solid biomass to generate electricity and useful "
                             "thermal energy (cogeneration).",
        "Ror": "(Run-of-the-river) System that harness the constant flow of water from rivers, without reservoirs."
               " A turbine is used to generate electricity. Different from *hydro*.",
        "Solar thermal": "Form of renewable energy that harnesses the power of the sun radiation to generate heat.",

    }
    for key, value in sorted(technos.items()):
        st.write(f"**{key}**: {value}")

with st.expander("**Countries considered**"):
    countries = ['AL', 'AT', 'BA', 'BE', 'BG', 'CH', 'CZ', 'DE', 'DK', 'EE', 'ES', 'FI', 'FR', 'GB', 'GR', 'HR', 'HU',
                 'IE', 'IT', 'LT', 'LU', 'LV', 'ME', 'MK', 'NL', 'NO', 'PL', 'PT', 'RO', 'RS', 'SE', 'SI', 'SK']
    countries_name = [pyc.countries.get(alpha_2=c).name for c in countries]
    st.markdown("The countries considered in this model are the following:")
    for c in sorted(countries_name):
        st.markdown(f"* {c}")

with st.expander("**Mapping of sectors**"):
    st.write(
        "To simplify the graphs, outputs are aggregated using a mapping table. This table is presented here bellow. The *sector* column is the one used to display the graphs.")
    st.dataframe(
        pd.read_csv(Path(network_path, Path(scenario_dict["Mix"]["path"]).parent, "sector_mapping.csv")),
        hide_index=True,
        use_container_width=True
    )

with st.expander("**Known issues**"):
    st.write(
        "* *Decentral heat production* value (in tab *Consumption Profiles*) contains both values of *central* and *decentral heat production*.")
