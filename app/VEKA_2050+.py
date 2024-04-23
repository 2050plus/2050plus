import streamlit as st

from st_common import st_page_config
from st_common import st_side_bar

st_page_config(layout="wide")
st_side_bar()

st.title("VEKA 2050+ results explorers")

st.markdown(
    """
    This notebook allows you to the explore the results of `VEKA 2050+` using a set of scenarios. Please refer to documentation for details.
    
    The following informations are available:
    - **Loads**: The total yearly load per energy carrier, year, country and subsector.
    - **Capacities**: The power production capacities installed per country, technologies and year.
    - **RES potentials**: The total potential for RES power production capacity considered by the model for the various technologies and countries.
    - **Consumption Profiles**: The load 3-hourly profiles for every carrier, year and subsector. This data is currently shown at system level (EU27+TYNDP) due to the very large quantity of data that needs to be handled for every country in the system.
    - **Production Profiles**: The same visualization as consumption for the supply side: 3-hourly production profiles by carrier and year.
    - **RES Profiles**: The same visualization as consumption for the supply side: 3-hourly production RES profiles by carrier and year. All countries are available.
    - **Imports & Exports**: Energy imports and exports between countries in the system, for all carriers, countries and years.
    - **Balancing**: (Under construction)
    - **Costs**: (Under construction)
    
    The following scenarios are available:
    - **1. Average**: Scenario *average* presented to the Steering Committee on 31/01/2024.
    - **1. Electrification**: Scenario *electrification* presented to the Steering Committee on 31/01/2024.
    - **2. Electrification**: Scenario *electrification* presented to the Steering Commitee on 29/04/2024.
""")
