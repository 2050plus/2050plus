import streamlit as st

from st_common import st_page_config
from st_common import st_side_bar

st_page_config(layout="wide")
st_side_bar()

st.title("Energie Studie 2050+ results explorers")

st.markdown(
    """
    This notebook allows you to the explore the results of `Energie Studie 2050+` using a set of scenarios. Please refer to documentation for further [details](https://2050plus-climact.readthedocs.io/en/latest/).
    
    The following informations are available:
    - **Loads**: The total yearly load per energy carrier, year, country and subsector.
    - **Capacities**: The power production capacities installed per country, technologies and year.
    - **RES potentials**: The total potential for RES power production capacity considered by the model for the various technologies and countries.
    - **Consumption**: The yearly load and 3-hourly profiles for every carrier, year and subsector. This data is currently shown at system level (ENTSO-E area), Belgium (BE), Flanders (FL) due to the very large quantity of data that needs to be handled for every country in the system.
    - **Production**: The same visualization as consumption for the supply side: yearly production and 3-hourly production profiles by carrier, year and subsector.
    - **Profiles**: The same visualization of 3-hourly profiles is presented. This view enables the comparison of profiles and leads to a better understanding of the dynamics.
    - **RES Profiles**: The same visualization as consumption for the supply side: 3-hourly production RES profiles by carrier and year. All countries are available.
    - **Curtailment**: The curtailment and capacity factors per technology, year and country.
    - **Imports & Exports**: Energy imports and exports between countries in the system, for all carriers, countries and years.
    - **Balancing**: (Under construction)
    - **Costs**: The total yearly cost of the modelled system per segments and year. 
    
    The following scenarios are available:
    - Mix, including the following sensitivities:
        - Nuclear extension
        - Nuclear cost
        - Pure optimisation
    - Electrification, including the following sensitivity:
        - Storage cost
    - Molecules, including the following sensitivity:
        - Molecules import
    - LSC (Least Structural Changes)

    \*The graphs show 2023 as the historical reference point. These data are not historical data as such, but rather a proxy created via PyPSA-Eur to give orders of magnitudes reference points.
""")
