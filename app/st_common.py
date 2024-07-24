from pathlib import Path

import pandas as pd
import streamlit as st
from PIL import Image

base_path = Path(__file__).parent
swoosh = Image.open(Path(base_path, "assets", "img", "swoosh.png"))
network_path = Path(base_path, "assets", "data")
scenario_dict = {
    "Central": {
        "path": "20240619/graph_extraction_st/central",
    },
    "*Central - Nuclear extension": {
        "path": "20240709/graph_extraction_st/nuc_extension",
    },
    "*Central - Nuclear cost": {
        "path": "20240709/graph_extraction_st/nuc_cost",
    },
    "*Central - Pure Optimisation": {
        "path": "20240709/graph_extraction_st/pure_opt",
    },
    "Electrification": {
        "path": "20240619/graph_extraction_st/electrification",
    },
    "*Electrification - Storage cost": {
        "path": "20240709/graph_extraction_st/storage_cost",
    },
    "Molecules": {
        "path": "20240619/graph_extraction_st/molecules",
    },
    "*Molecules - Molecules import": {
        "path": "20240709/graph_extraction_st/mol_import",
    },
    "LSC (Least Structural Changes)": {
        "path": "20240619/graph_extraction_st/lsc",
    },
}
CLIP_VALUE_TWH = 1e-1
COSTS_AREA = {"ENTSO-E area": "tot", "EU27": "eu27", "BE": "be", "FL": "fl", "DE": "de", "FR": "fr", "GB": "gb",
              "LU": "lu", "NL": "nl"}
PROFILES_AREA = ["ENTSO-E area", "BE", "FL"]
YEARS = ["2020", "2030", "2035", "2040", "2045", "2050"]


def st_page_config(layout=None):
    if layout is None:
        layout = "centered"
    st.set_page_config(
        page_title="VEKA 2050+",
        page_icon=swoosh,
        layout=layout,
        initial_sidebar_state="expanded",
        menu_items={
            "Get Help": "mailto:tgi@climact.com",
            "Report a bug": "mailto:tgi@climact.com",
            "About": "# Experimental Climact data explorer for PyPSA app",
        },
    )


@st.cache_data
def get_buses():
    return pd.read_csv(Path(network_path, scenario_dict["Central"]["path"], "buses.csv"), index_col=0)


# @st.cache_data
# def load_network(scenario, year):
#     return pypsa.Network(Path(network_path,
#                               scenario_dict[scenario]["path"],
#                               scenario_dict[scenario]["fn"].replace("YEAR", str(year))
#                               )
#                          )


def st_side_bar():
    with st.sidebar:
        scenario = st.selectbox(
            "Select your scenario",
            scenario_dict.keys(),
            index=0
        )
        st.markdown("\* Sensitivities and 2020 data are currently under active development")
    return scenario
