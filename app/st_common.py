from pathlib import Path

import streamlit as st
from PIL import Image

base_path = Path(__file__).parent
swoosh = Image.open(Path(base_path, "assets", "img", "swoosh.png"))
network_path = Path(base_path, "assets", "data")
scenario_dict = {
    "Central": {
        "path": "20240531/graph_extraction_st/central",
    },
    "Electrification": {
        "path": "20240531/graph_extraction_st/electrification",
    },
    "Molecule": {
        "path": "20240531/graph_extraction_st/molecule",
    },
    # "LSC": {
    #     "path": "20240531/graph_extraction_st/lsc",
    # },
}
CLIP_VALUE_TWH = 1e-1
GRAPH_AREA = {"EU27+TYNDP": "tot", "BE": "be", "EU27": "eu27"}
YEARS = ["2030", "2035", "2040", "2045", "2050"]


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
    return scenario
