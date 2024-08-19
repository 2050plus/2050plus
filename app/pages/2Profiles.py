from pathlib import Path

import pandas as pd
import plotly.express as px
import plotly.subplots as sp
import streamlit as st
from st_common import PROFILES_AREA
from st_common import get_years
from st_common import network_path
from st_common import scenario_dict
from st_common import st_page_config
from st_common import st_side_bar

st_page_config(layout="wide")
scenario = st_side_bar()

YEARS = get_years(scenario)

st.title("Production and consumption profiles")
st.markdown(
    "The 3-hourly consumption and production profiles for every carrier, year and subsector. This data is currently shown at system level (ENTSO-E area), Belgium (BE) and Flanders (FL) due to the very large quantity of data that needs to be handled for every country in the system. \n\n"
    "You can zoom on these interactive graphs for specific time windows and you can also select/deselect various categories if you want. We suggest to maximize the graph to improve the experience.")


@st.cache_data(show_spinner="Retrieving data ...")
def get_data_supply(scenario, year, selected_area):
    area = selected_area if selected_area != "ENTSO-E area" else ''
    return (
        pd.read_csv(
            Path(network_path, scenario_dict[scenario]["path"],
                 f"supply_temporal_{area.lower()}_{year}.csv".replace('__', '_')),
            header=[1, 2]
        )
        .set_index(("carrier", "sector"))
    )


@st.cache_data(show_spinner="Retrieving data ...")
def get_data_load(scenario, year, selected_area):
    area = selected_area if selected_area != "ENTSO-E area" else ''
    return (
        pd.read_csv(
            Path(network_path, scenario_dict[scenario]["path"],
                 f"load_temporal_{area.lower()}_{year}.csv".replace('__', '_')),
            header=[1, 2]
        )
        .set_index(("carrier", "sector"))
    )


col1, col2 = st.columns(2)
with col1:
    selected_area = st.selectbox('Choose area :', PROFILES_AREA)
with col2:
    year = st.selectbox('Choose the year:', YEARS, index=len(YEARS) - 1)

prod = "Production"
cons = "Consumption"
types = [prod, cons]
col1, col2, col3, col4 = st.columns(4)
with col1:
    carrier1 = st.selectbox('Choose your top carrier:',
                            get_data_supply(scenario, year, selected_area).columns.unique(0).sort_values(),
                            index=get_data_supply(scenario, year, selected_area).columns.unique(0).sort_values()
                            .get_loc("Electricity"))
unit1 = "GW" if "co2" not in carrier1.lower() else "Mt"
with col2:
    type1 = st.selectbox('Choose your top type:', types, index=0)
with col3:
    carrier2 = st.selectbox('Choose your bottom carrier:',
                            get_data_load(scenario, year, selected_area).columns.unique(0).sort_values(),
                            index=get_data_load(scenario, year, selected_area).columns.unique(0).sort_values()
                            .get_loc("Electricity"),)
unit2 = "GW" if "co2" not in carrier2.lower() else "Mt"
with col4:
    type2 = st.selectbox('Choose your bottom type:', types, index=1)

if type1 == prod:
    data1 = get_data_supply(scenario, year, selected_area)
elif type1 == cons:
    data1 = get_data_load(scenario, year, selected_area)
if type2 == prod:
    data2 = get_data_supply(scenario, year, selected_area)
elif type2 == cons:
    data2 = get_data_load(scenario, year, selected_area)

df = data1[carrier1]
fig1 = px.area(df)

df = data2[carrier2]
fig2 = px.area(df)

# For as many traces that exist per Express figure, get the traces from each plot and store them in an array.
# This is essentially breaking down the Express fig into it's traces
figure1_traces = []
figure2_traces = []
for trace in range(len(fig1["data"])):
    figure1_traces.append(fig1["data"][trace])
for trace in range(len(fig2["data"])):
    figure2_traces.append(fig2["data"][trace])

# Create a 1x2 subplot
this_figure = sp.make_subplots(rows=2, cols=1, shared_xaxes=True, shared_yaxes=True,
                               subplot_titles=[
                                   f"System {type1.lower()} profile for {carrier1} [{unit1}]",
                                   f"System {type2.lower()} profile for {carrier2} [{unit2}]",
                               ])

# Get the Express fig broken down as traces and add the traces to the proper plot within in the subplot
for traces in figure1_traces:
    this_figure.add_trace(traces.update(legendgrouptitle=dict(text=type1), legendgroup="1"), row=1, col=1)
for traces in figure2_traces:
    this_figure.add_trace(traces.update(legendgrouptitle=dict(text=type2), legendgroup="2"), row=2, col=1)

this_figure.update_layout(hovermode="x unified",
                          legend_title_text='Technologies',
                          legend_tracegroupgap=10,)

this_figure.update_traces(hovertemplate="%{y:,.0f}",
                          line=dict(width=0.1),
                          row=1, col=1)
this_figure.update_yaxes(title_text=f'{type1} [{unit1}]',
                         row=1, col=1)

this_figure.update_traces(hovertemplate="%{y:,.0f}",
                          line=dict(width=0.1),
                          row=2, col=1)
this_figure.update_yaxes(title_text=f'{type2} [{unit2}]',
                         row=2, col=1)
this_figure.update_xaxes(title_text='Timesteps',
                         row=2, col=1)

st.plotly_chart(
    this_figure
    , use_container_width=True
)
