from pathlib import Path

import pandas as pd
import plotly.express as px
import streamlit as st
from st_common import network_path
from st_common import scenario_dict
from st_common import st_page_config
from st_common import st_side_bar

st_page_config(layout="wide")
scenario = st_side_bar()

st.title("Renewable production per carrier")
st.markdown("The RES production 3-hourly profiles for every carrier, year and subsector. You can zoom on these interactive graphs for specific time windows and you can also select/deselect various categories if you want.")


@st.cache_data(show_spinner="Retrieving data ...")
def get_df(scenario, year):
    return (
        pd.read_csv(
            Path(network_path, scenario_dict[scenario]["path"], "graph_extraction_st", "res_temporal_" + year + ".csv"),
            header=0,
        )
    )


# %%Cell Name
# should be able to 
# - Display per carrier
# - 3h load profile
# - eventually per country
years = ['2030', '2040', '2050']

col1, col2 = st.columns(2)
with col1:
    year = st.selectbox('Choose the year:', years)
df = get_df(scenario, year)

with col2:
    country = st.selectbox('Choose your country:', ['EU27 + TYNDP'] + list(df["country"].unique()))
if country != 'EU27 + TYNDP':
    df = df.query("country in @country")
df = df.groupby(['carrier']).sum(numeric_only=True)

df_table = (
    (df.sum(axis=1) / 1e3  # TWh
     * 3)
    .rename(f"Annual production [TWh]")
    .to_frame()
    .style
    .format(precision=2, thousands=",", decimal='.')
)

st.subheader(f"Renewable annual production for {country}")
st.table(df_table)
st.subheader(f"Renewable production profiles for {country}")

carrier = st.selectbox('Choose your carrier:', list(df.index.unique()))
if carrier != 'all':
    df = df.query("carrier in @carrier")
df = df.sum(axis=0).rename(carrier).to_frame()

fig = px.area(
    df,
    title=f"{carrier} production profile for {country}  [GW]",
)
fig.update_traces(hovertemplate="%{y:,.0f}",
                  line=dict(width=0.1))
fig.update_layout(legend_traceorder="reversed",
                  hovermode="x unified",
                  legend_title_text='Technologies')
fig.update_yaxes(title_text='Production [GW]')
fig.update_xaxes(title_text='Timesteps')

st.plotly_chart(
    fig
    , use_container_width=True
)
