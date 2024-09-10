from pathlib import Path

import pandas as pd
import plotly.express as px
import streamlit as st
from st_common import get_years
from st_common import network_path
from st_common import scenario_dict
from st_common import st_page_config
from st_common import st_side_bar

st_page_config(layout="wide")
scenario = st_side_bar()

YEARS = get_years(scenario)

st.title("Renewable production profiles per carrier")
st.markdown(
    "The RES production 3-hourly profiles for every carrier, year and subsector. You can zoom on these interactive graphs for specific time windows and you can also select/deselect various categories if you want.")


@st.cache_data(show_spinner="Retrieving data ...")
def get_df(scenario, year):
    return (
        pd.read_csv(
            Path(network_path, scenario_dict[scenario]["path"], "res_temporal_" + year + ".csv"),
            header=0,
        )
    )


# %%Cell Name
# should be able to 
# - Display per carrier
# - 3h load profile
# - eventually per country

col1, col2 = st.columns(2)
with col1:
    year = st.selectbox('Choose the year:', YEARS, index= len(YEARS) - 1)
df = get_df(scenario, year)

with col2:
    country = st.selectbox('Choose your country:', ['ENTSO-E area'] + list(df["country"].unique()))
if country != 'ENTSO-E area':
    df = df.query("country in @country")
else:
    df = df.query("country != 'FL'")
df = df.groupby(['carrier']).sum(numeric_only=True)

df_table = (
    (df.sum(axis=1) / 1e3  # TWh
     * 3)
    .to_frame(name=year)
    .style
    .format(precision=2, thousands=",", decimal='.')
)

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

df_table.index.name = "Production [TWh]"
st.dataframe(df_table, use_container_width=True)