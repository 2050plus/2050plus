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

st.title("Loads per carrier")
st.markdown("The total yearly load per energy carrier, year, country and subsector.")


@st.cache_data(show_spinner="Retrieving data ...")
def get_data(scenario):
    df = (
        pd.read_csv(
            Path(network_path, scenario_dict[scenario]["path"], "graph_extraction_st", "supply_energy_df.csv"),
            header=0
        )
    )
    return df


# %%
df = get_data(scenario)

all = "EU27 + TYNDP"
col1, col2 = st.columns(2)
with col1:
    country = st.selectbox('Choose your country:', [all] + list(df["node"].unique()))
if country != all:
    df = df.query("node==@country").drop("node", axis=1)
else:
    df = df.groupby(by=["sector", "carrier"]).sum(numeric_only=True).reset_index()
with col2:
    carrier = st.selectbox('Choose your carrier:', df["carrier"].unique(), index=1)
df = df.query("carrier==@carrier").drop("carrier", axis=1)

df = df.groupby(by="sector").sum()

df_tot = pd.DataFrame(df.sum().rename("Total")).T
df = pd.concat([df, df_tot])

fig = px.bar(
    df,
    title=f"Load in {country} for {carrier} [TWh]",
    barmode="group",
    text_auto=".2s"
)

fig.update_traces(hovertemplate="%{y:,.0f}")
fig.update_layout(hovermode="x unified")
fig.update_yaxes(title_text='Consumption [TWh]')
fig.update_xaxes(title_text='Sectors')
fig.update_layout(legend_title_text='Technologies')

st.plotly_chart(
    fig
    , use_container_width=True
)

st.subheader(f"Annual load per sector for {carrier}")
st.table(df
         .rename(mapper=lambda x: x + " [TWh]", axis=1)
         .style
         .format(precision=2, thousands=",", decimal='.')
         )
