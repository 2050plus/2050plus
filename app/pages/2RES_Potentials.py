from pathlib import Path

import pandas as pd
import plotly.express as px
import streamlit as st
from st_common import network_path
from st_common import get_buses
from st_common import scenario_dict
from st_common import st_page_config
from st_common import st_side_bar

st_page_config(layout="wide")
scenario, compare = st_side_bar(show_compare=False)

st.title("Renewable production potentials")
st.markdown(
    "The total potential for RES power production capacity considered by the model for the various technologies and countries.")


@st.cache_data(show_spinner="Retrieving data ...")
def get_data(scenario):
    return (
        pd.read_csv(
            Path(network_path, scenario_dict[scenario]["path"], "res_potentials.csv"),
            header=0,
        )
    )


# %%
data = get_data(scenario)
df = data.copy()

st.header("Potentials per carrier")

df = df.groupby("country").sum()
carrier = st.multiselect('Choose your carrier:', list(df.columns.unique()), default=list(df.columns.unique())[2:4])
df = df.loc[:, carrier]
carrier_list = ' & '.join(list(map(str.capitalize, carrier)))
df.index.name = "Potential [GW]"

df_map = (
    df
    .rename_axis(index={"Potential [GW]": "country"})
    .sum(axis=1)
    .to_frame(name="Potential [GW]")
    .join(get_buses())
    .reset_index()
)

fig_map = px.scatter_mapbox(
    df_map,
    lat="lat",
    lon="lon",
    size="Potential [GW]",
    mapbox_style="carto-positron",
    zoom=2.6,
    height=700,
    hover_name="country",
    title=f"Cumulated potentials ({carrier_list}) [GW]",
    hover_data={"Potential [GW]": ":.2f"}
)
st.plotly_chart(fig_map, use_container_width=True)

fig = px.bar(
    df,
    title=f"{carrier_list} potentials [GW]",
    text_auto=".2s"
)

fig.update_yaxes(title_text='Potential [GW]')
fig.update_xaxes(title_text='Countries')
fig.update_traces(hovertemplate="%{y:,.0f}")
fig.update_layout(hovermode="x unified",
                  legend_title_text='Technologies')

st.plotly_chart(
    fig
    , use_container_width=True
)

st.dataframe(df
            .rename(mapper=lambda x: x.capitalize() + " [GW]", axis=1)
             .style
             .format(precision=2, thousands=",", decimal='.'),
             use_container_width=True
             )

st.divider()

st.header("Potentials per country")
df_tab = data.copy()
country = st.selectbox('Choose your country:', ['ENTSO-E area'] + list(df_tab.country.unique()))
if country != 'ENTSO-E area':
    df_tab = df_tab.set_index('country').loc[country]
else:
    df_tab = (
        df_tab
        .query("country != 'FL'")
        .sum(numeric_only=True)
    )

df_tab = (
    df_tab
    .rename(f"Potential for {country}")
    .to_frame().T
    .rename(mapper=lambda x: x.capitalize() + " [GW]", axis=1)
)

df_tab = (df_tab
          .style
          .format(precision=2, thousands=",", decimal='.')
          )

st.dataframe(df_tab, use_container_width=True)
