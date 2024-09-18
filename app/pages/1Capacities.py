from pathlib import Path

import pandas as pd
import plotly.express as px
import streamlit as st
from st_common import get_buses
from st_common import network_path
from st_common import scenario_dict
from st_common import st_page_config
from st_common import st_side_bar

st_page_config(layout="wide")
scenario, compare = st_side_bar()

st.title("Power production installed capacities")
st.markdown(
    "The power production capacities installed per country, technologies and year. Power production units are units able to supply electricity.")


@st.cache_data(show_spinner="Retrieving data ...")
def get_data(scenario):
    return (
        pd.read_csv(
            Path(network_path, scenario_dict[scenario]["path"], "power_capacities.csv"),
            header=0,
        )
    )


# %%
data = get_data(scenario)
if compare != '-':
    df_compare = get_data(compare)
    idx = ["country", "sector",]
    data = (
        (data.set_index(idx) - df_compare.set_index(idx))
        .reset_index()
    )
df = data.copy()

st.header("Installed capacities per country")

all = ['ENTSO-E area']
country = st.selectbox('Choose your country:', all + list(df.country.unique()))
if not ('ENTSO-E area' in country):
    df = df.query("country in @country")
else:
    df = df.query("country != 'FL'")

df = (
    df.drop(columns=['country'])
    .groupby(['sector'])
    .sum(numeric_only=True)
    .rename(index={"sector": "Technologies"})
)
df.index.name = "Capacity [GW]"

fig = px.bar(
    df,
    title="Power production installed capacities [GW]",
    barmode="group",
    text_auto=".2s"
)

fig.update_yaxes(title_text='Installed capacities [GW]')
fig.update_xaxes(title_text='Technologies')
fig.update_traces(hovertemplate="%{y:,.0f}")
fig.update_layout(hovermode="x unified",
                  legend_title_text='Years')

st.plotly_chart(
    fig
    , use_container_width=True
)

st.dataframe(df
             .style
             .format(precision=2, thousands=",", decimal='.'),
             use_container_width=True
             )

st.divider()

st.header("Split of capacities per country")

df_bar = data.copy()
df_bar_wind = (
    df_bar
    .query("sector.str.contains('wind')")
    .groupby(by="country").sum(numeric_only=True)
    .assign(sector="Wind (all)")
    .set_index("sector", append=True)
    .reset_index()
)
df_bar_vres = (
    df_bar
    .query("sector.str.contains('Solar|wind')")
    .groupby(by="country").sum(numeric_only=True)
    .assign(sector="Wind and Solar")
    .set_index("sector", append=True)
    .reset_index()
)
df_bar = pd.concat([df_bar, df_bar_wind, df_bar_vres])
technology = st.selectbox('Choose your technology:', list(df_bar.sector.sort_values().unique()))

df_bar = (df_bar
          .query("sector == @technology")
          .drop(columns=['sector'])
          .set_index('country')
          .rename_axis("Investment year")
          .fillna(0)
          )
df_bar.index.name = "Capacity [GW]"

df_map = (
    df_bar
    .join(get_buses())
    .reset_index()
    .rename(columns={"Capacity [GW]": "country"})
    .melt(id_vars=["country", "lat", "lon"], value_name="Capacity [GW]", var_name="year")
    .assign(absolute_size=lambda x: x["Capacity [GW]"].abs())
)
fig_map = px.scatter_mapbox(
    df_map,
    lat="lat",
    lon="lon",
    size="absolute_size",
    mapbox_style="carto-positron",
    zoom=2.6,
    height=700,
    hover_name="country",
    animation_frame="year",
    title=f"{technology.capitalize()} installed capacities [GW]",
    hover_data={"Capacity [GW]": ":.2f"}
)
fig_map.update_layout(sliders=[{"currentvalue": {"prefix": "Year: "}, "len": 0.8, "y": 0.07}])
fig_map.update_layout(updatemenus=[{"y": 0.07}])
st.plotly_chart(fig_map, use_container_width=True)

fig_bar = px.bar(
    df_bar,
    barmode="group",
)

fig_bar.update_yaxes(title_text='Installed capacities [GW]')
fig_bar.update_xaxes(title_text='Countries')
fig_bar.update_traces(hovertemplate="%{y:,.1f}", )
fig_bar.update_layout(hovermode="x unified",
                      legend_title_text='Years')
st.plotly_chart(fig_bar
                , use_container_width=True)

st.dataframe(df_bar.style.format(precision=2, thousands=",", decimal='.'), use_container_width=True)
