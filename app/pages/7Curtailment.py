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
scenario = st_side_bar()

st.title("Curtailment and capacity factors")
st.markdown("Curtailment and capacity factors per technology, year and country. Electricity units are only.")


@st.cache_data(show_spinner="Retrieving data ...")
def get_data(scenario, type):
    df = (
        pd.read_csv(
            Path(network_path, scenario_dict[scenario]["path"], f"{type}.csv"),
            header=0
        )
    )
    return df


df_raw_curt = get_data(scenario, "curtailment")

all = "ENTSO-E area"
country = st.selectbox('Choose your country:', [all] + list(df_raw_curt["country"].unique()))

# %%
st.header("Curtailment per country")

unit = st.selectbox('Choose your unit:', ["GWh", "%"], )

if country != all:
    df_tech = df_raw_curt.query("country==@country").drop("country", axis=1).copy()
else:
    df_tech = (
        df_raw_curt
        .query("country != 'FL'")
        .groupby(by=["sector", "year"]).sum(numeric_only=True)
        .reset_index().copy()
    )
df_tech[["Curtailment", "Available"]] = df_tech[["Curtailment", "Available"]].div(1e3)  # GWh

if unit == "%":
    df_tech = df_tech.assign(unit=unit, Curtailment=lambda x: x["Curtailment"] / x["Available"] * 100)

df_tech = df_tech.pivot(index="sector", columns="year", values="Curtailment")

fig = px.bar(
    df_tech,
    title=f"Curtailment in {country} [{unit}]",
    barmode="group",
    text_auto=".2s" if unit != "%" else ".2f"
)

fig.update_layout(hovermode="x unified", legend_title_text="Years")
fig.update_yaxes(title_text=f'Curtailment [{unit}]')
fig.update_xaxes(title_text='Technologies')
fig.update_traces(hovertemplate="%{y:,.0f}" if unit != "%" else "%{y:,.2f}")
fig.update_layout(hovermode="x unified")

st.plotly_chart(
    fig
    , use_container_width=True
)

st.dataframe(df_tech.rename_axis(f"Curtailment [{unit}]")
             .style
             .format(precision=2, thousands=",", decimal='.'),
             use_container_width=True
             )

st.header("Curtailment on map")

techs = sorted(df_raw_curt.sector.unique())
tech = st.selectbox('Choose your technology (only for map):', techs, index=techs.index("Solar"))

if unit == "%":
    df_map = df_raw_curt.copy().assign(unit=unit, Curtailment=lambda x: x["Curtailment"] / x["Available"] * 100)
else:
    df_map = df_raw_curt.copy()
    
df_map = (
    df_map
    .query("sector==@tech")
    .groupby(["country", "year"]).sum(numeric_only=True)
    .join(get_buses())
    .reset_index()
)
fig_map = px.scatter_mapbox(
    df_map,
    lat="lat",
    lon="lon",
    size="Curtailment",
    mapbox_style="carto-positron",
    zoom=2.6,
    height=700,
    hover_name="country",
    animation_frame="year",
    title=f"Curtailment for {tech} [{unit}]",
    hover_data={f"Curtailment": ":.2f"}
)
fig_map.update_layout(sliders=[{"currentvalue": {"prefix": "Year: "}, "len": 0.8, "y": 0.07}])
fig_map.update_layout(updatemenus=[{"y": 0.07}])
st.plotly_chart(fig_map, use_container_width=True)

st.divider()

# %%
st.header("Capacity factors per country")

df_raw_cfs = get_data(scenario, "cfs")

if country != all:
    df_tech = df_raw_cfs.query("country==@country").drop("country", axis=1).copy()
else:
    df_tech = (
        df_raw_cfs
        .query("country != 'FL'")
        .groupby(by=["sector", "year"]).mean(numeric_only=True)
        .reset_index().copy()
    )
df_tech["cfs"] = df_tech["cfs"].mul(100)  # %

techs = st.multiselect('Choose your technologies:', list(df_tech.sector.unique()),
                       default=[i for i in ["Solar", "Offwind", "Onwind", "Gas", "Gas CC", "Gas CHP", "Gas CHP CC",
                                            "Oil", "Coal/Lignite", "Fuel Cell", "Nuclear", "Hydroelectricity"]
                                if i in list(df_tech.sector.unique())])
df_tech = df_tech.query("sector in @techs")

df_tech = df_tech.pivot(index="sector", columns="year", values="cfs")

fig = px.bar(
    df_tech,
    title=f"Curtailment in {country} [{unit}]",
    barmode="group",
    text_auto=".2s" if unit != "%" else ".2f"
)

fig.update_layout(hovermode="x unified", legend_title_text="Years")
fig.update_yaxes(title_text='Capcity factors [%]')
fig.update_xaxes(title_text='Technologies')
fig.update_traces(hovertemplate="%{y:,.2f}")
fig.update_layout(hovermode="x unified")

st.plotly_chart(
    fig
    , use_container_width=True
)

st.dataframe(df_tech.rename_axis("Capacity factor [%]")
             .style
             .format(precision=2, thousands=",", decimal='.'),
             use_container_width=True
             )

st.header("Capacity factors on map")

techs = sorted(df_raw_cfs.sector.unique())
tech = st.selectbox('Choose your technology (only for map):', techs, index=techs.index("Gas"))

df_map = (
    df_raw_cfs
    .query("sector==@tech")
    .groupby(["country", "year"]).sum(numeric_only=True)
    .join(get_buses())
    .reset_index()
)
fig_map = px.scatter_mapbox(
    df_map,
    lat="lat",
    lon="lon",
    size="cfs",
    mapbox_style="carto-positron",
    zoom=2.6,
    height=700,
    hover_name="country",
    animation_frame="year",
    title=f"Capacity factor for {tech} [%]",
    hover_data={f"cfs": ":.2f"}
)
fig_map.update_layout(sliders=[{"currentvalue": {"prefix": "Year: "}, "len": 0.8, "y": 0.07}])
fig_map.update_layout(updatemenus=[{"y": 0.07}])
st.plotly_chart(fig_map, use_container_width=True)
