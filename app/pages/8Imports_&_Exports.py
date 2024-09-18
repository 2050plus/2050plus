from pathlib import Path

import numpy as np
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

st.title("Imports and exports per carrier")
st.markdown(
    "The energy imports and exports between countries in the system, for all carriers and countries. A negative value means that the area is exporting.")


@st.cache_data(show_spinner="Retrieving data ...")
def get_data(scenario):
    df = (
        pd.read_csv(
            Path(network_path, scenario_dict[scenario]["path"], "imports_exports.csv"),
            header=0
        )
    )
    return df


df = get_data(scenario)
if compare != '-':
    df_compare = get_data(compare)
    idx = ["imports_exports", "countries", "year", "carriers"]
    df = (
        (df.set_index(idx) - df_compare.set_index(idx))
        .reset_index()
    )


def query_imp_exp(df, carriers, country, imports_exports, year=None):
    if country != "BE":
        df_imp_exp = (
            df.query(""
                     "carriers == @carriers & "
                     "imports_exports == @imports_exports"
                     )
            .drop(["carriers", "imports_exports"], axis=1)
            .set_index(["year", "countries"])
            [country]
        )
        if year:
            df_imp_exp = df_imp_exp.loc[year]
    else:
        df_imp_exp = []
        regions = ["WL", "FL", "BX"]
        for c in regions:
            df_imp_exp.append(query_imp_exp(df, carriers, c, imports_exports, year))
        df_imp_exp = (
            pd.concat(df_imp_exp, axis=1)
            .query("countries not in @regions")
            .sum(axis=1)
            .rename("BE")
        )
        print("ok")
    return df_imp_exp


col1, col2 = st.columns(2)
with col1:
    country = st.selectbox('Choose your country:', np.sort(np.append(df.countries.unique(), "BE")), index=3)
with col2:
    carrier = st.selectbox('Choose your carrier:', df['carriers'].unique())

df_imp_x = query_imp_exp(df, carrier, country, 'imports').reset_index()
df_imp_x = df_imp_x[df_imp_x[country] != 0]
df_imp_x = df_imp_x.pivot_table(index="countries", values=country, columns="year")
df_imp_x.index.name = 'Annual import volume [TWh]'

df_exp_x = (-1 * query_imp_exp(df, carrier, country, 'exports')).reset_index()
df_exp_x = df_exp_x[df_exp_x[country] != 0]
df_exp_x = df_exp_x.pivot_table(index="countries", values=country, columns="year")
df_exp_x.index.name = 'Annual export volume [TWh]'

df_bar = pd.concat([df_imp_x.T, df_exp_x.T])

fig = px.bar(
    df_bar,
    title=f"Imports / Exports for {country} for {carrier} [TWh]",
    text_auto=".2s"
)

fig.update_traces(hovertemplate="%{y:,.0f}")
fig.update_yaxes(title_text='Annual exchange volume [TWh]', zeroline=True, zerolinewidth=3, zerolinecolor='black')
fig.update_xaxes(title_text='')
fig.update_layout(hovermode="closest",
                  legend_title_text='Exchange',
                  xaxis=dict(
                      tickmode='array',
                      tickvals=df_bar.index.unique().sort_values(),
                      ticktext=df_bar.index.unique().sort_values(),
                  )
                  )

st.plotly_chart(
    fig
    , use_container_width=True
)

buses = get_buses()
df_map = (
    pd.concat([df_imp_x.T, df_exp_x.T])
    .groupby("year").sum().T
    .join(buses)
    .reset_index()
    .rename(columns={'index': "country"})
    .melt(id_vars=["country", "lat", "lon"], value_name="Exchange [GW]", var_name="year")
)
df_map["Exchange abs [GW]"] = abs(df_map["Exchange [GW]"])
df_map["Balance"] = df_map["Exchange [GW]"].apply(lambda x: "Export" if x < 0 else "Import")
fig_map = px.scatter_mapbox(
    df_map,
    lat="lat",
    lon="lon",
    size="Exchange abs [GW]",
    color="Exchange [GW]",
    color_continuous_scale="RdYlBu",
    color_continuous_midpoint=0,
    mapbox_style="carto-positron",
    zoom=4,
    height=700,
    hover_name="country",
    animation_frame="year",
    title=f"Net balance of imports and exports for {country} for {carrier} [TWh]",
    hover_data={"Exchange [GW]": ":.2f", "Exchange abs [GW]": None},
)
fig_map.update_layout(sliders=[{"currentvalue": {"prefix": "Year: "}, "len": 0.8, "y": 0.07}])
fig_map.update_layout(updatemenus=[{"y": 0.07}])
fig_map.add_scattermapbox(
    lat=[buses.loc[country, "lat"]],
    lon=[buses.loc[country, "lon"]],
    name=country,
    hovertext=country,
    hovertemplate='',
    marker_color="lightgreen",
    marker_size=10,
)
st.plotly_chart(fig_map, use_container_width=True)

total_imp = (df_imp_x.sum())
try:
    df_imp_x.loc['Total'] = total_imp
except ValueError:
    pass
st.dataframe(df_imp_x
             .style
             .format(precision=2, thousands=",", decimal='.'),
             use_container_width=True)

total_exp = (df_exp_x.sum())
try:
    df_exp_x.loc['Total'] = total_exp
except ValueError:
    pass
st.dataframe(df_exp_x
             .style
             .format(precision=2, thousands=",", decimal='.'),
             use_container_width=True)
