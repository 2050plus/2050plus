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

st.title("Imports and exports per carrier")
st.markdown("The energy imports and exports between countries in the system, for all carriers and countries. A negative value means that the area is exporting.")


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


def query_imp_exp(df, carriers, country, imports_exports, year=None):
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
    return df_imp_exp


col1, col2 = st.columns(2)
with col1:
    country = st.selectbox('Choose your country:', df["countries"].unique(), index=12)
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

fig = px.bar(
    pd.concat([df_imp_x.T, df_exp_x.T]),
    title=f"Imports / Exports for {country} for {carrier} [TWh]",
    text_auto=".2s"
)

fig.update_traces(hovertemplate="%{y:,.0f}")
fig.update_yaxes(title_text='Annual exchange volume [TWh]', zeroline=True, zerolinewidth=3, zerolinecolor='black')
fig.update_xaxes(title_text='')
fig.update_layout(hovermode="closest",
                  legend_title_text='Exchange')

st.plotly_chart(
    fig
    , use_container_width=True
)

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
