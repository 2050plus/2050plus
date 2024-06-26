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
st.markdown("The energy imports and exports between countries in the system, for all carriers, countries and years.")


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


def query_imp_exp(df, carriers, country, year, imports_exports):
    df_imp_exp = (
        df.query(""
                 "carriers == @carriers & "
                 "year == @year & "
                 "imports_exports == @imports_exports"
                 )
        .drop(["carriers", "year", "imports_exports"], axis=1)
        .set_index('countries')
        [country]
    )
    return df_imp_exp


col1, col2, col3 = st.columns(3)
with col1:
    country = st.selectbox('Choose your country:', df["countries"].unique(), index=12)
with col2:
    carrier = st.selectbox('Choose your carrier:', df['carriers'].unique())
with col3:
    year = st.selectbox('Choose your year:', df["year"].unique())
df_imp_exp = (
    pd.concat([query_imp_exp(df, carrier, country, year, 'imports'),
               -1 * query_imp_exp(df, carrier, country, year, 'exports')],
              axis=1, keys=['imports', 'exports'])
)
df_imp_exp.rename(mapper=lambda x: x.capitalize(), axis=1, inplace=True)

fig = px.bar(
    df_imp_exp,
    title=f"Imports / Exports for {country} for {carrier} [TWh]",
    text_auto=".2s"
)

fig.update_traces(hovertemplate="%{y:,.0f}")
fig.update_yaxes(title_text='Annual exchange volume [TWh]')
fig.update_xaxes(title_text='Countries')
fig.update_layout(hovermode="x unified",
                  legend_title_text='Exchange')

st.plotly_chart(
    fig
    , use_container_width=True
)

df_imp_exp_ = df_imp_exp.drop(df_imp_exp.query('Imports == 0 and Exports ==0').index)
df_imp_exp_.rename(mapper=lambda x: x + " [TWh]", axis=1, inplace=True)

st.subheader(f"Annual {carrier} exchange volumes of {country} for {year} ")

st.dataframe(df_imp_exp_
             .style
             .format(precision=2, thousands=",", decimal='.'),
             use_container_width=True)
