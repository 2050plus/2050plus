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

st.title("Renewable production potentials")
st.markdown("The total potential for RES power production capacity considered by the model for the various technologies and countries.")


@st.cache_data(show_spinner="Retrieving data ...")
def get_df(scenario):
    return (
        pd.read_csv(
            Path(network_path, scenario_dict[scenario]["path"], "graph_extraction_st", "res_potentials.csv"),
            header=0,
        )
    )


# %%
data = get_df(scenario)
df = data.copy()

st.header("Potentials per carrier")

df = df.groupby("country").sum()
carrier = st.multiselect('Choose your carrier:', list(df.columns.unique()), default=list(df.columns.unique())[2:4])
df = df.loc[:, carrier]
carrier_list = ' & '.join(list(map(str.capitalize, carrier)))

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

st.divider()

st.header("Potentials per country")
df_tab = data.copy()
country = st.selectbox('Choose your country:', ['EU27 + TYNDP'] + list(df_tab.country.unique()))
if country != 'EU27 + TYNDP':
    df_tab = df_tab.set_index('country').loc[country]
else:
    df_tab = df_tab.sum(numeric_only=True)

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

st.table(df_tab)
