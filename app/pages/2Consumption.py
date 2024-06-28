from pathlib import Path

import pandas as pd
import plotly.express as px
import streamlit as st
from st_common import PROFILES_AREA
from st_common import YEARS
from st_common import network_path
from st_common import scenario_dict
from st_common import st_page_config
from st_common import st_side_bar

st_page_config(layout="wide")
scenario = st_side_bar()

st.title("Consumption per carrier")

st.markdown("The total energy consumption per year, country and subsector. This data is currently shown at system level (ENTSO-E area), Belgium (BE) and Flanders (FL) due to the very large quantity of data that needs to be handled for every country in the system.")


@st.cache_data(show_spinner="Retrieving data ...")
def get_data(scenario, year, selected_area):
    area = selected_area if selected_area != "ENTSO-E area" else ''
    return (
        pd.read_csv(
            Path(network_path, scenario_dict[scenario]["path"],
                 f"load_temporal_{area.lower()}_{year}.csv".replace('__', '_')),
            header=[1, 2]
        )
        .set_index(("carrier", "sector"))
    )


col1, col2 = st.columns(2)
with col1:
    selected_area = st.selectbox('Choose area :', PROFILES_AREA)

dfx = []
for y in YEARS:
    dfi = get_data(scenario, y, selected_area)
    dfi = ((dfi.sum() / 1e3
            * 8760 / len(dfi.axes[0]))
           .to_frame(name=y))
    dfx.append(dfi)
dfx = pd.concat(dfx, axis=1).fillna(0).sort_values(by=y, ascending=False)
dfx.index.name = 'Annual consumption [TWh]'

with col2:
    carrier = st.selectbox('Choose your carrier:', dfx.index.unique(0).sort_values(), index=6)
dfx = dfx.loc[carrier]


st.subheader(f"Annual {carrier} consumption per sector")

total = (dfx.sum())
dfx.loc['Total'] = total
dfx.index.name = "Annual consumption [TWh]"

fig = px.bar(
    dfx,
    title=f"Consumption in {selected_area} for {carrier} [TWh]",
    barmode="group",
    text_auto=".2s"
)

fig.update_traces(hovertemplate="%{y:,.0f}")
fig.update_layout(hovermode="x unified")
fig.update_yaxes(title_text='Consumption [TWh]')
fig.update_xaxes(title_text='Sectors')
fig.update_layout(legend_title_text='Years')

st.plotly_chart(
    fig
    , use_container_width=True
)

st.dataframe(
    dfx
    .style
    .format(precision=2, thousands=",", decimal='.'),
    use_container_width=True
)

st.divider()


st.subheader(f"Consumption profiles per carrier")

st.markdown(
    "The load 3-hourly profiles for every carrier, year and subsector. You can zoom on these interactive graphs for specific time windows and you can also select/deselect various categories if you want.")

year = st.selectbox('Choose the year:', YEARS)
data = get_data(scenario, year, selected_area)

df = data[carrier]

fig = px.area(
    df,
    title=f"System consumption profile for {carrier} [GW]",
)
fig.update_traces(hovertemplate="%{y:,.0f}",
                  line=dict(width=0.1))
fig.update_layout(legend_traceorder="reversed",
                  hovermode="x unified",
                  legend_title_text='Technologies')
fig.update_yaxes(title_text='Consumption [GW]')
fig.update_xaxes(title_text='Timesteps')

st.plotly_chart(
    fig
    , use_container_width=True
)
