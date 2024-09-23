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
# TODO : check if scenario is relevant
scenario, compare = st_side_bar()

st.title("Balancing capacities")
st.markdown(
    "The balancing capacities installed per country, technologies and year. Balancing units can help balancing the network, either by increasing production or by shifting consumption.")
with st.expander("**Why balancing the network is important ?**"):
    st.write(
        "Electricity as such cannot be stored. For the electricity grid to function, electricity consumption must always be equal to electricity production. As soon as an individual consumes 1 kilowatt, 1 kilowatt must simultaneously be produced by a generating unit on the grid. Thermal power stations (gas, coal, etc.), dams and nuclear power stations are called controllable: their output can be adjusted according to demand. They can be load following, i.e. they can adapt to fluctuations in demand. On the other hand, wind turbines and solar panels are called intermittent: their production depends on weather conditions and day/night cycles.")
    st.write(
        "Without controllable means, production would not be able to adapt to consumption. When consumption exceeds production, the frequency of the electricity grid (50 Hz) drops. This is the result of the slowing down of all the generators in the electricity network, which are under pressure to meet demand. Various measures, such as load shedding, are planned to reduce some of the consumption. If this is not sufficient, the frequency will continue to decrease until the safety shutdown of production units and the collapse of the electrical grid (example: 2003 blackout in Italy).")
    st.write(
        "In a 100% renewable energy system, there are no longer controllable units such as gas, coal and nuclear power stations that can adapt to fluctuations in demand. Other balancing capacity is therefore required. These can either increase production (e.g. hydro), shift consumption to a more favourable time when production is higher (e.g. heat pumps, BEV chargers), or convert electricity to store it during periods of surplus for redistribution during periods of shortage (e.g. PHS, battery chargers).")


@st.cache_data(show_spinner="Retrieving data ...")
def get_data(scenario, mode):
    return (
        pd.read_csv(
            Path(network_path, scenario_dict[scenario]["path"],
                 "balancing_" + mode + ".csv"),
            header=0,
        )
    )


# %%
data = get_data(scenario, "capacities")
if compare != '-':
    df_compare = get_data(compare, "capacities")
    idx = ["country", "carrier"]
    data = (
        (data.set_index(idx) - df_compare.set_index(idx))
        .reset_index()
    )
df = data.copy()
consumption = {'BEV charger', 'air heat pump', 'battery charger', 'ground heat pump', 'home battery charger',
               'water tanks charger'}

condition = df['carrier'].isin(consumption)
st.write()
df.loc[condition, df.columns[2:]] *= -1

st.header("Installed capacities per country")
st.markdown(
    'Negative values correspond to assets helping to balance the grid by consuming energy, positive value to assets helping to balance the grid by producing energy')

all = ['ENTSO-E area']
country = st.selectbox('Choose your country:', all + list(df.country.unique()))
if not ('ENTSO-E area' in country):
    df = df.query("country in @country")
else:
    df = df.query("country != 'FL'")

df = (
    df.drop(columns=['country'])
    .groupby(['carrier'])
    .sum(numeric_only=True)
    .rename_axis(index={"carrier": "Technologies"})
)

fig = px.bar(
    df,
    width=1000,
    height=500,
    title="Balancing installed capacities [GW]",
    barmode="group",
    text_auto=".2s"
)

fig.update_yaxes(title_text='Installed power [GW]')
fig.update_xaxes(title_text='Technologies')
fig.update_traces(hovertemplate="%{y:,.0f}")
fig.update_layout(hovermode="x unified",
                  legend_title_text='Technologies')

st.plotly_chart(fig, use_container_width=True)

# %%
st.header("Annual energy output per technology")

data2 = get_data(scenario, "supply")
if compare != '-':
    df_compare = get_data(compare, "supply")
    idx = ["country", "carrier"]
    data2 = (
        (data2.set_index(idx) - df_compare.set_index(idx))
        .reset_index()
    )
df2 = data2.copy()

technology = st.selectbox('Choose your technology:', list(df2.carrier.unique()))

df2 = (df2
       .query("carrier == @technology")
       .drop(columns=['carrier'])
       .set_index('country')
        .fillna(0)
       )
df2.index.name = "Energy output [GWh]"

df_map = (
    df2
    .join(get_buses())
    .reset_index()
    .rename(columns={"Energy output [GWh]": "country"})
    .melt(id_vars=["country", "lat", "lon"], value_name="Energy output [GWh]", var_name="year")
    .assign(absolute_size=lambda x: x["Energy output [GWh]"].abs())
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
    title="Energy output [GWh]",
    hover_data={"Energy output [GWh]": ":.2f"}
)
fig_map.update_layout(sliders=[{"currentvalue": {"prefix": "Year: "}, "len": 0.8, "y": 0.07}])
fig_map.update_layout(updatemenus=[{"y": 0.07}])
st.plotly_chart(fig_map, use_container_width=True)

fig_bar = px.bar(
    df2,
    height=500,
    barmode="group"
)

fig.update_yaxes(title_text='Energy [GWh]')
fig_bar.update_xaxes(title_text='Countries')
fig_bar.update_traces(hovertemplate="%{y:,.1f}")
fig_bar.update_layout(hovermode="x unified",
                      legend_title_text='Technologies')

st.plotly_chart(fig_bar, use_container_width=True)

st.dataframe(df2.style.format(precision=2, thousands=",", decimal='.'), use_container_width=True)

# # %%
# st.header("Focus on EV batteries (BEVs)")

# df3 = (data2.copy()
#        .query("carrier == 'BEV charger' | carrier == 'V2G'")
#        .drop(columns=['carrier'])
#        .set_index(['country'])
#        )

# # divides V2G charger by BEV
# df3 = df3 / df3.shift()
# df3 = df3[1::2]

# fig = px.bar(
#     df3,
#     width=1000,
#     height=600,
#     title="Fraction of the energy that the BEVs return to the grid",
#     barmode="group",
#     text_auto=".2s"
# )

# fig.update_layout(
#     xaxis_title="Countries",
#     hovermode="x unified",
#     legend_title_text="Technologies"
# )
# fig.update_traces(hovertemplate="%{y:,.2f}")

# st.plotly_chart(fig)

# # %%
# st.header("Focus on EV batteries (BEVs) - details for a country")

# # st.write("- **BEV charger** is the electrical consumption of the EV batteries from the grid")
# # st.write("- **V2G** is the electrical production from the EV batteries to the grid")


# df4 = data2.copy()

# country = st.selectbox('Choose your country:', list(df4.country.unique()), key="unique_key_by_widget")
# df4 = df4.query("country in @country")

# df4 = (df4
#        .query("carrier == 'BEV charger' | carrier == 'V2G'")
#        .drop(columns=['country'])
#        .set_index(['carrier'])
#        )

# df4 = df4.transpose()

# fig = px.bar(
#     df4,
#     title="Annual energy consumption/production [GWh]",
#     barmode="group",
#     text_auto=".2s"
# )

# fig.update_yaxes(title_text='Energy [GWh]')
# fig.update_xaxes(title_text='Years')
# fig.update_traces(hovertemplate="%{y:,.0f}")
# fig.update_layout(hovermode="x unified",
#                   legend_title_text='Technologies')

# st.plotly_chart(fig)
