from pathlib import Path

import pandas as pd
import plotly.express as px
import streamlit as st
from st_common import YEARS
from st_common import network_path
from st_common import scenario_dict
from st_common import st_page_config
from st_common import st_side_bar

AREAS = ["EU27+TYNDP", "BE"]

st_page_config(layout="wide")
scenario = st_side_bar()

st.title("Consumption profiles per carrier")
st.markdown("The load 3-hourly profiles for every carrier, year and subsector. This data is currently shown at system level (EU27+TYNDP) and Belgium (BE) due to the very large quantity of data that needs to be handled for every country in the system. You can zoom on these interactive graphs for specific time windows and you can also select/deselect various categories if you want.")

selected_area = st.selectbox('Choose area :', AREAS)

@st.cache_data(show_spinner="Retrieving data ...")
def get_data(scenario, year, selected_area):
    area = selected_area if selected_area != "EU27+TYNDP" else ''
    return (
        pd.read_csv(
            Path(network_path, scenario_dict[scenario]["path"],
                 f"load_temporal_{area.lower()}_{year}.csv".replace('__', '_')),
            header=[1, 2]
        )
        .set_index(("carrier", "sector"))
    )


# %%Cell Name
# should be able to 
# - Display per carrier
# - 3h load profile
# - eventually per country
# - eventuelly per subtype of supply

col1, col2 = st.columns(2)
with col1:
    year = st.selectbox('Choose the year:', YEARS)
data = get_data(scenario, year, selected_area)

with col2:
    carrier = st.selectbox('Choose your carrier:', data.columns.get_level_values(0).unique(), index=4)
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

st.subheader(f"Annual {carrier} consumption per sector for {year} ")

st.table(
    (df.sum() / 1e3
     * 8760 / len(df.axes[0]))
    .rename('Annual consumption [TWh]')
    .sort_values(ascending=False)
    .to_frame()
    .style
    .format(precision=2, thousands=",", decimal='.')
)
