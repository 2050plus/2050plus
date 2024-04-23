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

st.title("Production profiles per carrier")
st.markdown("The production 3-hourly profiles for every carrier, year and subsector. This data is currently shown at system level (EU27+TYNDP) due to the very large quantity of data that needs to be handled for every country in the system. You can zoom on these interactive graphs for specific time windows and you can also select/deselect various categories if you want.")


@st.cache_data(show_spinner="Retrieving data ...")
def get_data(scenario, year):
    return (
        pd.read_csv(
            Path(network_path, scenario_dict[scenario]["path"], "graph_extraction_st",
                 "supply_temporal_" + year + ".csv"),
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
years = ['2030', '2040', '2050']
col1, col2 = st.columns(2)
with col1:
    year = st.selectbox('Choose the year:', years)
data = get_data(scenario, year)

with col2:
    carrier = st.selectbox('Choose your carrier:', data.columns.get_level_values(0).unique(), index=4)
df = data[carrier]

fig = px.area(
    df,
    title=f"System production profile for {carrier} [GW]",
)
fig.update_traces(hovertemplate="%{y:,.0f}",
                  line=dict(width=0.1))
fig.update_layout(legend_traceorder="reversed",
                  hovermode="x unified",
                  legend_title_text='Technologies')
fig.update_yaxes(title_text='Production [GW]')
fig.update_xaxes(title_text='Timesteps')

st.plotly_chart(
    fig
    , use_container_width=True
)

st.subheader(f"Annual {carrier} production per technology for {year} ")

df_table = (
    (df.sum() / 1e3
     * 8760 / len(df.axes[0]))
    .rename('Annual production [TWh]')
    .sort_values(ascending=False)
    .to_frame()
    .style
    .format(precision=2, thousands=",", decimal='.')
)

st.table(df_table)
