from pathlib import Path

import pandas as pd
import plotly.express as px
import streamlit as st
from st_common import COSTS_AREA
from st_common import network_path
from st_common import scenario_dict
from st_common import st_page_config
from st_common import st_side_bar

AREAS = ["ENTSO-E area", "EU27", "BE"]

st_page_config(layout="wide")
scenario = st_side_bar()

st.title("Costs")


@st.cache_data(show_spinner="Retrieving data ...")
def get_data(scenario, path):
    df = (
        pd.read_csv(
            Path(network_path, scenario_dict[scenario]["path"], path),
            header=0
        )
    )
    return df


# %% Cost segment
st.header("Cost by unit segment")

df_cost_segments = get_data(scenario, "costs_segments.csv").set_index("config")
col1, col2 = st.columns([4, 4])
with col1:
    selected_cost_segment = st.selectbox("Choose your segment :",
                                         ["Total"] + list(df_cost_segments.cost_segment.unique()))
    if selected_cost_segment != "Total":
        df_cost_segments = df_cost_segments.query("cost_segment in @selected_cost_segment")
    else:
        df_cost_segments = df_cost_segments.query("cost_segment != 'Net_Imports'")
    st.text("A negative value means that the area is exporting and thus making a profit")
with col2:
    selected_area = st.selectbox("Choose area :", COSTS_AREA)
    df_cost_segments = df_cost_segments[df_cost_segments.index.str.endswith(COSTS_AREA[selected_area])]
    st.text("tot includes all modeled countries, so imports and exports = 0")

df_cost_segments = df_cost_segments.groupby(by="cost/carrier").sum().drop(columns=["cost_segment"])

df_cost_segments = df_cost_segments.div(1e9)

df_cost_segments.loc["Total"] = df_cost_segments.sum()

fig = px.bar(
    df_cost_segments,
    width=1000,
    height=400,
    title=f"{selected_cost_segment} costs for {selected_area} [Billion € / y]",
    barmode="group",
    text_auto=".3s"
)

fig.update_yaxes(title_text="Cost [Billion € / y]")
fig.update_xaxes(title_text="Cost type")
fig.update_traces(hovertemplate="%{y:,.3s}")
fig.update_layout(hovermode="x unified",
                  legend_title_text="Year")

st.plotly_chart(
    fig,
    use_container_width=True
)
df_cost_segments.index.set_names("Costs per type/carrier [B€/year]",inplace=True)
st.dataframe(df_cost_segments.style.format(precision=2, thousands=",", decimal='.'))

# %% Cost year
st.header("Cost by year")
df_cost_years = get_data(scenario, "costs_years.csv")

df_cost_years[["year", "area"]] = df_cost_years["config"].str.split('_', expand=True)
df_cost_years = df_cost_years.drop(columns=["config"])
df_cost_years = df_cost_years.set_index("cost_segment")

col1, col2 = st.columns([4, 4])
with col1:
    selected_year = st.selectbox("Choose year :", list(df_cost_years.year.unique()))
    df_cost_years = df_cost_years.query("year in @selected_year")
with col2:
    selected_area = st.selectbox("Choose area  :", list(COSTS_AREA))
    selected_area_ = COSTS_AREA[selected_area]
    df_cost_years = df_cost_years.query("area in @selected_area_")

df_cost_years = df_cost_years.drop(columns=["year", "area"])

df_cost_years = df_cost_years.div(1e9)

fig = px.bar(
    df_cost_years,
    width=1000,
    height=400,
    title=f"{selected_year} costs for {selected_area} [Billion € / y]",
    text_auto=".3s"
)

fig.update_yaxes(title_text="Cost [Billion € / y]")
fig.update_xaxes(title_text="Segment")
fig.update_traces(hovertemplate="%{y:,.3s}")
fig.update_layout(hovermode="x unified",
                  legend_title_text="Cost type")

st.plotly_chart(
    fig,
    use_container_width=True
)

df_cost_years.index.set_names("Costs per unit segment/type [B€/year]",inplace=True)
st.dataframe(df_cost_years.assign(total=df_cost_years.sum(axis=1)).style.format(precision=2, thousands=",", decimal='.'))


# %%
st.header("Marginal price of methane, electricity and hydrogen")
st.markdown(
    "The marginal price represents the cost of an additional unit of an energy carrier at one country at a given time. Please note that the shown value does not take into account taxes nor trading dynamic.")

st.subheader("Compare two areas over years")

df = get_data(scenario, "marginal_prices.csv")

col1, col2 = st.columns([4, 4])
with col1:
    country = st.selectbox("Choose your area:", list(df.countries.unique()))

    if not ("ENTSO-E area" in country):
        df1 = df.query("countries in @country").drop(columns="countries").set_index("carrier")
    else:
        raise Exception ("Fix me to remove Flanders")

    df1 = df1.rename(columns=lambda x: x + " Annual average " if not (x.endswith("_std")) else x.replace("_std",
                                                                                                         " Standard Deviation")).T
    df1 = df1.rename(columns=lambda x: x + " [€/MWh]")

    st.table(df1
             .style
             .format(precision=2, thousands=",", decimal='.'))
with col2:
    country2 = st.selectbox("Choose your country:", list(df.countries.unique()) + [''])

    if not ("ENTSO-E area" in country2):
        df2 = df.query("countries in @country2").drop(columns="countries").set_index("carrier")
    else:
        raise Exception ("Fix me to remove Flanders")

    df2 = df2.rename(columns=lambda x: x + " Annual average " if not (x.endswith("_std")) else x.replace("_std",
                                                                                                         " Standard Deviation")).T
    df2 = df2.rename(columns=lambda x: x + " [€/MWh]")

    st.table(df2
             .style
             .format(precision=2, thousands=",", decimal='.'))

# %% Box plot for costs

st.subheader("Compare variability through Europe")

df_t = get_data(scenario, "marginal_prices_t.csv")

col1, col2, col3 = st.columns([0.2, 0.2, .5])
with col1:
    carrier = st.selectbox("Choose areas to compare:", list(df_t.carrier.unique()), index=1)
with col2:
    year = st.selectbox("Choose year to select:", list(df_t.year.unique()))
with col3:
    countries_to_display = st.multiselect("Choose your carrier:", list(df.countries.unique()),
                                          default=["GB", "FL", "FR", "LU", "DE", "NL"])

df_t = (
    df_t.query("carrier == @carrier and year==@year and countries in @countries_to_display")
    .drop(columns=["carrier", "year"])
    .set_index(["countries"])
    .rename_axis(columns=["marginal price"])
    .T
)

fig = px.box(df_t)
fig.update_traces(hovertemplate="%{y:,.0f}",
                  line=dict(width=1))
fig.update_yaxes(title_text=f"Marginal price for {carrier} [€/MWh]")
fig.update_xaxes(title_text="Areas")

st.plotly_chart(fig)
