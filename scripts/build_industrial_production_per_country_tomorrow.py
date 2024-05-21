# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2020-2024 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT
"""
Build future industrial production per country.
"""

import pandas as pd
from _helpers import set_scenario_config
from prepare_sector_network import get

if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake("build_industrial_production_per_country_tomorrow")
    set_scenario_config(snakemake)

    params = snakemake.params.industry

    investment_year = int(snakemake.wildcards.planning_horizons)

    fn = snakemake.input.industrial_production_per_country
    production = pd.read_csv(fn, index_col=0)

    keys = ["Integrated steelworks", "Electric arc"]
    total_steel = production[keys].sum(axis=1)

    st_primary_fraction = get(params["St_primary_fraction"], investment_year)
    dri_h2_fraction = get(params["DRI_H2_fraction"], investment_year)
    dri_ch4_fraction = get(params["DRI_CH4_fraction"], investment_year)
    
    dri_fraction = dri_h2_fraction + dri_ch4_fraction
    
    int_steel = production["Integrated steelworks"].sum()
    fraction_persistent_primary = st_primary_fraction * total_steel.sum() / int_steel

    dri_h2 = (
        dri_h2_fraction * fraction_persistent_primary * production["Integrated steelworks"]
    )
    dri_ch4 = (
        dri_ch4_fraction * fraction_persistent_primary * production["Integrated steelworks"]
    )
    production.insert(2, "DRI H2 + Electric arc", dri_h2)
    production.insert(2, "DRI CH4 + Electric arc", dri_ch4)

    not_dri = 1 - dri_fraction
    production["Integrated steelworks"] = (
        not_dri * fraction_persistent_primary * production["Integrated steelworks"]
    )
    production["Electric arc"] = (
        total_steel
        - production["DRI H2 + Electric arc"]
        - production["DRI CH4 + Electric arc"]
        - production["Integrated steelworks"]
    )

    keys = ["Aluminium - primary production", "Aluminium - secondary production"]
    total_aluminium = production[keys].sum(axis=1)

    key_pri = "Aluminium - primary production"
    key_sec = "Aluminium - secondary production"

    al_primary_fraction = get(params["Al_primary_fraction"], investment_year)
    fraction_persistent_primary = (
        al_primary_fraction * total_aluminium.sum() / production[key_pri].sum()
    )

    production[key_pri] = fraction_persistent_primary * production[key_pri]
    production[key_sec] = total_aluminium - production[key_pri]

    production["HVC (mechanical recycling)"] = (
        get(params["HVC_mechanical_recycling_fraction"], investment_year)
        * production["HVC"]
    )
    production["HVC (chemical recycling)"] = (
        get(params["HVC_chemical_recycling_fraction"], investment_year)
        * production["HVC"]
    )

    production["HVC"] *= get(params["HVC_primary_fraction"], investment_year)

    fn = snakemake.output.industrial_production_per_country_tomorrow
    production.to_csv(fn, float_format="%.2f")
