..
  SPDX-FileCopyrightText: 2019-2023 The PyPSA-Eur Authors

  SPDX-License-Identifier: CC-BY-4.0

.. _veka_configurations:

##########################################
Model modifications
##########################################

To better match client needs, the model has been modified. To main action points have been implemented :
    - Addition of variable electric load in futur
    - Addtition of technology phase out capacity

Electric load
===========================
The historical electric load is disaggregated into historical share for :
	- Heat sector (where the historical electricity load is an intermediate vector for heat demand)
	- Industry sector (where the historical load is removed to be considered separately)
	- Residual electricity (all others sectors and loads)

The residual electricity is then dispatched over each country proportional to GDP and population.

As it appears that PyPSA-Eur uses the residual historical electricity load (which does not taken into account growth of appliances, which evolved between 27.5TWh in 2013 to 25.7TWh in 2022 up to  grow to 32.2TWh in 2030 )

Heat load
===========================

PyPSA-Eur considers JRC-IDEES historical load per country on an annual basis for hotwater and space heating purpose for residential and services subsectors.  SEE DOC FROM ENERGY DEMAND AND SUPPLY

Rules added
===========================

- :mod:`retrieve_load_futur`
- :mod:`build_country_profiles`
- :mod:`build_residual_load_profile`
- :mod:`build_future_load`
- :mod:`add_electricity_tomorrow`

External links
===========================

- Improve Gurobi using for `linopy` package (https://github.com/PyPSA/linopy/pull/162)
- Raised issue for `snakemake` package to better manage Gurobi licenses (https://github.com/snakemake/snakemake/issues/1801)
- Raised issue for `pulp` package to better manage Gurobi licenses (https://github.com/coin-or/pulp/issues/571)