..
  SPDX-FileCopyrightText: 2019-2023 The PyPSA-Eur Authors

  SPDX-License-Identifier: CC-BY-4.0

.. _veka_configurations:

##########################################
Model modifications
##########################################


To better match client needs, the model has been modified. To main action points have been implemented :

* Addition of variable electric load in futur
* Addition of technology phase out capacity

Electric load
===========================
By default, PyPSA uses the historical electric load of a reference year for each country for future planning horizons. 

From this historical load, the historical electricity used :

* In the heat sector is removed from the total load to let PyPSA optimize the supply of historical heat demand ;
* In the industry sector is modified to fit future industry energy demand.

The residual share of electricity is otherwise considered constant over future planning horizons, which however does not take into account the evolution of appliances consumption.

The model has hence been modified to allow the model to take into account an evolution of each sector's consumption. The annual sectorial electric demand projection for future horizons is extracted for each country from Climact's Pathways Explorer.

The variable electric load for future planning horizons is computed based on :

* A hourly demand per country for a reference year (from ENTSO-E) ;
* An annual demand per sector per country for a reference year (from PatEx) ;
* Hourly profiles per sector (representing which percentage if the sectorial annual lood is consumed at each hour) ;
* An annual demand per sector per country for future planning horizons (from PatEx).
	
The sectors considered by the PatEx annual electricity demand are :

* Industry
* Heat 
* Transport
* Power Supply
* Residual load (others)
	
By definition, the residual load is defined as what is left after substracting to the total load the different sectors, meaning no particular profile is defined for it. 

To build the variable electric load for future planning horizons it is necessary :

* To build a reference hourly demand for each sector, based on defined sector profiles and on annual demand ;
* To build a reference hourly profile for the Residual load sector based on the reference hourly demand and on the hourly sectorial demand ;
* To build future hourly demands based on defined and built sector profiles and on future annual demands
	

	
	
	

As it appears that PyPSA-Eur uses the residual historical electricity load growth of appliances is not taken into account even though it evolved over time

* 27.5TWh in 2013
* 25.7TWh in 2022
* 32.2TWh in 2030
In Belgium using according to Climact's previsions for Elia 

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

- Improve Gurobi usage for `linopy` package (https://github.com/PyPSA/linopy/pull/162)
- Raised issue for `snakemake` package to better manage Gurobi licenses (https://github.com/snakemake/snakemake/issues/1801)
- Raised issue for `pulp` package to better manage Gurobi licenses (https://github.com/coin-or/pulp/issues/571)