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

Electric load considered
===========================
By default, PyPSA uses the historical electric load of a reference year for each country for future planning horizons. 

From this historical load, the historical electricity used :

* In the heat sector is removed from the total load to let PyPSA optimize the supply of historical heat demand ;
* In the industry sector is modified to fit future industry energy demand.

The residual share of electricity is otherwise considered constant over future planning horizons, which however does not take into account the evolution of appliances consumption.

The model has hence been modified to take into account the evolution of each sector's consumption. The annual sectorial electric demand projection for future horizons is extracted for each country from Climact's Pathways Explorer.

The variable electric load for future planning horizons is computed based on :

* A hourly demand per country for a reference year (from ENTSO-E) ;
* An annual demand per sector per country for a reference year (from PatEx) ;
* Hourly profiles per sector (representing which percentage if the sectorial annual lood is consumed at each hour) ;
* An annual demand per sector per country for future planning horizons (from PatEx).

Annual future demand
---------------------------

The Pathways Explorer can provided the electrical annual demand per country according to country-specific scenarios for future years, with a granularity up to subsetcors.

    - Ajouter description PatEx + Lien vers le site + Slide https://climact.sharepoint.com/:p:/r/teams/POWER/Documents%20partages/General/PyPSA-Eur%20for%20starters.pptx?d=w2e4300e324874394b08a561826989c9b&csf=1&web=1&e=mLpADP&nav=eyJzSWQiOjQzNzgsImNJZCI6NDA1NzMwMzY0MH0 # ToDo


Profiles definition
---------------------------

    - Graphes ? # ToDo VLA

The sectors considered by the PatEx annual electricity demand are :

* Industry
* Heat 
* Transport
* Power Supply
* Residual load (others)
	
Each of those sectors are modeled with the exception of the residual load which, by definition, is defined as what is left after substracting to the total load the different sectors, meaning no particular profile is defined for it. 

* Industry 	: Industry electrical demand profile is considered as flat over the whole year
* Heat 		: Heat electrical demand profile is calculated similarly to PyPSA-methodology for space heating and hotwater demand :

  * An intraday hourly profile, depending on the sector (residential/service), the heat type (hotwater/space heating)and on week days/week-ends
  * An annual daily profile, considered flat for hotwater and spread accross the year according the daily average Heating Degree Day considering a threshold temeprature of 15°C
* Transport	: Transport electrical demand profiles are based on hourly profiles available at a week scale provided by the German Federal Highway Research Institute (BASt). Profiles for different types of vehicles are available ; the profile of all land transport types vehicles combined is considered as a proxy for electric rail, as no profile is available.
* Power supply : Power supply electrical demand profile (i.e. losses) is considered to be proportional to the total load for each time frame. Losses are assumed to be equal to represent 5% of the total load.


Methodology 
---------------------------
The future hourly electric demand for future planning horizons it build by :

* Building a reference hourly demand for each sector, based on defined sector profiles and on annual demand ;
* Building a reference hourly profile for the Residual load sector based on the reference hourly demand and on the hourly sectorial demand ;
* Building future hourly demands based on defined and built sector profiles and on future annual demands
	

As it appears that PyPSA-Eur uses the residual historical electricity load growth of appliances is not taken into account even though it evolved over time

* 27.5TWh in 2013
* 25.7TWh in 2022
* 32.2TWh in 2030
In Belgium using according to Climact's previsions for Elia 

- Lien vers l'étude ELIA # ToDo VLA

Technology phase out 
===========================
Some scenarios might want to explore what a future energy system would look like considering teh phase out of somme due to political choices (e.g. ban on coal powerplants by 2030).

A new option has hence been added to allow to phase out assets of specified conventional technologies : 

* Existing assets are adapted so that they are removed starting from the phase out date
* The lifetime of new assets is adapted to make sure they are removed at their phase out date. In case the lifetime of the new assets are reduced, annualized investment costs for new assets of technologies to phase out are adapted accordingly ; this is reflected through a higher annuity in the annualized capital cost calculation.





PyPSA-Eur considers JRC-IDEES historical load per country on an annual basis for hotwater and space heating purpose for residential and services subsectors.  SEE DOC FROM ENERGY DEMAND AND SUPPLY


Rules added
===========================

    - Update with https://climact.sharepoint.com/:p:/r/sites/PROSPECTIVE/Documents%20partages/VEKA/2022%20ENERGIE%20STUDIE%202050%20+/2022%20Energiestudie%202050%20shared/4.%20Workstream/5.%20PyPSA%20modelling/PyPSA%20modelling.pptx?d=w77e33e38d3714c63bbeefb693bc76110&csf=1&web=1&e=uxTqBd&nav=eyJzSWQiOjIxNDc0NzI0MzMsImNJZCI6MjQ5NjA3ODIxM30 # ToDo

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