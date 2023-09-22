..
  SPDX-FileCopyrightText: 2019-2023 The PyPSA-Eur Authors

  SPDX-License-Identifier: CC-BY-4.0

.. _veka_configurations:

##########################################
Configurations
##########################################

PyPSA-Eur is able to provide the energy supply of an energy system given :

* The configuration of the system by year-varying parameters (such as carbon budget or primary route share in steel production) and fixed parameters (such as maximum potential per renewable technologies or charging power of EVs) ;
* The list of technologies used ;
* Techno-economic parameters (such as investment costs, efficiency, FOM, VOM, lifetime, discount rate, etc.).

The quality of the optimization results depends on the database quality, as :

* The evolution trajectory of year-dependent parameters has an impact on the energy demand for some parameters (by dispatching the total demand over different vectors). Fixed parameters have an impact on some technologies (through potential, minimum capacity factor, etc.) ;
* The addition of a technology might lead to an energy system significantly cheaper throughout the optimization, due to otherwise non-existing or uninteresting interactions between technologies ;
* Some technologies might not be considered in the cost-optimal system because of too high CAPEX/OPEX, or might be massively installed because of too optimistic costs.

The way PyPSA deals with those different topics is explained in the following sections.

Technological assumptions
===========================

PyPSA-Eur optimization of the energy system is done by computing the cost-optimal sizing of each technology per geographical location.

A technology can be used for

* Generation of energy using energy carrier(s) to produce other energy carrier(s) :

  * i.e. Using coal to produce electricity, CO2 and captured CO2 for a coal powerplant with CC ;
  * i.e. Using hydrogen and captured CO2 to produce synthetic oil for the Fischer-Tropsh process;
* Storage of energy under a specific energy vector:

  * i.e. Centralized Thermal Energy Storage;
  * i.e. Hydrogen underground storage;
* Transmission of energy vector between two geographical locations:

  * i.e. AC and DC lines;
  * i.e. Methane pipelines;
  * i.e. CO2 pipelines;
  
Some technologies are added to the system only if an energy sector is considered in the optimization. An exhaustive list is given here below, sorted by module and with each energy carriers the technology uses.

Base technologies
---------------------------

.. csv-table::
   :header-rows: 1
   :file: configtables/base_techs.csv

Heat technologies
---------------------------

.. csv-table::
   :header-rows: 1
   :file: configtables/heat_techs.csv

Biomass technologies
---------------------------

.. csv-table::
   :header-rows: 1
   :file: configtables/biomass_techs.csv
   
Biomass and Heat technologies
---------------------------

.. csv-table::
   :header-rows: 1
   :file: configtables/biomass_heat_techs.csv
   
Industry technologies
---------------------------

.. csv-table::
   :header-rows: 1
   :file: configtables/ind_techs.csv
   
Transport technologies
---------------------------

.. csv-table::
   :header-rows: 1
   :file: configtables/transport_techs.csv


Techno-economic parameters
===========================

The definition of the technologies in PyPSA is done by retrieving data from a cost database and formatting it into the metrics used by PyPSA-Eur, namely :

* Annualized Capital cost 	(€/MW/year)
* Marginal cost 			(EUR/MWh)
* Lifetime 					(years)
* Efficiency(ies)			(MWhout/MWhin)
* CO2 intensity   			(tCO2/MWhout)
* Potential 				(MWhmax)
* Carrier(s)

The cost database (https://github.com/pypsa/technology-data) has a granularity of up to 5 years and is mostly based on the Danish Energy Agency (DEA) forecasts (March 2018 - August 2023).

It must be noted nonetheless that for some technologies, some techno-economic parameters are set from the configuration file instead of the cost database.

Configuration file
===========================

PyPSA-Eur optimization is mostly based on the choice of the technologies used and the techno-economic parameters from the cost database.

Some additional parameters can nonetheless be set from a separate configuration file. Those parameters can be grouped under different categories :

* On/off technology use : Levers (de)activating some technologies in PyPSA optimization

  * i.e. Conventional technologies to consider in future planning horizons;
  * i.e. Use of micro-CHP, solid biomass to liquid, etc.;
  * i.e. Considering distribution electric and/or gas networks;

* Technology parameters : techno-economic parameters that were not set from the cost database or that alter technologies

  * i.e. Potentials and correction factors for renewables;
  * i.e. Heat pump sink temperature;

* Demand-related parameters: share between different energy carriers of a given demand. They can be fixed over the explored time horizons or year-dependent

  * i.e. Share of primary route in steel production;
  * i.e. Share of EV/ICE/FC vehicles for land transport compared to today's demand;
  * i.e. Share of HVC routes compared to today's demand;
  * i.e. Year to consider for Eurostat reports;

* Simulation parameters : parameters impacting the optimization constraints and energy system definition

  * i.e. Temporal scale for the system optimization
  * i.e. Carbon budget per year (how much CO2 can be emitted annually);
  * i.e. Authorized expansion of AC/DC transmission lines (in terms of cost or transmission capacity);
  * i.e. Regionalized/copperplated ammonia at EU scale;
  * i.e. Emission pricing and sequestration costs per tCO2;
  * i.e. Locations where hydrogen storage is allowed;

Those additional parameters default values can be modified to match expert's best estimate.

Spatio-temporal specifications
---------------------------

PyPSA is technically able to define the energy supply down to a resolution of 1 hour and down to the spatial resolution of ENTSO-E transmission network. However, practically speaking, such a fine resolution (8760h on one year for ~8800 electrical nodes) is not feasible due to the huge computational burden linked to the optimization of such an energy system.

The system is hence clustered to a smaller number of equivalent electrical nodes  (i.e. clusters), small enough to allow acceptable runtimes but large enough to ensure a detailed representation of the energy system (power demand, renewable power generation, transmission infrastructures, etc.).

As mentioned in :cite:`frysztackiStrongEffect2021a`, we need to be especially be aware of the implications of those hypotheses. Model outputs are strongly influenced by network resolution. This is why we chose to take 37 clustered nodes into account while considering 180 renewables generation sites (onshore and offshore wind as well as utility-scale solar PV technologies). This gives a better estimation of the load factors for renewables without significantly increasing the computation time.

Temporal resolution has also been explored during the preliminary phase of the project. Two resolution techniques were proposed : time aggregation and time segmentation. Time aggregation averages timesteps on a given resolution (e.g.: 3h aggregation). Time segmentation use the `tsam` package (https://github.com/FZJ-IEK3-VSA/tsam). This package looks for typical periods using machine learning algorithms.  While having an impact on the computation time, we preferred a 3h time aggregation to be as close as possible to profiles. This choice also eases the interpretation of results.

More details about the spatial resolution are given in Section :ref:`spatial_resolution`.