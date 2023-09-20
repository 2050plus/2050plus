..
  SPDX-FileCopyrightText: 2019-2023 The PyPSA-Eur Authors

  SPDX-License-Identifier: CC-BY-4.0

.. _veka_overview:

##########################################
Overview
##########################################



High level model description
===========================

PyPSA-Eur is an open-source tool that allows the optimization of energy system which high temporal and spatial granularity for a large spectrum of technologies and energy vectors.
General description
---------------------------

PyPSA-Eur optimizes energy systems with : 

	- A temporal granularity up to 1 hour (which represents a maximum of 8760 timesteps per year) ;
	- A spatial granularity up to ENTSO-E transmission lines (which represents a maximum of 8800 electrical buses) ;
	- A technology list as long as 80 different technologies

PyPSA-Eur allows to connect different energy vectors (referred to as energy carriers in PyPSA-Eur) and technologies using them, hence allowing to tackle inter-sectorial and inter-carriers connections.

The system is optimized over the energy system cost (including OPEX and annualized CAPEX), for different consecutive planning horizons. The optimization is done over the sizing and dispatch of energy production and storage units and of transmission infrastructures to meet an inelastic demand for each of the time frames and geographical locations considered.

The cost optimization of the energy system for a given year is performed under different constraints, including an annual CO2 emission constraint for the given year and renewable potentials.

Optimization details
---------------------------
For each planning horizon, PyPSA-Eur provides the cost optimal energy dispacth of each asset at each node and for each time frame given an hourly demand at each node of each energy carrier, as well as the cost optimal sizing of production, storage and transmission units needed for such dispatch.

The optimization hence minimizes the annual total cost of the energy system for each planning horizon, which is defined as :

.. math::

    c = \sum{CAPEX_{n,s}}_{n,s} + \sum{CAPEX_{l,s}}_{l,s} + \sum{OPEX_{n,s,t}}_{n,s}

where :
	- $CAPEX_{n,s}$ is the annualized investment cost of units at node $n$ and for generator or storage asset $s$ ; 
	- $CAPEX_{l,s}$ is the annualized investment cost of infrastructure on line $l$ and for transmission asset $s$ ; 
	- $OPEX_{n,s,t}$ is the operational cost of units at node $n$, for generator or storage asset $s$ at time frame $t$; 
	
For each planning horizon, the capacity installed in previous planning horizon is taken into account and phased out assets are removed. Hence, the optimization is considered as "myopic" as it does not optimise the energy system over a continuous trajectory but rather planning horizon by planning horizon. 

Optimization constraints
---------------------------
The energy system is optimized under :
	- A carbon budget constraint, representing how much CO2 can be emitted over the operation of the energy system for the considered planning horizons
	- A potential constraint, representing the maximal capacity that can be installed per technology and per node. Currently, renewables are the only technologies for which this potential constraint applies.


Spatio-temporal specifications
---------------------------

PyPSA is technically able to define the energy supply down to a resolution of 1 hour and down to the spatial resolution of ENTSO-E transmission network.

However, pratically speaking, such a fine resolution (8760h for ~8800 electrical nodes) is not possible due to the enormous computational burden linked to the optimization of such an enormous energy system. 

The system is hence clustered to a smaller number of electrical nodes, small enough to allow acceptable runtimes but large enough to ensure a detailled representation of the energy system (power demand, renewable power generation, transmission infrastructures, etc).

Main limitations
===========================

As for 15 August 2023 v0.8.0
	- Hydrogen, CH4 and CO2 pipelines are considered to be lossless and free of electricity consumption
	- PyPSA could underestimate ren potentials because …
	- the following technologies are not yet supported by the model
		○ geothermal
		○ Industrial heat
	- Export are not taken into account in the current version of the model
	- Hydrogen can only be
		○ Used as a feedstock/energy carrier for industry and as an energy carrier for different
		○ Produced and is not imported
		○ Transported in a gaseous form through pipelines
Used as liquified hydrogen or hydrogen for shipping demand
