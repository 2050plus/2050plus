..
  SPDX-FileCopyrightText: 2019-2023 The PyPSA-Eur Authors

  SPDX-License-Identifier: CC-BY-4.0

.. _veka_overview:

##########################################
Overview
##########################################


This chapter aims to give a broad overview of the use of PyPSA in the context of the *Energie Studie 2050+* for the *Flemish Energy and Climate Agency* (VEKA, Belgium).

High level description
===========================
PyPSA-Eur is an open model dataset of the European energy system at the transmission network level that covers the full ENTSO-E area. Written in Python, it optimizes the energy system with high temporal and spatial resolution for a large spectrum of technologies and energy vectors. It covers demand and supply for all energy sectors.

General description
---------------------------

PyPSA-Eur optimizes energy systems with

- a temporal granularity up to 1 hour (which represents a maximum of 8760 timesteps each year),
- a spatial granularity up to ENTSO-E transmission lines (which represents a maximum of 8800 electrical buses), and
- a technology list reaching up to 80 different technologies.

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

* :math:`CAPEX_{n,s}` is the annualized investment cost of units at node *n* and for generator or storage asset *s* ;
* :math:`CAPEX_{l,s}` is the annualized investment cost of infrastructure on line *l* and for transmission asset *s* ;
* :math:`OPEX_{n,s,t} is the operational cost of units at node *n*, for generator or storage asset *s* at time frame *t*;
	
For each planning horizon, the capacity installed in previous planning horizon is taken into account and phased out assets are removed. Hence, the optimization is considered as "myopic" as it does not optimise the energy system over a continuous trajectory but rather planning horizon by planning horizon. 

Optimization constraints
---------------------------
The energy system is optimized under :

* A carbon budget constraint, representing how much CO2 can be emitted over the operation of the energy system for the considered planning horizons
* A potential constraint, representing the maximal capacity that can be installed per technology and per node. Currently, renewables are the only technologies for which this potential constraint applies.


Spatio-temporal specifications
---------------------------

PyPSA is technically able to define the energy supply down to a resolution of 1 hour and down to the spatial resolution of ENTSO-E transmission network. However, practically speaking, such a fine resolution (8760h on one year for ~8800 electrical nodes) is not feasible due to the huge computational burden linked to the optimization of such an energy system.

The system is hence clustered to a smaller number of equivalent electrical nodes  (i.e. clusters), small enough to allow acceptable runtimes but large enough to ensure a detailed representation of the energy system (power demand, renewable power generation, transmission infrastructures, etc).

As mentioned in :cite:`frysztackiStrongEffect2021a`, we need to be especially be aware of the implications of those hypothesis. Model outputs are strongly influenced by network resolution. This is why we chose to take 37 clustered nodes into account while considering 180 renewables generation sites (onshore and offshore wind as well as utility-scale solar PV technologies). This gives a better estimation of the load factors for renewables without significantly increasing the computation time.

Temporal resolution has also been explored during the preliminary phase of the project. Two resolution techniques were proposed : time aggregation and time segmentation. Time aggregation averages timesteps on a given resolution (e.g.: 3h aggregation). Time segmentation use the `tsam` package (https://github.com/FZJ-IEK3-VSA/tsam). This package looks for typical periods using machine learning algorithms.  While having an impact on the computation time, we preferred a 3h time aggregation to be as close as possible to profiles. This choice eases also the interpretation of results.

Main limitations
===========================

The current study has been developed on the version v0.8.0 of PyPSA (https://github.com/PyPSA/pypsa-eur/releases/tag/v0.8.0). During the development of the project, the following limitations have been identified :

* Hydrogen, CH4 and CO2 pipelines are considered to be lossless and free of electricity consumption.

* Existing conventional generators assets considered for the first planning horizons are extracted from powerplantmatching Python package, which maintained by TUBerlin.

* Existing renewable generators assets considered for the first planning horizons come from IRENA

* The following technologies are not yet supported by the model:

  * geothermal, and
  * industrial heat.

* Export are not taken into account in the current version of the model.

* Hydrogen can only be

  * used as a feedstock/energy carrier for industry and as an energy carrier for different,
  * produced and is not imported, or
  * transported in a gaseous form through pipelines, or
  * used as liquefied hydrogen or hydrogen for shipping demand.