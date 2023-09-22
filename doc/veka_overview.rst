..
  SPDX-FileCopyrightText: 2019-2023 The PyPSA-Eur Authors

  SPDX-License-Identifier: CC-BY-4.0

.. _veka_overview:

##########################################
Overview
##########################################


This chapter aims to give a broad overview of the use of PyPSA-Eur in the context of the *Energie Studie 2050+* for the *Flemish Energy and Climate Agency* (VEKA, Belgium).

High level description
===========================
PyPSA-Eur is an open model dataset of the European energy system at the transmission network level that covers the full ENTSO-E area. Written in Python, it optimizes the energy system with high temporal and spatial resolution for a large spectrum of technologies and energy vectors. It covers demand and supply for all energy sectors.

General description
---------------------------

PyPSA-Eur optimizes energy systems with

- a temporal granularity up to 1 hour (which represents a maximum of 8760 time steps each year),
- a spatial granularity up to ENTSO-E transmission lines (which represents a maximum of 8800 electrical buses), and
- a technology list reaching up to 80 different technologies (see :ref:`Technological assumptions`).

PyPSA-Eur allows connecting different energy vectors (referred to as energy carriers in PyPSA-Eur) and technologies using them, hence allowing to tackle inter-sectorial and inter-carriers connections.

The system is optimized over the energy system cost (including OPEX and annualized CAPEX), for different consecutive planning horizons. The optimization is done over the sizing and dispatch of energy production and storage units and of transmission infrastructures to meet an inelastic demand for each of the time frames and geographical locations considered.

The cost optimization of the energy system for a given year is performed under different constraints, including an annual CO2 emission constraint for the given year and renewable potentials.

Optimization details
---------------------------
For each planning horizon, PyPSA-Eur provides the cost optimal energy dispatch of each asset at each node and for each time frame given an hourly demand at each node of each energy carrier, as well as the cost optimal sizing of production, storage and transmission units needed for such dispatch (see `Power System Optimization <https://pypsa.readthedocs.io/en/latest/optimal_power_flow.html#power-system-optimization>`_)

The optimization hence minimizes the annual total cost of the energy system for each planning horizon, which is defined as :

.. math::

    c = \sum_{n,s}{CAPEX_{n,s}} + \sum_{l}{CAPEX_{l}} + \sum_{t}{w_t \cdot \left( \sum_{n,s}OPEX_{n,s,t}\right)} + \sum_{t} \left[suc_{n,s,t} + sdc_{n,s,t} \right]

where :

* :math:`CAPEX_{n,s}` is the annualized investment cost of units at node *n* for generator or storage asset *s*
* :math:`CAPEX_{l}` is the annualized investment cost of infrastructure on line *l*
* :math:`OPEX_{n,s,t}` is the operational cost of units at node *n*, for generator or storage asset *s* at time frame *t*
* :math:`w_{t}` is the weighting of time *t* in the objective function (e.g. multiple hours)
* :math:`suc_{n,s,t}` is the start-up cost if generator with unit commitment is started at time :math:`t`
* :math:`sdc_{n,s,t}` is the shut-down cost if generator with unit commitment is shut down at time :math:`t`

If *myopic* optimization is configured, the solver does not optimize the energy system over a continuous trajectory but rather planning horizon by planning horizon. For each planning horizon, the capacity installed in the previous planning horizon is taken into account and phased out assets are removed.

Given the spatial and temporal resolution used, the use of a commercial solver is required to compute the solution in a reasonable amount of time. Therefore, a `Gurobi <https://www.gurobi.com/>`_ license, which is the standard choice in the PyPSA community, has been bought.

Optimization constraints
---------------------------
The energy system is optimized under various constrains, among which :

* A carbon budget constraint, representing how much CO2 can be emitted over the operation of the energy system for the considered planning horizons
* A potential constraint, representing the maximal capacity that can be installed per technology and per node. Currently, renewables are the only technologies for which this potential constraint applies.


Main limitations
===========================

The current study has been developed on the version v0.8.0 of PyPSA (https://github.com/PyPSA/pypsa-eur/releases/tag/v0.8.0). Please refer to the `release notes <https://pypsa-eur.readthedocs.io/en/latest/release_notes.html>`_ for further details on fixes and improvements since then. A `limitations <https://pypsa-eur.readthedocs.io/en/latest/limitations.html>`_ list is also maintained by the PyPSA-Eur team. During the development of this project, the following attention points have been identified :

* Hydrogen, CH4 and CO2 pipelines are considered to be lossless and free of electricity consumption.

* Existing conventional generators assets considered for the first planning horizons are extracted from *powerplantmatching* Python package (https://github.com/PyPSA/powerplantmatching), which maintained by TUBerlin.

* Existing renewable generators assets considered for the first planning horizons come from `IRENA <http://pxweb.irena.org/pxweb/en/IRENASTAT>`_.

* The following technologies are not yet supported by the model:

  * geothermal, and
  * industrial heat.

* Export are not taken into account in the current version of the model. The only way to model it is exogenous.

* Hydrogen can only be

  * used as a feedstock/energy carrier for industry and as an energy carrier for different,
  * produced and is not imported, or
  * transported in a gaseous form through pipelines, or
  * used as liquefied hydrogen or hydrogen for shipping demand.