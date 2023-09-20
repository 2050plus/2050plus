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

PyPSA-Eur optimizes energy systems over 

- Temporal granularity up to 1 hour (which represents a maximum of 8760 timesteps per year)
- Spatial granularity can go up to ENTSO-E transmission lines (which represents a maximum of 8800 electrical buses)
- Technology list can be as long as 80 different technologies

PyPSA-Eur allows to connect different energy vectors (referred to as energy carriers in PyPSA-Eur) and technologies using them, hence allowing to tackle inter-sectorial and inter-carriers connections.


Optimization details
---------------------------

The system is optimized over the energy system cost (including OPEX and annualized CAPEX), for different consecutive planning horizons. The optimization is done over the sizing and dispatch of energy production and storage units and of transmission infrastructures to meet an inelastic demand for each of the time frames and geographical locations considered.

The cost optimization of the energy system for a given year is performed under different constraints, including an annual CO2 emission constraint for the given year and renewable potentials.

Spatio-temporal specifications
---------------------------

PyPSA is technically able to define the energy supply down to a resolution of 1 hour and down to the spatial resolution of ENTSO-E transmission network.

However, pratically speaking, such a fine resolution (8760h for ~8000 electrical nodes) is not possible due to the enormous computational burden linked to the optimization of such an enormous energy system. 

The system is hence clustered to a smaller number of electrical nodes, small enough to allow acceptable runtimes but large enough to ensure a detailled representation of the energy system (power demand, renewable power generation, transmission infrastructures, etc).





Main limitations
===========================


The current study has been developed on the version v0.8.0 of PyPSA (https://github.com/PyPSA/pypsa-eur/releases/tag/v0.8.0). During the development of the project, the following limitations have been identified :

* Hydrogen, CH4 and CO2 pipelines are considered to be lossless and free of electricity consumption.

* PyPSA could underestimate renewables potentials because

* The following technologies are not yet supported by the model:

  * geothermal, and
  * industrial heat.

* Export are not taken into account in the current version of the model.

* Hydrogen can only be

  * used as a feedstock/energy carrier for industry and as an energy carrier for different,
  * produced and is not imported, or
  * transported in a gaseous form through pipelines, or
  * used as liquefied hydrogen or hydrogen for shipping demand.
