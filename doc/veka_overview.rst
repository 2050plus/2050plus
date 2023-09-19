..
  SPDX-FileCopyrightText: 2019-2023 The PyPSA-Eur Authors

  SPDX-License-Identifier: CC-BY-4.0

.. _veka_overview:

##########################################
Overview
##########################################


PyPSA-Eur is an open-source tool that allows the optimization of energy system which high temporal and spatial granularity for a large spectrum of technologies.

	- Temporal granularity can go up from 1year to 1hour (which represents a maximum of 8760 timesteps per year)
	- Spatial granularity can go up to ENTSO-E transmission lines (which represents a maximum of 8800 electrical buses)
	- Technology list can be as long as 80 different technologies

The optimization is done over the sizing and dispatch of production, storage and transmission energy units to meet an inelastic demand for each of the time frames and geographical locations considered.

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

Spacio-temporal specifications
===========================

PyPSA is technically able to define the energy supply down to a resolution of 1hour and down to a spatial resolution of ENTSO-E transmission network. However, pratically speaking, such a fine resolution (8760h for ~6000 electrical nodes) is not possible due to the enormous computational burden linked to the optimization of such an enormous energy system. A trade-off has hence to be found.


High level model description
===========================

Lorem ipsum dolor sit amet. Eos excepturi veniam et illo vitae qui nihil libero. Ut modi commodi et porro doloremque id molestiae laborum.

General description
---------------------------

Lorem ipsum dolor sit amet. Eos excepturi veniam et illo vitae qui nihil libero. Ut modi commodi et porro doloremque id molestiae laborum.

Optimization details
---------------------------

Lorem ipsum dolor sit amet. Eos excepturi veniam et illo vitae qui nihil libero. Ut modi commodi et porro doloremque id molestiae laborum.