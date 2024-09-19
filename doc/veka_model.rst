..
  SPDX-FileCopyrightText: 2019-2023 The PyPSA-Eur Authors

  SPDX-License-Identifier: CC-BY-4.0

.. _veka_modifications:

##########################################
Model modifications
##########################################


To better match client needs, various modifications have been added to PyPSA-Eur.
To implement this, two approaches have been considered :

1. when feasible, the request is generalized and submitted to the community as pull request (see :ref:`Submitted pull requests`) ;
2. if not, a custom modification is added to the model (see :ref:`Custom modifications`).

Submitted pull requests
===========================

Here is the list of pull requests shared with the community. Depending on the pull request, the modification has been integrated or not to the official repository of PyPSA-Eur.

- Define methanol energy demand for industry (https://github.com/PyPSA/pypsa-eur/pull/1068) ;
- Improve run_gurobi to wait for available token (https://github.com/PyPSA/linopy/pull/281) ;
- Fix non steel related coal demand during transition (using sector_ratios_fraction_future) (https://github.com/PyPSA/pypsa-eur/pull/1047) ;
- Add calculate_nodal_supply_energy in make summary (https://github.com/PyPSA/pypsa-eur/pull/1046) ;
- Fix gas network retrofit in brownfield (https://github.com/PyPSA/pypsa-eur/pull/1036) ;
- Improve `agg_p_nom_limits` configuration (https://github.com/PyPSA/pypsa-eur/pull/1023) ;
- Fix small issues in `add_land_use_constraint_m` (https://github.com/PyPSA/pypsa-eur/pull/1022) ;
- Clarify suffix usage in add existing baseyear (https://github.com/PyPSA/pypsa-eur/pull/1017) ;
- Fix custom busmap read in cluster network (https://github.com/PyPSA/pypsa-eur/pull/1008) ;
- Fix typo (https://github.com/PyPSA/pypsa-eur/pull/1005, https://github.com/PyPSA/pypsa-eur/pull/1045) ;
- Fix fill missing in industry sector ratios intermediate (https://github.com/PyPSA/pypsa-eur/pull/1004) ;
- Fix index for existing capacities in add_existing_baseyear (https://github.com/PyPSA/pypsa-eur/pull/992) ;
- Fix grouping year reference in add_land_use_constraint_m (https://github.com/PyPSA/pypsa-eur/pull/991) ;
- Fix error with symbol of buses in simplify_network (https://github.com/PyPSA/pypsa-eur/pull/987) ;
- Fix type error in cluster_network with "m" configuration (https://github.com/PyPSA/pypsa-eur/pull/986) ;
- Fix duplicated years in `add_land_use_constraint_m` (https://github.com/PyPSA/pypsa-eur/pull/968) ;
- Add decommissioning of renewables assets (https://github.com/PyPSA/pypsa-eur/pull/959) ;
- Add warning when negative bev availability profile values (https://github.com/PyPSA/pypsa-eur/pull/858) ;
- Fix typo in buses definition for oil boilers in `add_industry` in `prepare_sectors_networks` (https://github.com/PyPSA/pypsa-eur/pull/812) ;
- Fix nodal fraction with distributed generators (https://github.com/PyPSA/pypsa-eur/pull/798) ;
- Add option for SMR CC (https://github.com/PyPSA/pypsa-eur/pull/757) ;
- Add rule to update IRENA renewables capacities (https://github.com/PyPSA/pypsa-eur/pull/756) ;
- Improve Gurobi usage for `linopy` package (https://github.com/PyPSA/linopy/pull/162) ;
- Raised issue for `snakemake` package to better manage Gurobi licenses (https://github.com/snakemake/snakemake/issues/1801) ;
- Raised issue for `pulp` package to better manage Gurobi licenses (https://github.com/coin-or/pulp/issues/571).

Custom modifications
===========================

When no suitable for a pull request, the customisation is added on our fork of the model. Techno-economic assumptions have been thoroughly reviewed as described in :ref:`Techno-economic parameters`.


- Add project specific configurations ;
- Develop a `Streamlit <https://streamlit.io/>`_ app available at https://climact-veka-2050plus.streamlit.app/ ;
- Develop a data management pipeline to produce files required by the app ;
- Define a set of configurations and techno-economic assumptions, needed for scenario management ;
- Add custom configuration to model Flanders as a node ;
- Add custom configuration to model nuclear powerplants ;
- Fix PV shares in Belgium to reflect current imbalance between regions ;
- Split HVC primary route into three types of routes (NSC, NSC CC and MTO) to include VLAIO scenarios for Flanders ;
- Split DRI route of steel in two routes (DRI CH4 and DRI H2), needed for scenario management of industry. DRI CH4 is model using `MIDREX <https://www.sciencedirect.com/science/article/pii/S0973082623002132>`_ process ;
- Fix BEV dsm restriction time to take time aggregation and time zones into account ;
- Add emission prices for each planning horizon ;
- Update IRENA data to 2023 ;
- Update status of TYNDP links to the latest version available ;
- Limit offshore wind capacity in Belgium to a realistic value ;
- Add a proxy to model distribution costs of gas boilers ;
- Update OCGT, CCGT and coal CC efficiencies to include electricity demand of CC ;
- Add CC and H2 gas turbines ;
- Add hydrogen import terminals ;
- Add a marginal costs for gas and hydrogen pipelines to model operating costs ;
- Add a rule to produce map with capacities ;
- Add a rule to manage multiple set of costs, needed for scenario management ;
- Generalize `agg_p_nom_minmax` constraint to handle nuclear and all types of offshore wind ;
- Fix config provider to work as expected with reference scenario ;
- Stick to version 0.5.11 of https://github.com/PyPSA/powerplantmatching to ensure stability  ;
- Fix various minor issues.


