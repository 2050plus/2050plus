..
  SPDX-FileCopyrightText: 2019-2023 The PyPSA-Eur Authors

  SPDX-License-Identifier: CC-BY-4.0

.. _veka_configurations:

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

- Update OCGT, CCGT and coal CC efficiencies to include electricity demand of CC
- Add hydrogen import terminals
- Limit offshore wind capacity in Belgium to a realistic value
- Update status of TYNDP links to the latest version available
- Update industry modelling to include VLAIO scenarios for Flanders
- Fix various issues
- Generalize `agg_p_nom_minmax` constraint to handle nuclear (https://github.com/2050plus/2050plus/commit/13de195661be946ced134bdc0d846b7307b3b886) ;
- Add `plot_capacities` to summary rules (https://github.com/2050plus/2050plus/commit/f840947c0dbd505660fabc9e6d37156ec09cf626) ;
- Stick to version 0.5.11 of https://github.com/PyPSA/powerplantmatching to ensure stability (https://github.com/2050plus/2050plus/commit/704c187df614aa1fca6b6d2e0e2b59384996433a) ;
- Develop a data management pipeline to produce files required by the app (https://github.com/2050plus/2050plus/commit/e789cd4a9677118f6e9172799cafd56b6870739b, https://github.com/2050plus/2050plus/commit/5ea351e270edf61097ca004a90948657e9977d9f) ;
- Develop a `Streamlit <https://streamlit.io/>`_ app available at https://climact-veka-2050plus.streamlit.app/ (https://github.com/2050plus/2050plus/commit/20d2f4138aec1fa8392f93f261828551dbdde1dc, https://github.com/2050plus/2050plus/commit/ac62b2d6bd48c58891467b1cad2191160013674e, https://github.com/2050plus/2050plus/commit/2cf00459fbf469b2d77d2d653f777e36f2d28273, https://github.com/2050plus/2050plus/commit/c089dca3dda4a7db095e3a924afdebdb1fecab11, https://github.com/2050plus/2050plus/commit/e2a231732c15973e8e79e365d6b382c1b95f2486, https://github.com/2050plus/2050plus/commit/c2fe92168d4b3e02307c4a38e021a4b9294313af) ;
- Add emission prices for each planning horizon (https://github.com/2050plus/2050plus/commit/1964af211dc70afd5a583eb2de8d317d982f5494) ;
- Split DRI route of steel in two routes (DRI CH4 and DRI H2), needed for scenario management of industry. DRI CH4 is model using `MIDREX <https://www.sciencedirect.com/science/article/pii/S0973082623002132>`_ process. (https://github.com/2050plus/2050plus/commit/ddb5e77c53957416d95390a3fdf4a055f94bf708) ;
- Add CC and H2 gas turbines (https://github.com/2050plus/2050plus/commit/0d831f810cc22df9e7a6d3b9107532a1b7b7da53) ;
- Add a marginal costs for gas and hydrogen pipelines to model operating costs (https://github.com/2050plus/2050plus/commit/17c64e008455aad988038f2034615e35da2ff295, https://github.com/2050plus/2050plus/commit/780612e1a779fe74c3c38f46400cc955ccbc2915) ;
- Define a set of configurations and techno-economic assumptions, needed for scenario management (https://github.com/2050plus/2050plus/commit/98d36d779c9561d04984ddc46c14738fa479a274, https://github.com/2050plus/2050plus/commit/e5f151f5ca39633e7c0a45614e29a859158117e4) ;
- Add custom configuration to model Flanders as a node (https://github.com/2050plus/2050plus/commit/209aaca2986c6549505faad5875ab985076a42f9) ;
- Add custom configuration to model nuclear powerplants (https://github.com/2050plus/2050plus/commit/822069e7e86b4a061013c27cd2056bd97719c4c8) ;
- Add a rule to manage multiple set of costs, needed for scenario management (https://github.com/2050plus/2050plus/commit/1c2bf6b7bdcf35d87b801c5123d95b84dff6383e) ;
- Add default configuration for project (https://github.com/2050plus/2050plus/commit/ed216fccaa2b1d55451c5b548d9bf4ec09bfff46, https://github.com/2050plus/2050plus/commit/e9f6e6505577c16fa42ad2970085e3822e9e7faf).


