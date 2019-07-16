
import sys

sys.path = ["/home/vres/data/tom/lib/pypsa"] + sys.path

import numpy as np
import pandas as pd
import logging
logger = logging.getLogger(__name__)
import gc
import os

import pypsa

from pypsa.descriptors import free_output_series_dataframes

# Suppress logging of the slack bus choices
pypsa.pf.logger.setLevel(logging.WARNING)

from vresutils.benchmark import memory_logger





#First tell PyPSA that links can have multiple outputs by
#overriding the component_attrs. This can be done for
#as many buses as you need with format busi for i = 2,3,4,5,....
#See https://pypsa.org/doc/components.html#link-with-multiple-outputs-or-inputs


override_component_attrs = pypsa.descriptors.Dict({k : v.copy() for k,v in pypsa.components.component_attrs.items()})
override_component_attrs["Link"].loc["bus2"] = ["string",np.nan,np.nan,"2nd bus","Input (optional)"]
override_component_attrs["Link"].loc["bus3"] = ["string",np.nan,np.nan,"3rd bus","Input (optional)"]
override_component_attrs["Link"].loc["efficiency2"] = ["static or series","per unit",1.,"2nd bus efficiency","Input (optional)"]
override_component_attrs["Link"].loc["efficiency3"] = ["static or series","per unit",1.,"3rd bus efficiency","Input (optional)"]
override_component_attrs["Link"].loc["p2"] = ["series","MW",0.,"2nd bus output","Output"]
override_component_attrs["Link"].loc["p3"] = ["series","MW",0.,"3rd bus output","Output"]



def patch_pyomo_tmpdir(tmpdir):
    # PYOMO should write its lp files into tmp here
    import os
    if not os.path.isdir(tmpdir):
        os.mkdir(tmpdir)
    from pyutilib.services import TempfileManager
    TempfileManager.tempdir = tmpdir

def prepare_network(n, solve_opts=None):
    if solve_opts is None:
        solve_opts = snakemake.config['solving']['options']

    if 'clip_p_max_pu' in solve_opts:
        for df in (n.generators_t.p_max_pu, n.storage_units_t.inflow):
            df.where(df>solve_opts['clip_p_max_pu'], other=0., inplace=True)

    if solve_opts.get('load_shedding'):
        n.add("Carrier", "Load")
        n.madd("Generator", n.buses.index, " load",
               bus=n.buses.index,
               carrier='load',
               sign=1e-3, # Adjust sign to measure p and p_nom in kW instead of MW
               marginal_cost=1e2, # Eur/kWh
               # intersect between macroeconomic and surveybased
               # willingness to pay
               # http://journal.frontiersin.org/article/10.3389/fenrg.2015.00055/full
               p_nom=1e9 # kW
        )

    if solve_opts.get('noisy_costs'):
        for t in n.iterate_components():
            #if 'capital_cost' in t.df:
            #    t.df['capital_cost'] += 1e1 + 2.*(np.random.random(len(t.df)) - 0.5)
            if 'marginal_cost' in t.df:
                t.df['marginal_cost'] += 1e-2 + 2e-3*(np.random.random(len(t.df)) - 0.5)

        for t in n.iterate_components(['Line', 'Link']):
            t.df['capital_cost'] += (1e-1 + 2e-2*(np.random.random(len(t.df)) - 0.5)) * t.df['length']

    if solve_opts.get('nhours'):
        nhours = solve_opts['nhours']
        n.set_snapshots(n.snapshots[:nhours])
        n.snapshot_weightings[:] = 8760./nhours

    return n

def add_opts_constraints(n, opts=None):
    if opts is None:
        opts = snakemake.wildcards.opts.split('-')

    if 'BAU' in opts:
        mincaps = snakemake.config['electricity']['BAU_mincapacities']
        def bau_mincapacities_rule(model, carrier):
            gens = n.generators.index[n.generators.p_nom_extendable & (n.generators.carrier == carrier)]
            return sum(model.generator_p_nom[gen] for gen in gens) >= mincaps[carrier]
        n.model.bau_mincapacities = pypsa.opt.Constraint(list(mincaps), rule=bau_mincapacities_rule)

    if 'SAFE' in opts:
        peakdemand = (1. + snakemake.config['electricity']['SAFE_reservemargin']) * n.loads_t.p_set.sum(axis=1).max()
        conv_techs = snakemake.config['plotting']['conv_techs']
        exist_conv_caps = n.generators.loc[n.generators.carrier.isin(conv_techs) & ~n.generators.p_nom_extendable, 'p_nom'].sum()
        ext_gens_i = n.generators.index[n.generators.carrier.isin(conv_techs) & n.generators.p_nom_extendable]
        n.model.safe_peakdemand = pypsa.opt.Constraint(expr=sum(n.model.generator_p_nom[gen] for gen in ext_gens_i) >= peakdemand - exist_conv_caps)

def add_lv_constraint(n):
    line_volume = getattr(n, 'line_volume_limit', None)
    if line_volume is not None and not np.isinf(line_volume):
        n.model.line_volume_constraint = pypsa.opt.Constraint(
            expr=((sum(n.model.passive_branch_s_nom["Line",line]*n.lines.at[line,"length"]
                        for line in n.lines.index[n.lines.s_nom_extendable]) +
                    sum(n.model.link_p_nom[link]*n.links.at[link,"length"]
                        for link in n.links.index[(n.links.carrier=='DC') &
                                                    n.links.p_nom_extendable]))
                    <= line_volume)
        )

def add_eps_storage_constraint(n):
    if not hasattr(n, 'epsilon'):
        n.epsilon = 1e-5
    fix_sus_i = n.storage_units.index[~ n.storage_units.p_nom_extendable]
    n.model.objective.expr += sum(n.epsilon * n.model.state_of_charge[su, n.snapshots[0]] for su in fix_sus_i)

def add_battery_constraints(network):

    nodes = list(network.buses.index[network.buses.carrier == "battery"])

    def battery(model, node):
        return model.link_p_nom[node + " charger"] == model.link_p_nom[node + " discharger"]*network.links.at[node + " charger","efficiency"]

    network.model.battery = pypsa.opt.Constraint(nodes, rule=battery)



def add_chp_constraints(network):

    options = snakemake.config["sector"]

    if hasattr(network.links.index,"str") and network.links.index.str.contains("CHP").any():

        #AC buses with district heating
        urban_central = n.buses.index[n.buses.carrier == "urban central heat"]
        if not urban_central.empty:
            urban_central = urban_central.str[:-len(" urban central heat")]


        def chp_nom(model,node):
            return network.links.at[node + " urban central CHP electric","efficiency"]*options['chp_parameters']['p_nom_ratio']*model.link_p_nom[node + " urban central CHP electric"] == network.links.at[node + " urban central CHP heat","efficiency"]*options['chp_parameters']['p_nom_ratio']*model.link_p_nom[node + " urban central CHP heat"]


        network.model.chp_nom = pypsa.opt.Constraint(urban_central,rule=chp_nom)


        def backpressure(model,node,snapshot):
            return options['chp_parameters']['c_m']*network.links.at[node + " urban central CHP heat","efficiency"]*model.link_p[node + " urban central CHP heat",snapshot] <= network.links.at[node + " urban central CHP electric","efficiency"]*model.link_p[node + " urban central CHP electric",snapshot]

        network.model.backpressure = pypsa.opt.Constraint(urban_central,list(network.snapshots),rule=backpressure)


        def top_iso_fuel_line(model,node,snapshot):
            return model.link_p[node + " urban central CHP heat",snapshot] + model.link_p[node + " urban central CHP electric",snapshot] <= model.link_p_nom[node + " urban central CHP electric"]

        network.model.top_iso_fuel_line = pypsa.opt.Constraint(urban_central,list(network.snapshots),rule=top_iso_fuel_line)




def fix_branches(n, lines_s_nom=None, links_p_nom=None):
    if lines_s_nom is not None and len(lines_s_nom) > 0:
        for l, s_nom in lines_s_nom.iteritems():
            n.model.passive_branch_s_nom["Line", l].fix(s_nom)
        if isinstance(n.opt, pypsa.opf.PersistentSolver):
            n.opt.update_var(n.model.passive_branch_s_nom)

    if links_p_nom is not None and len(links_p_nom) > 0:
        for l, p_nom in links_p_nom.iteritems():
            n.model.link_p_nom[l].fix(p_nom)
        if isinstance(n.opt, pypsa.opf.PersistentSolver):
            n.opt.update_var(n.model.link_p_nom)

def solve_network(n, config=None, solver_log=None, opts=None):
    if config is None:
        config = snakemake.config['solving']
    solve_opts = config['options']

    solver_options = config['solver'].copy()
    if solver_log is None:
        solver_log = snakemake.log.solver
    solver_name = solver_options.pop('name')

    def run_lopf(n, allow_warning_status=False, fix_zero_lines=False, fix_ext_lines=False):
        free_output_series_dataframes(n)

        if not hasattr(n, 'opt') or not isinstance(n.opt, pypsa.opf.PersistentSolver):
            pypsa.opf.network_lopf_build_model(n, formulation=solve_opts['formulation'])
            add_opts_constraints(n, opts)
            add_lv_constraint(n)
            #add_eps_storage_constraint(n)
            add_battery_constraints(n)
            add_chp_constraints(n)

            pypsa.opf.network_lopf_prepare_solver(n, solver_name=solver_name)

        if fix_zero_lines:
            fix_lines_b = (n.lines.s_nom_opt == 0.) & n.lines.s_nom_extendable
            fix_links_b = (n.links.carrier=='DC') & (n.links.p_nom_opt == 0.) & n.links.p_nom_extendable
            fix_branches(n,
                         lines_s_nom=pd.Series(0., n.lines.index[fix_lines_b]),
                         links_p_nom=pd.Series(0., n.links.index[fix_links_b]))

        if fix_ext_lines:
            fix_branches(n,
                         lines_s_nom=n.lines.loc[n.lines.s_nom_extendable, 's_nom_opt'],
                         links_p_nom=n.links.loc[(n.links.carrier=='DC') & n.links.p_nom_extendable, 'p_nom_opt'])


        if not fix_ext_lines and hasattr(n.model, 'line_volume_constraint'):

            def extra_postprocessing(n, snapshots, duals):
                index = list(n.model.line_volume_constraint.keys())
                cdata = pd.Series(list(n.model.line_volume_constraint.values()),
                                  index=index)
                n.line_volume_limit_dual =  -cdata.map(duals).sum()
                print("line volume limit dual:",n.line_volume_limit_dual)

        else:
            extra_postprocessing = None

        # Firing up solve will increase memory consumption tremendously, so
        # make sure we freed everything we can
        gc.collect()

        #from pyomo.opt import ProblemFormat
        #print("Saving model to MPS")
        #n.model.write('/home/ka/ka_iai/ka_kc5996/projects/pypsa-eur/128-B-I.mps', format=ProblemFormat.mps)
        #print("Model is saved to MPS")
        #sys.exit()


        status, termination_condition = \
        pypsa.opf.network_lopf_solve(n,
                                     solver_logfile=solver_log,
                                     solver_options=solver_options,
                                     formulation=solve_opts['formulation'],
                                     extra_postprocessing=extra_postprocessing
                                     #free_memory={'pypsa'}
                                     )

        assert status == "ok" or allow_warning_status and status == 'warning', \
            ("network_lopf did abort with status={} "
             "and termination_condition={}"
             .format(status, termination_condition))

        return status, termination_condition

    lines_ext_b = n.lines.s_nom_extendable
    if lines_ext_b.any():
        # puh: ok, we need to iterate, since there is a relation
        # between s/p_nom and r, x for branches.
        msq_threshold = 0.01
        lines = pd.DataFrame(n.lines[['r', 'x', 'type', 'num_parallel']])

        lines['s_nom'] = (
            np.sqrt(3) * n.lines['type'].map(n.line_types.i_nom) *
            n.lines.bus0.map(n.buses.v_nom)
        ).where(n.lines.type != '', n.lines['s_nom'])

        lines_ext_typed_b = (n.lines.type != '') & lines_ext_b
        lines_ext_untyped_b = (n.lines.type == '') & lines_ext_b

        def update_line_parameters(n, zero_lines_below=10, fix_zero_lines=False):
            if zero_lines_below > 0:
                n.lines.loc[n.lines.s_nom_opt < zero_lines_below, 's_nom_opt'] = 0.
                n.links.loc[(n.links.carrier=='DC') & (n.links.p_nom_opt < zero_lines_below), 'p_nom_opt'] = 0.

            if lines_ext_untyped_b.any():
                for attr in ('r', 'x'):
                    n.lines.loc[lines_ext_untyped_b, attr] = (
                        lines[attr].multiply(lines['s_nom']/n.lines['s_nom_opt'])
                    )

            if lines_ext_typed_b.any():
                n.lines.loc[lines_ext_typed_b, 'num_parallel'] = (
                    n.lines['s_nom_opt']/lines['s_nom']
                )
                logger.debug("lines.num_parallel={}".format(n.lines.loc[lines_ext_typed_b, 'num_parallel']))

            if isinstance(n.opt, pypsa.opf.PersistentSolver):
                n.calculate_dependent_values()

                assert solve_opts['formulation'] == 'kirchhoff', \
                    "Updating persistent solvers has only been implemented for the kirchhoff formulation for now"

                n.opt.remove_constraint(n.model.cycle_constraints)
                del n.model.cycle_constraints_index
                del n.model.cycle_constraints_index_0
                del n.model.cycle_constraints_index_1
                del n.model.cycle_constraints

                pypsa.opf.define_passive_branch_flows_with_kirchhoff(n, n.snapshots, skip_vars=True)
                n.opt.add_constraint(n.model.cycle_constraints)

        iteration = 1

        lines['s_nom_opt'] = lines['s_nom'] * n.lines['num_parallel'].where(n.lines.type != '', 1.)
        status, termination_condition = run_lopf(n, allow_warning_status=True)

        def msq_diff(n):
            lines_err = np.sqrt(((n.lines['s_nom_opt'] - lines['s_nom_opt'])**2).mean())/lines['s_nom_opt'].mean()
            logger.info("Mean square difference after iteration {} is {}".format(iteration, lines_err))
            return lines_err

        min_iterations = solve_opts.get('min_iterations', 2)
        max_iterations = solve_opts.get('max_iterations', 999)

        while msq_diff(n) > msq_threshold or iteration < min_iterations:
            if iteration >= max_iterations:
                logger.info("Iteration {} beyond max_iterations {}. Stopping ...".format(iteration, max_iterations))
                break

            update_line_parameters(n)
            lines['s_nom_opt'] = n.lines['s_nom_opt']
            iteration += 1

            status, termination_condition = run_lopf(n, allow_warning_status=True)

        update_line_parameters(n, zero_lines_below=100)

        logger.info("Starting last run with fixed extendable lines")

        # Not really needed, could also be taken out
        # if 'snakemake' in globals():
        #     fn = os.path.basename(snakemake.output[0])
        #     n.export_to_netcdf('/home/vres/data/jonas/playground/pypsa-eur/' + fn)

    status, termination_condition = run_lopf(n, fix_ext_lines=True)

    # Drop zero lines from network
    # zero_lines_i = n.lines.index[(n.lines.s_nom_opt == 0.) & n.lines.s_nom_extendable]
    # if len(zero_lines_i):
    #     n.mremove("Line", zero_lines_i)
    # zero_links_i = n.links.index[(n.links.p_nom_opt == 0.) & n.links.p_nom_extendable]
    # if len(zero_links_i):
    #     n.mremove("Link", zero_links_i)


    return n

if __name__ == "__main__":
    # Detect running outside of snakemake and mock snakemake for testing
    if 'snakemake' not in globals():
        from vresutils.snakemake import MockSnakemake, Dict
        snakemake = MockSnakemake(
            wildcards=dict(network='elec', simpl='', clusters='45', lv='1.25', opts='Co2L-3H-T-H'),
            input=["networks/{network}_s{simpl}_{clusters}_lv{lv}_{opts}.nc"],
            output=["results/networks/s{simpl}_{clusters}_lv{lv}_{opts}-test.nc"],
            log=dict(gurobi="logs/{network}_s{simpl}_{clusters}_lv{lv}_{opts}_gurobi-test.log",
                     python="logs/{network}_s{simpl}_{clusters}_lv{lv}_{opts}_python-test.log")
        )


    tmpdir = snakemake.config['solving'].get('tmpdir')
    if tmpdir is not None:
        patch_pyomo_tmpdir(tmpdir)

    logging.basicConfig(filename=snakemake.log.python,
                        level=snakemake.config['logging_level'])

    with memory_logger(filename=getattr(snakemake.log, 'memory', None), interval=30.) as mem:
        n = pypsa.Network(snakemake.input[0],
                          override_component_attrs=override_component_attrs)

        n = prepare_network(n)
        n = solve_network(n)

        n.export_to_netcdf(snakemake.output[0])

    logger.info("Maximum memory usage: {}".format(mem.mem_usage))
