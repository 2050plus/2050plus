# SPDX-FileCopyrightText: : 2017-2023 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT

from os.path import normpath, exists
from shutil import copyfile, move, rmtree

from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider

HTTP = HTTPRemoteProvider()

from snakemake.utils import min_version

min_version("7.7")

conf_file = os.path.join(workflow.current_basedir, "config/config.yaml")
conf_default_file = os.path.join(workflow.current_basedir, "config/config.default.yaml")
if not exists(conf_file) and exists(conf_default_file):
    copyfile(conf_default_file, conf_file)


configfile: "config/config.yaml"


COSTS = f"data/costs_{config['costs']['year']}.csv"
ATLITE_NPROCESSES = config["atlite"].get("nprocesses", 4)

run = config.get("run", {})
RDIR = run["name"] + "/" if run.get("name") else ""
CDIR = RDIR if not run.get("shared_cutouts") else ""

LOGS = "logs/" + RDIR
BENCHMARKS = "benchmarks/" + RDIR
RESOURCES = "resources/" + RDIR if not run.get("shared_resources") else "resources/"
RESULTS = "results/" + RDIR


localrules:
    purge,


wildcard_constraints:
    simpl="[a-zA-Z0-9]*",
    clusters="[0-9]+(m|c)?|all",
    ll="(v|c)([0-9\.]+|opt)",
    opts="[-+a-zA-Z0-9\.]*",
    sector_opts="[-+a-zA-Z0-9\.\s]*",


include: "rules/common.smk"
include: "rules/collect.smk"
include: "rules/retrieve.smk"
include: "rules/build_electricity.smk"
include: "rules/build_sector.smk"
include: "rules/solve_electricity.smk"
include: "rules/postprocess.smk"
include: "rules/validate.smk"


if config["foresight"] == "overnight":

    include: "rules/solve_overnight.smk"


if config["foresight"] == "myopic":

    include: "rules/solve_myopic.smk"


if config["foresight"] == "perfect":

    include: "rules/solve_perfect.smk"


rule all:
    input:
        RESULTS + "graphs/costs.pdf",
        expand(
            RESULTS + "maps/power-network-s{simpl}-{clusters}.pdf",
            **config["scenario"]
        ),
        expand(
            RESULTS
            + "maps/elec_s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}-costs-all_{planning_horizons}.pdf",
            **config["scenario"]
        ),
        expand(
            RESULTS
            + "maps/elec_s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}-h2_network_{planning_horizons}.pdf",
            **config["scenario"]
        ),
        expand(
            RESULTS
            + "maps/elec_s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}-ch4_network_{planning_horizons}.pdf",
            **config["scenario"]
        ),
    default_target: True


rule purge:
    run:
        import builtins

        do_purge = builtins.input(
            "Do you really want to delete all generated resources, \nresults and docs (downloads are kept)? [y/N] "
        )
        if do_purge == "y":
            rmtree("resources/", ignore_errors=True)
            rmtree("results/", ignore_errors=True)
            rmtree("doc/_build", ignore_errors=True)
            print("Purging generated resources, results and docs. Downloads are kept.")
        else:
            raise Exception(f"Input {do_purge}. Aborting purge.")


rule dag:
    message:
        "Creating DAG of workflow."
    output:
        dot=RESOURCES + "dag.dot",
        pdf=RESOURCES + "dag.pdf",
        png=RESOURCES + "dag.png",
    conda:
        "envs/environment.yaml"
    shell:
        """
        snakemake --rulegraph all | sed -n "/digraph/,\$p" > {output.dot}
        dot -Tpdf -o {output.pdf} {output.dot}
        dot -Tpng -o {output.png} {output.dot}
        """


rule doc:
    message:
        "Build documentation."
    output:
        directory("doc/_build"),
    shell:
        "make -C doc html"


rule sync:
    params:
        cluster=f"{config['remote']['ssh']}:{config['remote']['path']}",
    shell:
        """
        rsync -uvarh --ignore-missing-args --files-from=.sync-send . {params.cluster}
        rsync -uvarh --no-g {params.cluster}/resources . || echo "No resources directory, skipping rsync"
        rsync -uvarh --no-g {params.cluster}/results . || echo "No results directory, skipping rsync"
        rsync -uvarh --no-g {params.cluster}/logs . || echo "No logs directory, skipping rsync"
        """
