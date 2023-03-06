# -*- coding: utf-8 -*-

from shutil import copy
from pathlib import Path

import yaml

files = {
    "config.yaml": "config.yaml",
    "Snakefile": "Snakefile",
    "scripts/solve_network.py": "solve_network.py",
    "scripts/prepare_sector_network.py": "prepare_sector_network.py",
}

if __name__ == "__main__":
    if "snakemake" not in globals():
        from helper import mock_snakemake

        snakemake = mock_snakemake("copy_config")

    basepath = Path(f"results/{snakemake.params.RDIR}configs/")

    for f, name in files.items():
        copy(f, basepath / name)

    with open(basepath / "config.snakemake.yaml", "w") as yaml_file:
        yaml.dump(
            snakemake.config,
            yaml_file,
            default_flow_style=False,
            allow_unicode=True,
            sort_keys=False,
        )
