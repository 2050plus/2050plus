# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2024 Climact for The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT

"""
Create data ready to present (main file).
See also :
    - graph_extraction_extract
    - graph_extraction_transform
    - graph_extraction_load
    - graph_extraction_utils
"""
import logging
from itertools import repeat
from pathlib import Path

from joblib import Parallel
from joblib import delayed
from tqdm import tqdm

from scripts.graph_extraction_extract import extract_data
from scripts.graph_extraction_load_st import load_data_st
from scripts.graph_extraction_transform import transform_data
from scripts.graph_extraction_utils import load_config

logger = logging.getLogger(__name__)


def compute_scenario_data(run, scenario):
    logger.info(f"Start processing of {scenario} ({run})")

    config_file = "config/config.veka.yaml"
    analysis_path = Path("analysis", run)
    dir_export = "graph_data"

    # Configuration
    config = load_config(config_file, analysis_path, dir_export, scenario=scenario)

    config["eu27_countries"] = ["AT", "BG", "BE", "CY", "CZ", "DE", "DK", "EE", "GR", "ES", "FI", "FR", "HR",
                                "HU", "IE", "IT", "LT", "LU", "LV", "MT", "NL", "PL", "PT", "RO", "SE", "SI", "SK", ]
    config["eu27_countries"] = list(set(config["eu27_countries"]).intersection(set(config["countries"])))
    # global variables for which to do work
    config["countries"] = {"tot": None, "eu27": config["eu27_countries"], "be": ["BE"], "fl": ["FL"],
                           "fr": ["FR"], "de": ["DE"], "gb": ["GB"], "nl": ["NL"], "lu": ["LU"]}
    config["imp_exp_carriers"] = ["elec", "gas", "H2"]

    # Extract data
    n, n_ext, context = extract_data(config, )

    # Transform data
    transform_data(config, n, n_ext)

    # Load data
    load_data_st(config, context)

    logger.info(f"Done for {scenario} ({run})")


def main():
    scenarios = list(zip(repeat("20240619"), ["central", "electrification", "molecules", "lsc"]))
    sensitivities = list(
        zip(repeat("20240709"), ["mol_import", "nuc_cost", "nuc_extension", "storage_cost", "pure_opt"]))

    runs = scenarios + sensitivities

    Parallel(n_jobs=4)(delayed(compute_scenario_data)(r, s) for r, s in runs)


if __name__ == "__main__":
    main()
