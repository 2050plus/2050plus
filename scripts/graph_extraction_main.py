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
from pathlib import Path

from scripts.graph_extraction_extract import extract_data
from scripts.graph_extraction_load_st import load_data_st
from scripts.graph_extraction_load_xl import load_data_xl
from scripts.graph_extraction_transform import transform_data
from scripts.graph_extraction_utils import load_config

logger = logging.getLogger(__name__)


def main():
    logger.info("Start processing")

    config_file = "config/config.veka.yaml"
    analysis_path = Path("analysis", "20240616")
    scenario = "central"
    dir_export = "graph_data"

    # Configuration
    config = load_config(config_file, analysis_path, dir_export, scenario=scenario)

    config["eu27_countries"] = ["AT", "BE", "BG", "CY", "CZ", "DE", "DK", "EE", "GR", "ES", "FI", "FR", "HR", "HU",
                                "IE", "IT", "LT", "LU", "LV", "MT", "NL", "PL", "PT", "RO", "SE", "SI", "SK"]
    config["eu27_countries"] = list(set(config["eu27_countries"]).intersection(set(config["countries"])))
    # global variables for which to do work
    config["countries"] = {"tot": None, "be": ["BE"], "fl": ["FL"], "eu27": config["eu27_countries"]}
    config["imp_exp_carriers"] = ["elec", "gas", "H2"]
    config["carriers_to_plot"] = ["electricity", "heat"]

    # Extract data
    n, n_ext = extract_data(config, )

    # Transform data
    transform_data(config, n, n_ext)

    # Load data
    load_data_xl(config)
    load_data_st(config)

    logger.info("Done")


if __name__ == "__main__":
    main()
