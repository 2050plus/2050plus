..
  SPDX-FileCopyrightText: 2019-2023 The PyPSA-Eur Authors

  SPDX-License-Identifier: CC-BY-4.0

.. _veka_starter:

##########################################
Starter kit
##########################################

To start using the model and its dependencies, multiple scripts and tool should be used.

1. Get the model
    * Clone the git repository from `GitHub <https://github.com/2050plus/2050plus>`_.
    * Install your environment using `envs` folder (latest data were produced with `envs/environment.aws_r6a.12xlarge.yaml`). A detailed installation guide is available in the `PyPSA-Eur documentation <https://pypsa-eur.readthedocs.io/en/latest/installation.html>`_.
    * Optimisation has been done using Gurobi (commercial solver). HiGHS is a good open source alternative. Other alternatives are listed in the `PyPSA-Eur documentation <https://pypsa-eur.readthedocs.io/en/latest/installation.html#install-a-solver>`_.

2. Run the model
    * Use snakemake using the following command

    .. code:: bash

         snakemake -call all --configfile config/config.veka.yaml -n

    * Model can be ran locally or on a remote server to improve performances (should work on a AWS EC2 r6a.12xlarge).
    * Edit the configurations files to improve the scenarios : `config/config.veka.yaml` and `config/scenarios.veka.yaml`.

3. Extract data from the PyPSA-Eur outputs
    * Copy desired files from `results` and `resources` folders to `analysis/{run}` folder. If you use a remote server, use SFTP to download data from the server.
    * Configure the correct run versions in `scripts/graph_extraction_main.py` for scenarios, sensitivities and reference runs.
    * Run the data extraction pipeline (locally).
    * The output data (`analysis/{run}/graph_extraction_st`) are created in `analysis` folder and should be copied in the `app/assets/data` folder if you want to add them to the streamlit.
4. Ship the data in streamlit
    * Streamlit platform is free to use when the data are available on Github publicly. Currently our `website <https://climact-veka-2050plus.streamlit.app/>`_ is linked to `main-v10` branch on the `Github of Climact <https://github.com/Climact/2050plus-climact/tree/main-v10>`_.
    * To use the new data, reconfigure `app/st_common.py` with the new runs paths (change `scenario_dict` values).
    * Merge the updated data with the GitHub branch `main-v10` to ship the data online.
    * Streamlit will automatically update the webpage. If you want to be sure that cached data are updated, you can reboot the app on the Streamlit dashboard.
    * The webpage can be debugged locally using:

    .. code:: bash

        streamlit run app/VEKA_2050+.py
    * Since Streamlit can also be run as a container, we have added `app/Dockerfile_demo` as an example.

5. The documentation is shared online through `Readthedocs <https://2050plus-climact.readthedocs.io/en/latest/>`_.
    * Pushing new updates on the `main-v10` branch on GitHub will also automatically update the page.
    * Documentation is stored in the `doc` folder.
    * The documentation can be tested locally using:

    .. code:: bash

        cd doc
        make clean
        make html
        open _build/html/index.html

    * The documentation is also configured to generate a PDF file using:

    .. code:: bash

        cd doc
        make clean
        # Need twice for LaTeX references
        make latexpdf; make latexpdf
        open _build/latex/PyPSA-Eur.pdf
