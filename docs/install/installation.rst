.. _installation:

*********************************
CWAS-Plus installation
*********************************


Requirements
###############

For CWAS-Plus to run, the users need to (1) install conda and Ensembl Variant Effect Predictor (VEP) and (2) download a few databases for annotation.

**1. Required installations**

- **Conda**: CWAS-Plus install packages using conda environment. Please install conda or mamba.
- **VEP**: VEP is used for variant annotation. Please refer to the `reference <https://ensembl.org/info/docs/tools/vep/script/vep_download.html>`_ and install VEP.


**2. VEP resources**

- **Missense database**: Missense database containing scores to classify damaging missense variants are required. In CWAS-Plus, MPC database is available. If users want to use other database, they can download their own database and customize the ``configuration.txt`` to use it.

  .. code-block:: solidity
    
    git clone https://github.com/joonan-lab/cwas-dataset.git
    
    
- **loftee plugin**: loftee will be used to classify protein-truncating variants. Be aware of which branch you are cloning.

  .. code-block:: solidity
    
    cd $HOME/.vep/Plugins
    git clone -b grch38 https://github.com/konradjk/loftee.git
    
    
- **gerp bigwig**: This file will be used for loftee plugin.

  .. code-block:: solidity

    wget https://personal.broadinstitute.org/konradk/loftee_data/GRCh38/gerp_conservation_scores.homo_sapiens.GRCh38.bw

    
- **Human ancestor fasta**: This file will be used for loftee plugin.

  .. code-block:: solidity
    
    wget https://personal.broadinstitute.org/konradk/loftee_data/GRCh38/human_ancestor.fa.gz
    wget https://personal.broadinstitute.org/konradk/loftee_data/GRCh38/human_ancestor.fa.gz.fai
    wget https://personal.broadinstitute.org/konradk/loftee_data/GRCh38/human_ancestor.fa.gz.gzi


    
- **Conservation file**: This file will be used for loftee plugin.

  .. code-block:: solidity
    
    wget https://personal.broadinstitute.org/konradk/loftee_data/GRCh38/loftee.sql.gz
    gunzip -k loftee.sql.gz




Install CWAS-Plus
####################


To install CWAS-Plus, git clone the repository from github or download the package.
We recommend using a conda environment with python installed to avoid global installation. Use ``-f`` option when overwriting the package.


.. code-block:: solidity
    
    conda create -n cwas python=3.10
    conda activate cwas
    git clone https://github.com/joonan-lab/cwas.git
    cd cwas
    pip install .


To start CWAS-Plus, type the command below. This will create a workspace (``.cwas``) for CWAS-Plus in home directory. You can specify the directory that will be used as a working directory. As a default, ``$HOME/.cwas`` will be set. If you have a pre-installed VEP, this process will find it automatically and type it to the configuration file.

- -w: Path to the CWAS working directory. All default CWAS processes will save their output here if no specific output directory is given. By default, the directory is set to ``$HOME/.cwas``.

.. code-block:: solidity

    cwas start -w .cwas_wd



