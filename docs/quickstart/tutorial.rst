=======================
Tutorial for CWAS-Plus
=======================

This is a quick tutorial for CWAS-Plus. Specific description of arguments are described in the page of each step.



1. :ref:`Install CWAS-Plus <installation>`

  The users can install CWAS-Plus through Github. Create a conda environment and install the package.
  
  ``cwas start`` command creates a working directory along with a configuration file.

  .. code-block:: solidity
    
    git clone https://github.com/joonan-lab/cwas.git
    cd cwas
    conda env create -f environment.yml -n cwas
    conda activate cwas
    python setup.py install
    cwas start -w .cwas_wd

2. :ref:`Configuration <configuration>`

  Inside the CWAS working directory, there is a configuration file (``configuration.txt``). Fill in the file with paths of the required tools and data.

  After filling the configuration file, ``cwas configuration`` command will create symlinks of annotation datasets into the working directory and fill the ``.cwas_env`` file in the home directory for storing environmental variables.

  .. code-block:: solidity

    cwas configuration

3. :ref:`Prepare annotation datasets <data-prep-label>`

  Gather and merge functional annotations and scores into a single bed file.

  .. code-block:: solidity

    cwas preparation -p 8

4. :ref:`Annotation <annotation>`

  Annotate the input VCF file with VEP and bed custom annotation algorithm.

  .. code-block:: solidity

    cwas annotation -v INPUT.vcf -o_dir OUTPUT_DIR -n 8

5. :ref:`Categorization <categorization>`

  Categorize variants into groups based on the annotation datasets. A single category is a combination of five domains (i.e., variant type, gene biotype, gene list, functional annotation and functional score). Details are provided in the :ref:`Overview of annotation datasets <overview>`.

  .. code-block:: solidity

    cwas categorization -i INPUT.annotated.vcf -o_dir OUTPUT_DIR -p 8

6. :ref:`Burden test <burdentest>`

  Calculate the burden of each category by comparing the number of variants per case and control. Two types of tests are used for p-value calculation: binomial test and permutation test.
   
  - Binomial test

     .. code-block:: solidity

        cwas binomial_test -i INPUT.categorization_result.txt.gz -o_dir OUTPUT_DIR -s SAMPLE_LIST.txt -a ADJUST_FACTOR.txt

  - Permutation test
   
     .. code-block:: solidity

        cwas permutation_test -i INPUT.categorization_result.txt.gz -o_dir OUTPUT_DIR -s SAMPLE_LIST.txt -a ADJUST_FACTOR.txt -n 10000 -p 8 -b



