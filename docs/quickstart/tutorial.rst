=======================
Tutorial for CWAS-Plus
=======================

This is a quick tutorial for CWAS-Plus. Specific description of arguments are described in the page of each step.

1. Install CWAS-Plus

  .. code-block:: solidity
    
    git clone https://github.com/joonan-lab/cwas.git
    cd cwas
    conda env create -f environment.yml -n cwas
    conda activate cwas
    python setup.py install
    cwas start -w .cwas_wd

2. Configuration

  .. code-block:: solidity

    cwas configuration

3. Prepare annotation datasets

  .. code-block:: solidity

    cwas preparation -p 8

4. Annotation

  .. code-block:: solidity

    cwas annotation -v INPUT.vcf -o_dir OUTPUT_DIR -n 8

5. Categorization

  .. code-block:: solidity

    cwas categorization -i INPUT.annotated.vcf -o_dir OUTPUT_DIR -p 8

6. Burden test
   
  - Binomial test

     .. code-block:: solidity

        cwas binomial_test -i INPUT.categorization_result.txt.gz -o_dir OUTPUT_DIR -s SAMPLE_LIST.txt -a ADJUST_FACTOR.txt

  - Permutation test
   
     .. code-block:: solidity

        cwas permutation_test -i INPUT.categorization_result.txt.gz -o_dir OUTPUT_DIR -s SAMPLE_LIST.txt -a ADJUST_FACTOR.txt -n 10000 -p 8 -b



