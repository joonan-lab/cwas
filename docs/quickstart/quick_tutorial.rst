*********************************
Quick tutorial for CWAS-Plus
*********************************

This is a quick tutorial for CWAS-Plus. Specific descriptions of arguments are described in the page of each step.



1. :ref:`Install CWAS-Plus <installation>`
###########################################

  The users can install CWAS-Plus through Github. Create a conda environment and install the package.
  
  ``cwas start`` command creates a working directory (``-w``) along with a configuration file.

  .. code-block:: solidity
    
    git clone https://github.com/joonan-lab/cwas.git
    cd cwas
    conda env create -f environment.yml -n cwas
    conda activate cwas
    python setup.py install
    cwas start -w .cwas_wd

2. :ref:`Configuration <configuration>`
###########################################

  Inside the CWAS working directory, there is a configuration file (``configuration.txt``). Fill in the file with paths of the required tools and data.

  After filling the configuration file, ``cwas configuration`` command will create symlinks of annotation datasets into the working directory and fill the ``.cwas_env`` file in the home directory for storing environmental variables.

  .. code-block:: solidity

    cwas configuration

3. :ref:`Prepare annotation datasets <data-prep-label>`
###########################################

  The parameters of the command are as below:

   - p: The number of processors.

  Gather and merge functional annotations and scores into a single bed file.

  .. code-block:: solidity

    cwas preparation -p 8

4. :ref:`Annotation <annotation>`
###########################################

  The parameters of the command are as below:

   - -v, --vcf_file: Path to the input vcf file. This file could be gzipped or not.
   - -n, --num_cores: Number of worker processes that will be used for the annotation process. By default, 1.
   - -o_dir, --output_directory: Path to the directory where the output files will be saved. By default, outputs will be saved at ``$CWAS_WORKSPACE``.


  Annotate the input VCF file with VEP and bed custom annotation algorithm.

  .. code-block:: solidity

    cwas annotation -v INPUT.vcf -o_dir OUTPUT_DIR -n 8

5. :ref:`Categorization <categorization>`
###########################################

  The parameters of the command are as below:

   - -i, --input_file: Path to the annotated VCF, resulted from annotation process. This file contains a specific pattern of ``.annotated.vcf`` in the file name. This file could be gzipped or not.
   - -o_dir, --output_directory: Path to the directory where the output files will be saved. By default, outputs will be saved at ``$CWAS_WORKSPACE``.
   - -p, --num_proc: Number of worker processes that will be used for the categorization process. To prevent crashes caused by insufficient RAM when processing large input VCF files (e.g., over 10 million variants) using multiple cores, using small number of cores and monitoring the memory usage are recommended. By default, 1.
   - -m, --matrix: Generate a correlation matrix and a matrix with intersected number of variants (or samples) between every two categories. Available options are ``variant`` or ``sample``. By default, False.

     - variant: Use the intersected number of variants between two categories.
     - sample: Use the intersected number of samples between two categories.


  Categorize variants into groups based on the annotation datasets. A single category is a combination of five domains (i.e., variant type, gene biotype, gene list, functional annotation and functional score). Details are provided in the :ref:`Overview of annotation datasets <overview>`.

  .. code-block:: solidity

    cwas categorization -i INPUT.annotated.vcf -o_dir OUTPUT_DIR -p 8

6. :ref:`Burden test <burdentest>`
######################################

  Calculate the burden of each category by comparing the number of variants per case and control. Two types of tests are used for p-value calculation: binomial test and permutation test.
   
  - Binomial test

    - -i, --input_file: Path to the categorized txt file, resulted from categorization process. This file could be gzipped or not.
    - -o_dir, --output_directory: Path to the directory where the output files will be saved. By default, outputs will be saved at ``$CWAS_WORKSPACE``.
    - -s, --sample_info: Path to the txt file containing the sample information for each sample. This file must have three columns (``SAMPLE``, ``FAMILY``, ``PHENOTYPE``) with the exact name.
    - -a, --adjustment_factor: Path to the txt file containing the adjust factors for each sample. This is optional. With this option, CWAS-Plus multiplies the number of variants (or carriers, in -u option) with the adjust factor per sample.
    - -u, --use_n_carrier: Enables the use of the number of samples with variants in each category for burden test instead of the number of variants. With this option, CWAS-Plus counts the number of samples that carry at least one variant of each category.

     .. code-block:: solidity

        cwas binomial_test -i INPUT.categorization_result.txt.gz -o_dir OUTPUT_DIR -s SAMPLE_LIST.txt -a ADJUST_FACTOR.txt

  - Permutation test

    - -i, --input_file: Path to the categorized txt file, resulted from categorization process. This file could be gzipped or not.
    - -o_dir, --output_directory: Path to the directory where the output files will be saved. By default, outputs will be saved at ``$CWAS_WORKSPACE``.
    - -s, --sample_info: Path to the txt file containing the sample information for each sample. This file must have three columns (``SAMPLE``, ``FAMILY``, ``PHENOTYPE``) with the exact name.
    - -a, --adjustment_factor: Path to the txt file containing the adjust factors for each sample. This is optional. With this option, CWAS-Plus multiplies the number of variants (or carriers, in -u option) with the adjust factor per sample.
    - -n, --num_perm: Number of permutations for label-swapping. By default, 10000.
    - -p, --num_proc: Number of worker processes that will be used for the permutation process. By default, 1.
    - -b, --burden_shift: Generates an output file containing binomial p-values for each label-swapped permutation. By default, False.
    - -rr, --perm_rr: Generates an output file containing relative risks for each label-swapped permutation. By default, False.
    - -u, --use_n_carrier: Enables the use of the number of samples with variants in each category for burden test instead of the number of variants. With this option, CWAS-Plus counts the number of samples that carry at least one variant of each category.

     .. code-block:: solidity

        cwas permutation_test -i INPUT.categorization_result.txt.gz -o_dir OUTPUT_DIR -s SAMPLE_LIST.txt -a ADJUST_FACTOR.txt -n 10000 -p 8 -b


7. :ref:`Caculate the correlation matrix <categorization>`
#############################################################

  Caculate the correlation matrix from the intersected number of variants (or samples) between every two categories.

  The parameters of the command are as below:

    - -i, --input_file: Path to the concatenated z-scores.
    - -if, --input_format: Specify the format of the input file. Available options are ``corr`` or ``inter``. By default, ``corr`` will be used. Each format refers to the following:

      - corr: A matrix with correlation values between categories.
      - inter: A matrix with intersected number of variants (or samples) between categories.

    - -o_dir, --output_directory: Path to the directory where the output files will be saved. By default, outputs will be saved at ``$CWAS_WORKSPACE``.
    - -n, --num_sim: Number of eigen values to use in calculating the number of effective tests. The maximum number is equivalent to the number of categories. By default, 10000.
    - -s, --sample_info: Path to the txt file containing the sample information for each sample. This file must have three columns (``SAMPLE``, ``FAMILY``, ``PHENOTYPE``) with the exact name. Required only when input format is set to ``inter``. By default, None.
    - -t, --tag: Tag used for the name of the output files. By default, None.
    - -c, --category_set_path: Path to a text file containing categories for eigen decomposition. If not specified, all of the categories in the z-score file will be used. This file must contain ``Category`` column with the name of categories to be used.
    - -ef, --eff_num_test: Calculate the effective number of tests. By default, False.

    +-------------------------------------------------------+
    |Category                                               |
    +=======================================================+
    |All_CHD8Common_All_IntergenicRegion_EarlyCREMicro      |
    +-------------------------------------------------------+
    |All_CHD8Common_phastCons46way_PromoterRegion_EarlyCREL4|
    +-------------------------------------------------------+
    |All_DDD_All_PromoterRegion_EarlyCREOligo               |
    +-------------------------------------------------------+

  .. code-block:: solidity

    cwas categorization -i INPUT.annotated.vcf -o_dir OUTPUT_DIR -p 8 -m variant


8.  :ref:`Calculate the number of effective tests <effnumtest>`
#################################################################

  From correlation matrix, compute eigen values and vectors. Based on these outputs, users can calculate the number of effective tests.

  The parameters of the command are as below:
  
  - -i, --input_file: Path to the categorized txt file, resulted from categorization process. This file could be gzipped or not.
  - -o_dir, --output_directory: Path to the directory where the output files will be saved. By default, outputs will be saved at ``$CWAS_WORKSPACE``.
  - -s, --sample_info: Path to the txt file containing the sample information for each sample. This file must have three columns (``SAMPLE``, ``FAMILY``, ``PHENOTYPE``) with the exact name.
  - -a, --adjustment_factor: Path to the txt file containing the adjust factors for each sample. This is optional. With this option, CWAS-Plus multiplies the number of variants (or carriers, in -u option) with the adjust factor per sample.
  - -c, --category_set_path: Path to a text file containing categories for training. If not specified, all of the categories categorization file will be used. This file must contain ``Category`` column with the name of categories to be used.
  - -t, --tag: Tag used for the name of the output files. By default, None.
  - -u, --use_n_carrier: Enables the use of the number of samples with variants in each category for burden test instead of the number of variants. With this option, CWAS-Plus counts the number of samples that carry at least one variant of each category.
  - -thr, --threshold: The number of variants in controls (or the number of control carriers) used to select rare categories. For example, if set to 3, categories with less than 3 variants in controls will be used for training. By default, 3.
  - -tf, --train_set_fraction: The fraction of the training set. For example, if set to 0.7, 70% of the samples will be used as training set and 30% will be used as test set. By default, 0.7.
  - -n_reg, --num_regression: Number of regression trials to calculate a mean of R squares. By default, 10.
  - -f, --fold: Number of folds for cross-validation.
  - -l, --logistic: (hold) Make a logistic model with L1 penalty. By default, False.
  - -n, --n_permute: The number of permutations used to calculate the p-value. By default, 1,000.
  - --predict_only: If set, only predict the risk score and skip the permutation process. By default, False.
  - -p, --num_proc: Number of worker processes that will be used for the permutation process. By default, 1.


  .. code-block:: solidity

    cwas effective_num_test -i INPUT.correlation_matrix.pkl -o_dir OUTPUT_DIR -t test -c CATEGORY_SET.txt -ef


9.  :ref:`Risk score analysis <riskscore>`
##############################################

  Identify the best predictor of the phenotype by training Lasso regression model with the number of variants within each category across samples.

  .. code-block:: solidity

    cwas risk_score -i INPUT.categorization_result.txt.gz \
    -o_dir OUTPUT_DIR \
    -s SAMPLE_LIST.txt \
    -a ADJUST_FACTOR.txt \
    -c CATEGORY_SET.txt \
    -thr 3 \
    -tf 0.7 \
    -n_reg 10 \
    -f 5 \
    -n 1000 \
    -p 8


10.  :ref:`Burden shift analysis <riskscore>`
##############################################

  Identify the overrepresented domains associated to the phenotype.

  The parameters of the command are as below:

  - -i_dir, --input_directory: Path to the directory where the input files are stored. This directory must include three required files.

    - Eigen vector file: This is the output file from :ref:`calculation of effective number of tests <effnumtest>`. The file name must have pattern ``*eig_vecs*.txt.gz``.
    - Category correlation matrix file: This is the output file from :ref:`categorization <categorization>`. The file name must have pattern ``*correlation_matrix*.pkl``.
    - Permutation test file: This is the output file from :ref:`burden test <permtest>`. The file name must have pattern ``*permutation_test*.txt.gz``.

  - -o_dir, --output_directory: Path to the directory where the output files will be saved. By default, outputs will be saved at ``$CWAS_WORKSPACE``.
  - -r, --range: Range (i.e., (start,end)) to find optimal K for k-means clustering. It must contain two integers that are comma-separated. The first integer refers to the start number and must be above 1. The second integer refers to the end.
  - -k, --k_val: K for K-means clustering. With this argument, users can determine K manually. ``-r`` and ``-k`` arguments are mutually exclusive. If ``-k`` is given, ``-r`` will be ignored.
  - -s, --seed: Seed value for t-SNE. Same seed will generate same results for the same inputs.
  - -t, --tag: Tag used for the name of the output files. By default, None.
  - -c, --category_set_path: Path to a text file containing categories for training. If not specified, all of the categories categorization file will be used. This file must contain ``Category`` column with the name of categories to be used.
  - -c_count, --cat_count
  - -CT, --count_threshold: The treshold of variant (or sample) counts. The least amount of variants a category should have.
  - -CR, --corr_threshold: The threshold of correlation values between clusters. Computed by the mean value of correlation values of categories within a cluster.
  - -S, --size_threshold: The threshold of the number of categories per cluster. The least amount of categories a cluster should have.
  - -p, --num_proc: Number of worker processes that will be used for the DAWN analysis. By default, 1.


  .. code-block:: solidity


11.  :ref:`DAWN analysis <dawn>`
####################################

  Investigate the relationship between categories and identify the specific type of categories clustered within the network.

  .. code-block:: solidity
  
      cwas dawn -i_dir INPUT_DIR \
      -o_dir OUTPUT_DIR \
      -r 2,500 \
      -s 123 \
      -t test \
      -c CATEGORY_SET.txt \
      -c_count CATEGORY_COUNTS.txt \
      -CT 2 \
      -CR 0.7 \
      -S 20 \
      -p 8


