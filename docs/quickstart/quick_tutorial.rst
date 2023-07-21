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
    cwas start


  Download example input.

  .. code-block:: solidity

    cd $HOME
    git clone https://github.com/joonan-lab/cwas-input-example.git


2. :ref:`Configuration <configuration>`
###########################################

  Inside the CWAS working directory, there is a configuration file (``configuration.txt``). Fill in the file with paths of the required tools and data.

  Download required resources (annotation datasets and VEP resources) for CWAS-Plus. This might take a few hours depending on the speed of the network.

  .. code-block:: solidity

    cd $HOME
    git clone https://github.com/joonan-lab/cwas-dataset.git
    cd cwas-dataset
    tar -zxvf functional_annotations.tar.gz # Decompress bed files
    mv functional_annotations/* . # Move bed files to the parent directory
    sh download_vep_resources.sh

  Copy the ``configuration.txt`` in the ``cwas-dataset`` to the CWAS-Plus working directory (by default, ``$HOME/.cwas``).

  After copying, modify the path of the *VEP*, *ANNOTATION_DATA_DIR* and *VEP_CACHE_DIR* to the exact path from the user's environment.

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
   - -p, --num_proc: Number of worker processes that will be used for the annotation process. By default, 1.
   - -o_dir, --output_directory: Path to the directory where the output files will be saved. By default, outputs will be saved at ``$CWAS_WORKSPACE``.


  Annotate the input VCF file with VEP and bed custom annotation algorithm.

  .. code-block:: solidity

    cwas annotation -v $HOME/cwas-input-example/de_novo_variants.vcf -o_dir $HOME/cwas_output -p 8

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

    cwas categorization -i $HOME/cwas_output/de_novo_variants.annotated.vcf -o_dir $HOME/cwas_output -p 8 -m variant

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
        
        cwas binomial_test -i $HOME/cwas_output/de_novo_variants.categorization_result.txt.gz -o_dir $HOME/cwas_output -s $HOME/cwas-input-example/samples.txt -a $HOME/cwas-input-example/adj_factors.txt

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
        
        cwas permutation_test -i $HOME/cwas_output/de_novo_variants.categorization_result.txt.gz -o_dir $HOME/cwas_output -s $HOME/cwas-input-example/samples.txt -a $HOME/cwas-input-example/adj_factors.txt -n 10000 -p 8 -b


7.  :ref:`Calculate the number of effective tests <effnumtest>`
#################################################################

  From correlation matrix, compute eigen values and vectors. Based on these outputs, users can calculate the number of effective tests.

  The parameters of the command are as below:

    - -i, --input_file: Path to the concatenated z-scores.
    - -if, --input_format: Specify the format of the input file. Available options are ``corr`` or ``inter``. By default, ``corr`` will be used. Each format refers to the following:

      - corr: A matrix with correlation values between categories.
      - inter: A matrix with intersected number of variants (or samples) between categories.

    - -o_dir, --output_directory: Path to the directory where the output files will be saved. By default, outputs will be saved at ``$CWAS_WORKSPACE``.
    - -n, --num_sim: Number of eigen values to use in calculating the number of effective tests. The maximum number is equivalent to the number of categories. By default, 10000.
    - -s, --sample_info: Path to the txt file containing the sample information for each sample. This file must have three columns (``SAMPLE``, ``FAMILY``, ``PHENOTYPE``) with the exact name. Required only when input format is set to ``inter``. By default, None.
    - -t, --tag: Tag used for the name of the output files. By default, None.
    - -c, --category_set: Path to a text file containing categories for eigen decomposition. If not specified, all of the categories in the z-score file will be used. This file must contain ``Category`` column with the name of categories to be used.
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


  Create a category set with categories with more than 7 variants.

  .. code-block:: solidity
    
    zcat $HOME/cwas_output/de_novo_variants.category_counts.txt.gz | head -1 > $HOME/cwas_output/subset_categories.v2.txt
    zcat $HOME/cwas_output/de_novo_variants.category_counts.txt.gz | awk '$2 > 7' >> $HOME/cwas_output/subset_categories.v2.txt

  Now run the below command.

  .. code-block:: solidity
    
    cwas effective_num_test -i $HOME/cwas_output/de_novo_variants.correlation_matrix.pkl -o_dir $HOME/cwas_output -ef -if corr -n 10000 -c $HOME/cwas_output/subset_categories.v2.txt



8.  :ref:`Risk score analysis <riskscore>`
##############################################
 

  Identify the overrepresented domains associated to the phenotype.

  The parameters of the command are as below:
  
  - -i, --input_file: Path to the categorized txt file, resulted from categorization process. This file could be gzipped or not.
  - -o_dir, --output_directory: Path to the directory where the output files will be saved. By default, outputs will be saved at ``$CWAS_WORKSPACE``.
  - -s, --sample_info: Path to the txt file containing the sample information for each sample. This file must have three columns (``SAMPLE``, ``FAMILY``, ``PHENOTYPE``) with the exact name.
  - -a, --adjustment_factor: Path to the txt file containing the adjust factors for each sample. This is optional. With this option, CWAS-Plus multiplies the number of variants (or carriers, in -u option) with the adjust factor per sample.
  - -c, --category_set: Path to a text file containing categories for training. If not specified, all of the categories categorization file will be used. This file must contain ``Category`` column with the name of categories to be used.
  - -t, --tag: Tag used for the name of the output files. By default, None.
  - -u, --use_n_carrier: Enables the use of the number of samples with variants in each category for burden test instead of the number of variants. With this option, CWAS-Plus counts the number of samples that carry at least one variant of each category.
  - -thr, --threshold: The number of variants in controls (or the number of control carriers) used to select rare categories. For example, if set to 3, categories with less than 3 variants in controls will be used for training. By default, 3.
  - -tf, --train_set_fraction: The fraction of the training set. For example, if set to 0.7, 70% of the samples will be used as training set and 30% will be used as test set. By default, 0.7.
  - -n_reg, --num_regression: Number of regression trials to calculate a mean of R squares. By default, 10.
  - -f, --fold: Number of folds for cross-validation.
  - -l, --logistic:  Make a logistic model with L1 penalty (Lasso model). By default, False.
  - -n, --n_permute: The number of permutations used to calculate the p-value. By default, 1,000.
  - --predict_only: If set, only predict the risk score and skip the permutation process. By default, False.
  - -p, --num_proc: Number of worker processes that will be used for the permutation process. By default, 1.


  Create a category set with noncoding categories.

  .. code-block:: solidity
    
    zcat $HOME/cwas_output/de_novo_variants.category_info.txt.gz | head -1 > $HOME/cwas_output/subset_categories.txt
    zcat $HOME/cwas_output/de_novo_variants.category_info.txt.gz | awk '$12 == 1 && $6 == "EncodeTFBS"' >> $HOME/cwas_output/subset_categories.txt

  Now run the below command.

  .. code-block:: solidity
    
    cwas risk_score -i $HOME/cwas_output/de_novo_variants.categorization_result.txt.gz \
    -o_dir $HOME/cwas_output \
    -s $HOME/cwas-input-example/samples.txt \
    -a $HOME/cwas-input-example/adj_factors.txt \
    -c $HOME/cwas_output/subset_categories.txt \
    -thr 3 \
    -tf 0.7 \
    -n_reg 10 \
    -f 5 \
    -n 1000 \
    -p 8



9.  :ref:`Burden shift analysis <burdenshift>`
################################################

  Identify the overrepresented domains associated to the phenotype.

  The parameters of the command are as below:

  - -i, --input_file: Path to the input file which is the result of binomial burden test (\*.burden_test.txt.gz).
  - -b, --burden_res: Path to the result of burden shift from permutation test (\*.binom_pvals.txt.gz).
  - -o_dir, --output_directory: Path to the directory where the output files will be saved. By default, outputs will be saved at ``$CWAS_WORKSPACE``.
  - -c, --category_info: Path to the category information file from binomial burden test (\*.category_info.txt.gz).
  - -c_count, --cat_count: Path of the categories counts file from binomial burden test (\*.category_counts.txt.gz).
  - -t, --tag: Tag used for the name of the output files. By default, None.
  - -c_cutoff, --count_cutoff: The number of cutoff for category counts. It must be positive value. By default, 7.
  - --pval: P-value threshold. By default, 0.05.

  .. code-block:: solidity
    
    cwas burden_shift -i $HOME/cwas_output/de_novo_variants.burden_test.txt.gz \
    -b $HOME/cwas_output/de_novo_variants.binom_pvals.txt.gz \
    -o_dir $HOME/cwas_output \
    -c $HOME/cwas_output/de_novo_variants.category_info.txt.gz \
    -c_count $HOME/cwas_output/de_novo_variants.category_counts.txt.gz \
    -c_cutoff 7 \
    --pval 0.05




10.   :ref:`DAWN analysis <dawn>`
####################################

  Investigate the relationship between categories and identify the specific type of categories clustered within the network.

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
  - -c, --category_set: Path to a text file containing categories for training. If not specified, all of the categories categorization file will be used. This file must contain ``Category`` column with the name of categories to be used.
  - -c_count, --cat_count: Path of the categories counts file from burden test
  - -CT, --count_threshold: The treshold of variant (or sample) counts. The least amount of variants a category should have.
  - -CR, --corr_threshold: The threshold of correlation values between clusters. Computed by the mean value of correlation values of categories within a cluster.
  - -S, --size_threshold: The threshold of the number of categories per cluster. The least amount of categories a cluster should have.
  - -p, --num_proc: Number of worker processes that will be used for the DAWN analysis. By default, 1.


  .. code-block:: solidity
  
      cwas dawn -i_dir $HOME/cwas_output \
      -o_dir $HOME/cwas_output \
      -r 2,500 \
      -s 123 \
      -t test \
      -c $HOME/cwas_output/subset_categories.txt \
      -c_count $HOME/cwas_output/de_novo_variants.category_counts.txt.gz \
      -CT 2 \
      -CR 0.7 \
      -S 20 \
      -p 8


