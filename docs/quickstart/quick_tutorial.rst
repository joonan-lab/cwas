*********************************
Quick tutorial for CWAS-Plus
*********************************

This is a quick tutorial for CWAS-Plus. Specific descriptions of arguments are described in the page of each step.



1. :ref:`Install CWAS-Plus <installation>`
###########################################

  Users can install CWAS-Plus through pip or github. We recommend installing under conda environment to avoid global installation.
  
  ``cwas start`` command creates a working directory (``-w``) along with a configuration file.

  - Github

  .. code-block:: solidity
    
    conda create -n cwas python=3.10 r-base=4.2.2
    conda activate cwas
    git clone https://github.com/joonan-lab/cwas.git
    pip install cwas
    cwas start

  - pip

  .. code-block:: solidity
    
    conda create -n cwas python=3.10 r-base=4.2.2
    conda activate cwas
    pip install cwas
    cwas start

  The installation of R package **glmnet** is also required for risk score analysis.

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
    git lfs pull
    tar -zxvf functional_annotations.tar.gz # Decompress bed files
    mv functional_annotations/* . # Move bed files to the parent directory
    sh download_vep_resources.sh

  Copy the ``configuration.txt`` in the ``cwas-dataset`` to the CWAS-Plus working directory (by default, ``$HOME/.cwas``).

  .. code-block:: solidity
    
    cp $HOME/cwas-dataset/configuration.txt $HOME/.cwas/

  After copying, modify the path of the *VEP*, *ANNOTATION_DATA_DIR* and *VEP_CACHE_DIR* to the exact path from the user's environment.

  When preparing the ANNOTATION_KEY_CONFIG yaml file, please avoid using underscores ('_') in the annotation name. Underscores are used for distinguishing different domains within a single category.

  For example, check below.

  .. code-block:: solidity

    functional_score:
      bed1.bed.gz: annot1
      bed2.bed.gz: annot2
    functional_annotation:
      bed3.bed.gz: annot_3 # Do not use underscores like this. Users can use 'annot3' instead.
      bed4.bed.gz: annot4

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

   - -v, --vcf_file: Path to the input vcf file. This file could be bgzipped or not.
   - -p, --num_proc: Number of worker processes that will be used for the annotation process. By default, 1.
   - -o_dir, --output_directory: Path to the directory where the output files will be saved. By default, outputs will be saved at ``$CWAS_WORKSPACE``.


  Annotate the input VCF file with VEP and bed custom annotation algorithm.

  .. code-block:: solidity

    cwas annotation -v $HOME/cwas-input-example/de_novo_variants.vcf -o_dir $HOME/cwas_output -p 8

  Below are the output files generated.

  .. code-block:: solidity

    $HOME/cwas_output
    ...
    ├── de_novo_variants.vep.vcf.gz
    ├── de_novo_variants.vep.vcf.gz.tbi
    ├── de_novo_variants.annotated.vcf.gz
    ...

5. :ref:`Categorization <categorization>`
###########################################

  The parameters of the command are as below:

   - -i, --input_file: Path to the annotated VCF, resulted from annotation process. This file contains a specific pattern of ``.annotated.vcf`` in the file name. This file could be gzipped or not.
   - -o_dir, --output_directory: Path to the directory where the output files will be saved. By default, outputs will be saved at ``$CWAS_WORKSPACE``.
   - -p, --num_proc: Number of worker processes that will be used for the categorization process. To prevent crashes caused by insufficient RAM when processing large input VCF files (e.g., over 10 million variants) using multiple cores, using small number of cores and monitoring the memory usage are recommended. By default, 1.

  Categorize variants into groups based on the annotation datasets. A single category is a combination of five domains (i.e., variant type, gene biotype, gene list, functional annotation and functional score). Details are provided in the :ref:`Overview of annotation datasets <overview>`.

  .. code-block:: solidity

    cwas categorization -i $HOME/cwas_output/de_novo_variants.annotated.vcf.gz -o_dir $HOME/cwas_output -p 8

  .. code-block:: solidity

    $HOME/cwas_output
    ...
    ├── de_novo_variants.categorization_result.zarr
    ├── de_novo_variants.intersection_matrix.zarr
    ├── de_novo_variants.correlation_matrix.zarr
    ...


6. :ref:`Burden test <burdentest>`
######################################

  Calculate the burden of each category by comparing the number of variants per case and control. Two types of tests are used for p-value calculation: binomial test and permutation test.
   
  - Binomial test

    - -i, --input_file: Path to the categorized zarr directory, resulted from categorization process.
    - -o_dir, --output_directory: Path to the directory where the output files will be saved. By default, outputs will be saved at ``$CWAS_WORKSPACE``.
    - -s, --sample_info: Path to the txt file containing the sample information for each sample. This file must have three columns (``SAMPLE``, ``FAMILY``, ``PHENOTYPE``) with the exact name.
    - -a, --adjustment_factor: Path to the txt file containing the adjust factors for each sample. This is optional. With this option, CWAS-Plus multiplies the number of variants (or carriers, in -u option) with the adjust factor per sample.
    - -u, --use_n_carrier: Enables the sample-level analysis (use of the number of samples with variants in each category for burden test instead of the number of variants). With this option, CWAS-Plus counts the number of samples that carry at least one variant of each category.

     .. code-block:: solidity
        
        cwas binomial_test -i $HOME/cwas_output/de_novo_variants.categorization_result.zarr -o_dir $HOME/cwas_output -s $HOME/cwas-input-example/samples.txt -a $HOME/cwas-input-example/adj_factors.txt

  Below are the output files generated.

  .. code-block:: solidity

    $HOME/cwas_output
    ...
    ├── de_novo_variants.burden_test.volcano_plot.pdf
    ├── de_novo_variants.burden_test.txt
    ├── de_novo_variants.category_counts.txt
    ├── de_novo_variants.category_info.txt
    ...

  - Permutation test

    - -i, --input_file: Path to the categorized zarr directory, resulted from categorization process.
    - -o_dir, --output_directory: Path to the directory where the output files will be saved. By default, outputs will be saved at ``$CWAS_WORKSPACE``.
    - -s, --sample_info: Path to the txt file containing the sample information for each sample. This file must have three columns (``SAMPLE``, ``FAMILY``, ``PHENOTYPE``) with the exact name.
    - -a, --adjustment_factor: Path to the txt file containing the adjust factors for each sample. This is optional. With this option, CWAS-Plus multiplies the number of variants (or carriers, in -u option) with the adjust factor per sample.
    - -n, --num_perm: Number of permutations for label-swapping. By default, 10000.
    - -p, --num_proc: Number of worker processes that will be used for the permutation process. By default, 1.
    - -b, --burden_shift: Generates an output file containing binomial p-values for each label-swapped permutation. By default, False.
    - -u, --use_n_carrier: Enables the sample-level analysis (use of the number of samples with variants in each category for burden test instead of the number of variants). With this option, CWAS-Plus counts the number of samples that carry at least one variant of each category.

     .. code-block:: solidity
        
        cwas permutation_test -i $HOME/cwas_output/de_novo_variants.categorization_result.zarr -o_dir $HOME/cwas_output -s $HOME/cwas-input-example/samples.txt -a $HOME/cwas-input-example/adj_factors.txt -n 10000 -p 8 -b

  Below are the output files generated.

  .. code-block:: solidity

    $HOME/cwas_output
    ...
    ├── de_novo_variants.permutation_test.txt.gz
    ├── de_novo_variants.binom_pvals.txt.gz
    ...



7.  :ref:`Generate correlation matrix <correlation>`
#################################################################

For (1) calculating study-wide significance threshold and (2) generating DAWN analysis input, correlation values between every two CWAS categories are required.

In this step, users can generate two matrices, (1) a matrix that contains the number of variants (or samples, with --use_carrier option) that intersect between categories, (2) a matrix that contains correlation values between categories. The correlation matrix is computed from the intersected matrix (1). The users can choose one of the matrices for calculating the number of effective tests and DAWN analysis.

The parameters of the command are as below:

- -i, input_file: Path to the categorized zarr directory, resulted from categorization process.
- -v, --annotated_vcf: Path to the annotated VCF, resulted from annotation process. Required for variant-level correlation matrix (`--cm variant`).
- -o_dir, --output_directory: Path to the directory where the output files will be saved. By default, outputs will be saved at ``$CWAS_WORKSPACE``.
- -p, --num_proc: Number of worker processes that will be used for the categorization process. To prevent crashes caused by insufficient RAM when processing large input VCF files (e.g., over 10 million variants) using multiple cores, using small number of cores and monitoring the memory usage are recommended. By default, 1.
- -cm, --corr_matrix: Generate a correlation matrix between every two categories. Available options are ``variant`` or ``sample``. By default, False.

  - variant: Use the intersected number of variants between two categories.
  - sample: Use the intersected number of samples between two categories.

- -im, --intersection_matrix: Generate a matrix with intersected number of variants (or samples with variants) bewteen categories.
- -c_info, --category_info: Path to a text file with category information (`*.category_info.txt`).
- -d, --domain_list: Domain list to filter categories based on GENCODE domain. By default, `all`.

.. code-block:: solidity

    cwas correlation -i INPUT.categorization_result.zarr -v INPUT.annotated.vcf.gz -o_dir OUTPUT_DIR -c_info OUTPUT.category_info.txt -p 8 -cm variant -im

Example run:

.. code-block:: solidity
    
    cwas correlation -i $HOME/cwas_output/de_novo_variants.categorization_result.zarr -v $HOME/cwas_output/de_novo_variants.annotated.vcf.gz -o_dir $HOME/cwas_output -c_info $HOME/cwas_output/de_novo_variants.category_info.txt -p 8 -cm variant -im

Below are the output files generated.

.. code-block:: solidity

  $HOME/cwas_output
  ...
  ├── de_novo_variants.intersection_matrix.zarr
  ├── de_novo_variants.correlation_matrix.zarr
  ...


8.  :ref:`Calculate the number of effective tests <effnumtest>`
#################################################################

  From correlation matrix, compute eigen values and vectors. Based on these outputs, users can calculate the number of effective tests.

  The parameters of the command are as below:

    - -i, --input_file: Path to a matrix of correlation or intersected number of variants between two categories.
    - -if, --input_format: Specify the format of the input file. Available options are ``corr`` or ``inter``. By default, ``corr`` will be used. Each format refers to the following:

      - corr: A matrix with correlation values between categories.
      - inter: A matrix with intersected number of variants (or samples) between categories.

    - -o_dir, --output_directory: Path to the directory where the output files will be saved. By default, outputs will be saved at ``$CWAS_WORKSPACE``.
    - -n, --num_sim: Number of eigen values to use in calculating the number of effective tests. The maximum number is equivalent to the number of categories. By default, 10000.
    - -s, --sample_info: Path to the txt file containing the sample information for each sample. This file must have three columns (``SAMPLE``, ``FAMILY``, ``PHENOTYPE``) with the exact name. Required only when input format is set to ``inter`` or ``-thr`` is not given. By default, None.
    - -c_count, --cat_count: Path of the categories counts file from binomial burden test (\*.category_counts.txt).
    - -t, --tag: Tag used for the name of the output files. By default, None.
    - -c_set, --category_set: Path to a text file containing categories for eigen decomposition. If not specified, all of the categories (surpassing the cutoff) will be used. This file must contain ``Category`` column with the name of categories to be used.

    +-------------------------------------------------------+
    |Category                                               |
    +=======================================================+
    |All_CHD8Common_All_IntergenicRegion_EarlyCREMicro      |
    +-------------------------------------------------------+
    |All_CHD8Common_phastCons46way_PromoterRegion_EarlyCREL4|
    +-------------------------------------------------------+
    |All_DDD_All_PromoterRegion_EarlyCREOligo               |
    +-------------------------------------------------------+

    - -ef, --eff_num_test: Calculate the effective number of tests. For calculation, the users should use all categories (with the number of variants/samples≥cutoff). By default, False.
    - -thr, --threshold: The number of variants (or samples) to filter categories. By default, None.


  Now run the below command.

  .. code-block:: solidity
    
    cwas effective_num_test -i $HOME/cwas_output/de_novo_variants.correlation_matrix.zarr -o_dir $HOME/cwas_output -ef -thr 8 -if corr -n 10000 -c_count $HOME/cwas_output/de_novo_variants.category_counts.txt

    cat $HOME/cwas_output/de_novo_variants.category_info.txt | head -1 > $HOME/cwas_output/subset_categories.txt
    cat $HOME/cwas_output/de_novo_variants.category_info.txt | awk '$12 == 1 && $6 == "EncodeTFBS"' >> $HOME/cwas_output/subset_categories.txt

    cwas effective_num_test -i $HOME/cwas_output/de_novo_variants.correlation_matrix.zarr -o_dir $HOME/cwas_output -thr 8 -if corr -t TFBS -n 10000 -c_set $HOME/cwas_output/subset_categories.txt -c_count $HOME/cwas_output/de_novo_variants.category_counts.txt


  Below are the output files generated.

  .. code-block:: solidity

    $HOME/cwas_output
    ...
    ├── de_novo_variants.neg_lap.pickle
    ├── de_novo_variants.eig_vals.pickle
    ├── de_novo_variants.eig_vecs.txt.gz
    ├── de_novo_variants.neg_lap.TFBS.pickle
    ├── de_novo_variants.eig_vals.TFBS.pickle
    ├── de_novo_variants.eig_vecs.TFBS.txt.gz
    ...

  The number of effective tests will be shown like below.

  .. code-block:: solidity
    
    [RESULT] The number of effective tests is 2596.


9.  :ref:`Risk score analysis <riskscore>`
##############################################
 

  Identify the overrepresented domains associated to the phenotype.

  The parameters of the command are as below:
  
  - -i, --input_file: Path to the categorized zarr directory, resulted from categorization process.
  - -o_dir, --output_directory: Path to the directory where the output files will be saved. By default, outputs will be saved at ``$CWAS_WORKSPACE``.
  - -s, --sample_info: Path to the txt file containing the sample information for each sample. This file must have three columns (``SAMPLE``, ``FAMILY``, ``PHENOTYPE``) with the exact name.
  - -a, --adjustment_factor: Path to the txt file containing the adjust factors for each sample. This is optional. With this option, CWAS-Plus multiplies the number of variants (or carriers, in -u option) with the adjust factor per sample.
  - -c_info, --category_info: Path to a text file category information (`*.category_info.txt`).
  - -d, --domain_list: Domain list to filter categories based on GENCODE domain. If 'run_all' is given, all available options will be tested. Available options are `run_all,all,coding,noncoding,ptv,missense,damaging_missense,promoter,noncoding_wo_promoter,intron,intergenic,utr,lincRNA`. By default, all.
  - -t, --tag: Tag used for the name of the output files. By default, None.
  - --do_each_one: Use each annotation from functional annotation to calculate risk score. By default, False.
  - --leave_one_out: Calculate risk score while excluding one annotation from functional annotation. This option is not used when the `--do_each_one` flag is enabled. By default, False.
  - -u, --use_n_carrier: Enables the sample-level analysis (the use of the number of samples with variants in each category for burden test instead of the number of variants). With this option, CWAS-Plus counts the number of samples that carry at least one variant of each category.
  - -thr, --threshold: The number of variants in controls (or the number of control carriers) used to select rare categories. For example, if set to 3, categories with less than 3 variants in controls will be used for training. By default, 3.
  - -tf, --train_set_fraction: The fraction of the training set. For example, if set to 0.7, 70% of the samples will be used as training set and 30% will be used as test set. By default, 0.7.
  - -n_reg, --num_regression: Number of regression trials to calculate a mean of R squares. By default, 10.
  - -f, --fold: Number of folds for cross-validation.
  - -n, --n_permute: The number of permutations used to calculate the p-value. By default, 1,000.
  - --predict_only: If set, only predict the risk score and skip the permutation process. By default, False.
  - -S, --seed: Seed of random state. By default, 42.
  - -p, --num_proc: Number of worker processes that will be used for the permutation process. By default, 1.


  Now run the below command. The below command calculates risk scores for noncoding domain categories.

  .. code-block:: solidity
    
    cwas risk_score -i $HOME/cwas_output/de_novo_variants.categorization_result.zarr \
    -o_dir $HOME/cwas_output \
    -s $HOME/cwas-input-example/samples.txt \
    -a $HOME/cwas-input-example/adj_factors.txt \
    -c_info $HOME/cwas_output/de_novo_variants.category_info.txt \
    -d noncoding \
    -thr 3 \
    -tf 0.7 \
    -n_reg 10 \
    -f 5 \
    -n 1000 \
    -p 8

  Below are the output files generated.

  .. code-block:: solidity

    $HOME/cwas_output
    ...
    ├── de_novo_variants.lasso_coef_thres_3.txt
    ├── de_novo_variants.lasso_null_models_thres_3.txt
    ├── de_novo_variants.lasso_results_thres_3.txt
    ├── de_novo_variants.lasso_histogram_thres_3.pdf
    ...


10.   :ref:`Burden shift analysis <burdenshift>`
################################################

  Identify the overrepresented domains associated to the phenotype.

  The parameters of the command are as below:

  - -i, --input_file: Path to the input file which is the result of binomial burden test (\*.burden_test.txt).
  - -b, --burden_res: Path to the result of burden shift from permutation test (\*.binom_pvals.txt.gz).
  - -c_info, --category_info: Path to a text file with category information (`*.category_info.txt`).
  - -o_dir, --output_directory: Path to the directory where the output files will be saved. By default, outputs will be saved at ``$CWAS_WORKSPACE``.
  - -c_set, --cat_set: Path to the category information file from binomial burden test (\*.category_info.txt).
  - -c_count, --cat_count: Path of the categories counts file from binomial burden test (\*.category_counts.txt).
  - -t, --tag: Tag used for the name of the output files. By default, None.
  - -c_cutoff, --count_cutoff: The number of cutoff for category counts. It must be positive value. By default, 7.
  - --pval: P-value threshold. By default, 0.05.

  .. code-block:: solidity
    
    cwas burden_shift -i $HOME/cwas_output/de_novo_variants.burden_test.txt \
    -b $HOME/cwas_output/de_novo_variants.binom_pvals.txt.gz \
    -o_dir $HOME/cwas_output \
    -c_info $HOME/cwas_output/de_novo_variants.category_info.txt \
    -c_count $HOME/cwas_output/de_novo_variants.category_counts.txt \
    -c_cutoff 7 \
    --pval 0.05

  Below are the output files generated.

  .. code-block:: solidity

    $HOME/cwas_output
    ...
    ├── de_novo_variants.burdenshift_p0.05_cutoff7.dist_plot.pdf
    ├── de_novo_variants.burdenshift_p0.05_cutoff7.result_plot.pdf
    ├── de_novo_variants.burdenshift_p0.05_cutoff7.txt
    ...



11.    :ref:`DAWN analysis <dawn>`
####################################

  Investigate the relationship between categories and identify the specific type of categories clustered within the network.

  The parameters of the command are as below:

  - -e, --eig_vector: Eigen vector file. This is the output file from :ref:`calculation of effective number of tests <effnumtest>`. The file name must have pattern ``*eig_vecs*.txt.gz``.
  - -c, --corr_mat: Category correlation matrix file. This is the output file from :ref:`categorization <categorization>`. The file name must have pattern ``*correlation_matrix*.zarr``.
  - -P, --permut_test: Permutation test file. This is the output file from :ref:`burden test <permtest>`. The file name must have pattern ``*permutation_test*.txt.gz``.
  - -o_dir, --output_directory: Path to the directory where the output files will be saved. By default, outputs will be saved at ``$CWAS_WORKSPACE``.
  - -r, --range: Range (i.e., (start,end)) to find optimal K for k-means clustering. It must contain two integers that are comma-separated. The first integer refers to the start number and must be above 1. The second integer refers to the end.
  - -k, --k_val: K for K-means clustering. With this argument, users can determine K manually. ``-r`` and ``-k`` arguments are mutually exclusive. If ``-k`` is given, ``-r`` will be ignored.
  - -s, --seed: Seed value for t-SNE. Same seed will generate same results for the same inputs.
  - -t, --tag: Tag used for the name of the output files. By default, None.
  - -c_count, --cat_count: Path of the categories counts file from burden test.
  - -C, --count_threshold: The treshold of variant (or sample) counts. The least amount of variants a category should have.
  - -R, --corr_threshold: The threshold of correlation values between clusters. Computed by the mean value of correlation values of categories within a cluster.
  - -S, --size_threshold: The threshold of the number of categories per cluster. The least amount of categories a cluster should have.
  - -p, --num_proc: Number of worker processes that will be used for the DAWN analysis. By default, 1.


  .. code-block:: solidity
  
      cwas dawn \
      -e $HOME/cwas_output/de_novo_variants.eig_vecs.TFBS.txt.gz \
      -c $HOME/cwas_output/de_novo_variants.correlation_matrix.zarr \
      -P $HOME/cwas_output/de_novo_variants.permutation_test.txt.gz \
      -o_dir $HOME/cwas_output \
      -r 2,200 \
      -s 123 \
      -t TFBS.exact \
      -c_count $HOME/cwas_output/de_novo_variants.category_counts.txt \
      -C 20 \
      -R 0.12 \
      -S 2 \
      -p 8

  Below are the output files generated.

  .. code-block:: solidity

    $HOME/cwas_output
    ...
    ├── TFBS.exact.cluster_annotation.csv
    ├── TFBS.exact.graph_layout.csv
    ├── TFBS.exact.iplot.igraph.pdf
    ├── TFBS.exact.iplot.igraph_with_community.pdf
    ├── TFBS.exact.iplot.igraph_with_number.pdf
    ├── TFBS.exact.ipvalue_fdr.txt
    ├── TFBS.exact.ipvalue_fdr_igraph.csv
    ├── TFBS.exact.ipvalue_fdr_ipvalue_risk.csv
    ├── TFBS.exact_choose_K_silhouette_score_plot.pdf
    ...

