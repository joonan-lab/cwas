.. _effnumtest:

********************************************
Calculate the effective number of tests
********************************************

For multiple comparisons, the effective number of tests is required. Each category is equal to each test, but categories could be correlated and these highly correlated categories can be defined as a single effective test. By assessing the correlation between categories, the users can gain the number of effective tests.

The output files from eigen decomposition will also be used for DAWN analysis. For preparing the input for DAWN analysis, calculating eigen values and vectors for categories of interest is recommended. For example, if the uesr is interested in promoter categories, only substracting promoter categories for eigen decomposition is required.

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



1. Find the number of effective tests

  - Only categories with a value (number of variants or samples) greater than or equal to cutoff are used. The cutoff is used to select informative significant tests with a sufficient number of variants (or samples).
        
    - With specified cutoff: Categories with a value (number of variants or samples) greater than or equal to **specified cutoff** are used.

    .. code-block:: solidity
            
        cwas effective_num_test -i INPUT.correlation_matrix.zarr -o_dir OUTPUT_DIR -if corr -n 10000 -ef -thr 8 -c_count INPUT.category_counts.txt

    - Without specified cutoff: The cutoff is automatically calculated and applied to filter categories with a value (number of variants or samples) greater than or equal to cutoff. The cutoff represents the minimum number of variants (or samples) required for a one-sided binomial test with p\<0.05, assuming the null hypothesis is a Binomial(m, No. cases/No. total samples) distribution with 1 mutation in controls and m-1 mutations in cases.

    .. code-block:: solidity
        
        cwas effective_num_test -i INPUT.correlation_matrix.zarr -o_dir OUTPUT_DIR -if corr -n 10000 -ef -c_count INPUT.category_counts.txt

2. Generate inputs for DAWN analysis

  - Use the identical cutoff for the number of variants (or samples) as in ``1. Find the number of effective tests``, while focusing only on specific categories relevant to the users' domains of interest. For example, users can exclusively use intergenic categories. Aditionally, when generating DAWN analysis inputs, **omit the** ``-ef`` **argument**, as the number of effective tests calculated for this subset of categories of interest will not be used elsewhere.

    .. code-block:: solidity
            
        cwas effective_num_test -i INPUT.correlation_matrix.zarr -o_dir OUTPUT_DIR -if corr -n 10000 -c CATEGORY_SET.txt -c_count INPUT.category_counts.txt

