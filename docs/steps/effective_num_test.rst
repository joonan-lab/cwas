.. _effnumtest:

========================================
Calculate the effective number of tests
========================================

For multiple comparisons, the effective number of tests is required. Each category is equal to each test, but categories could be correlated and these highly correlated categories can be defined as a single effective test. By assessing the correlation between categories, the users can gain the number of effective tests.

The output files from eigen decomposition will also be used for DAWN analysis. For preparing the input for DAWN analysis, calculating eigen values and vectors for categories of interest is recommended. For example, if the uesr is interested in promoter categories, only substracting promoter categories for eigen decomposition is required.

The parameters of the command are as below:

- -i, --input_file: Path to the concatenated z-scores.
- -if, --input_format: Specify the format of the input file. Available options are ``corr``, ``inter``, and ``zscores``. By default, ``corr`` will be used. Each format refers to the following:

  - corr: A matrix with correlation values between categories.
  - inter: A matrix with intersected number of variants (or samples) between categories.
  - zscores: A matrix with concatenated z-scores.

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

- -ef, --eff_num_test: Calculate the effective number of tests. By default, false.



.. code-block:: solidity

    cwas effective_num_test -i INPUT.zscores.txt.gz -o_dir OUTPUT_DIR -t test -c CATEGORY_SET.txt -ef


