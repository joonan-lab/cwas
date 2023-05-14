.. _effnumtest:

========================================
Calculate the effective number of tests
========================================

To obtain the measure of category-wide significance, the effective number of tests is required. Each category is equal to each test, but categories could be correlated and these highly correlated categories can be defined as a single effective test. By assessing the correlation between categories, the users can gain the number of effective tests.

The output files from eigen decomposition will also be used for DAWN analysis. For preparing the input for DAWN analysis, calculating eigen values and vectors for categories of interest is recommended. For example, if the uesr is interested in promoter categories, only substracting promoter categories for eigen decomposition is required.

The parameters of the command are as below:

- -i, --input_file : Path to the concatenated z-scores.
- -o_dir, --output_directory : Path to the directory where the output files will be saved. By default, outputs will be saved at ``$CWAS_WORKSPACE``.
- -t, --tag : Tag used for the name of the output files. By default, None.
- -c, --category_set_path : Path to a text file containing categories for eigen decomposition. If not specified, all of the categories in the z-score file will be used. This file must contain ``Category`` column with the name of categories to be used.

  +-------------------------------------------------------+
  |Category                                               |
  +=======================================================+
  |All_CHD8Common_All_IntergenicRegion_EarlyCREMicro      |
  +-------------------------------------------------------+
  |All_CHD8Common_phastCons46way_PromoterRegion_EarlyCREL4|
  +-------------------------------------------------------+
  |All_DDD_All_PromoterRegion_EarlyCREOligo               |
  +-------------------------------------------------------+

- -ef, --eff_num_test : Calculate the effective number of tests. By default, false.

.. code-block:: solidity

    cwas effective_num_test -i INPUT.zscores.txt.gz -o_dir OUTPUT_DIR -t test -c CATEGORY_SET.txt -ef