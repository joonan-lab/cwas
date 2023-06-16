.. _concatzscores:

======================
Concatenate z-scores
======================

Gather p-values of categories from each simulation and convert p-values to z-scores. With this command, z-scores will be saved into a single txt.gz file.

The parameters of the command are as below:

- -i_dir, --input_directory: Path to the directory where burden test results are saved.
- -o_dir, --output_directory: Path to the directory where the output files will be saved. By default, outputs will be saved at ``$CWAS_WORKSPACE``.
- -s , --sample_info: Path to the txt file containing the sample information for each sample. This file must have three columns (``SAMPLE``, ``FAMILY``, ``PHENOTYPE``) with the exact name.
- -t, --tag: Tag used for the name of the output file. By default, None.
- -c, --category_set_path: Path to a txt file containing categories. Only z-scores of these categories will be concatenated. By default, None. This file must contain ``Category`` column with the name of categories to be used.

  +-------------------------------------------------------+
  |Category                                               |
  +=======================================================+
  |All_CHD8Common_All_IntergenicRegion_EarlyCREMicro      |
  +-------------------------------------------------------+
  |All_CHD8Common_phastCons46way_PromoterRegion_EarlyCREL4|
  +-------------------------------------------------------+
  |All_DDD_All_PromoterRegion_EarlyCREOligo               |
  +-------------------------------------------------------+

- -p , --num_proc: Number of worker processes that will be used for the concatenation process. By default, 1.



.. code-block:: solidity
  
    cwas concat_zscore -i_dir INPUT_DIR -o_dir OUTPUT_DIR -s SAMPLE_LIST.txt -c CATEGORY_SET.txt -p 8

