.. _dawn:

###############################
DAWN analysis
###############################

The users can investigate the relationship between categories and identify the specific type of categories that are clustered within the network of categories of interest.


- -i, --input_file: Path to the categorized txt file, resulted from categorization process. This file could be gzipped or not.
- -o_dir, --output_directory: Path to the directory where the output files will be saved. By default, outputs will be saved at ``$CWAS_WORKSPACE``.



- -c, --category_set_path: Path to a text file containing categories for training. If not specified, all of the categories categorization file will be used. This file must contain ``Category`` column with the name of categories to be used.
- -t, --tag: Tag used for the name of the output files. By default, None.
- -u, --use_n_carrier: Enables the use of the number of samples with variants in each category for burden test instead of the number of variants. With this option, CWAS-Plus counts the number of samples that carry at least one variant of each category.

- -p, --num_proc: Number of worker processes that will be used for the permutation process. By default, 1.



.. code-block:: solidity
  
    cwas dawn -i INPUT.categorization_result.txt.gz \
    -o_dir OUTPUT_DIR


