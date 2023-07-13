.. _dawn:

*********************************
DAWN analysis
*********************************

The users can investigate the relationship between categories and identify the specific type of categories that are clustered within the network of categories of interest.


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
- -c_count, --cat_count
- -CT, --count_threshold: The treshold of variant (or sample) counts. The least amount of variants a category should have.
- -CR, --corr_threshold: The threshold of correlation values between clusters. Computed by the mean value of correlation values of categories within a cluster.
- -S, --size_threshold: The threshold of the number of categories per cluster. The least amount of categories a cluster should have.
- -p, --num_proc: Number of worker processes that will be used for the DAWN analysis. By default, 1.


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



