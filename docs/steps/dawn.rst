.. _dawn:

*********************************
DAWN analysis
*********************************

The users can investigate the relationship between categories and identify the specific type of categories that are clustered within the network of categories of interest.

  - -e, --eig_vector: Eigen vector file. This is the output file from :ref:`calculation of effective number of tests <effnumtest>`. The file name must have pattern ``*eig_vecs*.zarr``.
  - -c, --corr_mat: Category correlation matrix file. This is the output file from :ref:`categorization <categorization>`. The file name must have pattern ``*correlation_matrix*.pkl``.
  - -P, --permut_test: Permutation test file. This is the output file from :ref:`burden test <permtest>`. The file name must have pattern ``*permutation_test*.txt.gz``.
  - -o_dir, --output_directory: Path to the directory where the output files will be saved. By default, outputs will be saved at ``$CWAS_WORKSPACE``.
  - -r, --range: Range (i.e., (start,end)) to find optimal K for k-means clustering. It must contain two integers that are comma-separated. The first integer refers to the start number and must be above 1. The second integer refers to the end.
  - -k, --k_val: K for K-means clustering. With this argument, users can determine K manually. ``-r`` and ``-k`` arguments are mutually exclusive. If ``-k`` is given, ``-r`` will be ignored.
  - -s, --seed: Seed value for t-SNE. Same seed will generate same results for the same inputs.
  - -T, --tsen_method: Gradient calculation algorithm for t-SNE, which is used in TSNE of sklearn. If the dataset is large, 'barnes_hut' is recommended. By default, exact.
  - -t, --tag: Tag used for the name of the output files. By default, None.
  - -c_count, --cat_count: Path of the categories counts file from burden test.
  - -C, --count_threshold: The treshold of variant (or sample) counts. The least amount of variants a category should have.
  - -R, --corr_threshold: The threshold of correlation values between clusters. Computed by the mean value of correlation values of categories within a cluster.
  - -S, --size_threshold: The threshold of the number of categories per cluster. The least amount of categories a cluster should have.
  - -p, --num_proc: Number of worker processes that will be used for the DAWN analysis. By default, 1.


.. code-block:: solidity
  
    cwas dawn -e INPUT_EIG_VEC \
    -c INPUT_CORR_MATRIX \
    -P INPUT_PERMUATION_RESULT \
    -o_dir OUTPUT_DIR \
    -r 2,500 \
    -s 123 \
    -t test \
    -c_count CATEGORY_COUNTS.txt \
    -C 20 \
    -R 0.7 \
    -S 2 \
    -p 8



