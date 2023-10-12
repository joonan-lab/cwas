.. _correlation:

*********************************
Calculate correlation matrix
*********************************

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
- c_info, --category_info: Path to a text file with category information (`*.category_info.txt`).
- -d, --domain_list: Domain list to filter categories based on GENCODE domain. By default, `all`.

.. code-block:: solidity

    cwas correlation -v INPUT.annotated.vcf.gz -o_dir OUTPUT_DIR -c_info OUTPUT.category_info.txt -p 8 -cm variant -im