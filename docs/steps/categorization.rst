.. _categorization:

*********************************
Categorization
*********************************

The variants are categorized according to the annotation results. To count the number of variants that falls into a category in each sample, CWAS-Plus checks whether a variant is annotated to an annotation term (for example, whether a variant is annotated to a gene that is one of disease risk genes). As these annotation terms are classified into five major groups, the combination of terms from each group results in a single category. While categorization, CWAS-Plus excludes the redundant categories that share the exact same variants (such as, missense variants that fall into a 'single-nucleotide variant (SNV)' term and 'All (SNV and insertion-deletion)' term are the same, as all missense variants are SNVs).

In this step, users can also generate two matrices, (1) a matrix that contains the number of variants (or samples, with --use_carrier option) that intersect between categories, (2) a matrix that contains correlation values between categories. The correlation matrix is computed from the intersected matrix (1). The users can choose one of the matrices for calculating the number of effective tests and DAWN analysis.

The parameters of the command are as below:

- -i, --input_file: Path to the annotated VCF, resulted from annotation process. This file contains a specific pattern of ``.annotated.vcf`` in the file name. This file could be gzipped or not.
- -o_dir, --output_directory: Path to the directory where the output files will be saved. By default, outputs will be saved at ``$CWAS_WORKSPACE``.
- -p, --num_proc: Number of worker processes that will be used for the categorization process. To prevent crashes caused by insufficient RAM when processing large input VCF files (e.g., over 10 million variants) using multiple cores, using small number of cores and monitoring the memory usage are recommended. By default, 1.
- -m, --matrix: Generate a correlation matrix and a matrix with intersected number of variants (or samples) between every two categories. Available options are ``variant`` or ``sample``. By default, False.

  - variant: Use the intersected number of variants between two categories.
  - sample: Use the intersected number of samples between two categories.

.. code-block:: solidity

    cwas categorization -i INPUT.annotated.vcf.gz -o_dir OUTPUT_DIR -p 8

    cwas categorization -i INPUT.annotated.vcf.gz -o_dir OUTPUT_DIR -p 8 -m variant