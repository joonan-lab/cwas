###############################
Categorization
###############################

The variants are categorized according to the annotation results. To count the number of variants that falls into a category in each sample, CWAS-Plus checks whether a variant is annotated to an annotation term (for example, whether a variant is annotated to a gene that is one of disease risk genes). As these annotation terms are classified into five major groups, the combination of terms from each group results in a single category. While categorization, CWAS-Plus excludes the redundant categories that share the exact same variants (such as, missense variants that fall into a 'single-nucleotide variant (SNV)' term and 'All (SNV and insertion-deletion)' term are the same, as all missense variants are SNVs).

The parameters of the command are as below:

- -i, --input_file : Path to the annotated VCF, resulted from annotation process. This file could be gzipped or not.
- -o_dir, --output_directory : Path to the directory where the output files will be saved. By default, outputs will be saved at ``$CWAS_WORKSPACE``.
- -p, --num_proc : Number of worker processes that will be used for the categorization process. By default, 1.

.. code-block:: solidity

    cwas categorization -i INPUT.annotated.vcf -o_dir OUTPUT_DIR -p 8