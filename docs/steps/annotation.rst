.. _annotation:

*********************************
Annotation
*********************************

CWAS-Plus annotated variants using VEP and its own custom annotation algorithm. The parameters of the command are as below:

- -v, --vcf_file: Path to the input vcf file. This file could be gzipped or not.
- -p, --num_proc: Number of worker processes that will be used for the annotation process. By default, 1.
- -o_dir, --output_directory: Path to the directory where the output files will be saved. By default, outputs will be saved at ``$CWAS_WORKSPACE``.


.. code-block:: solidity

  cwas annotation -v INPUT.vcf -o_dir OUTPUT_DIR -p 8