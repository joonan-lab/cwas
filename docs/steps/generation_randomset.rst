===============================
Generation of random variants
===============================

To quantify the relationship between categories, CWAS-Plus generates multiple sets of random variants. With this command, a number of random variants were generated, with each set containing the same number of as the input variants. The number of variants in each family is conserved, while the generated variants are distributed to members of the family randomly.

The parameters of the command are as below:

- -i, --input_file : Path to the annotated VCF file, resulted from annotation process. This file could be gzipped or not.
- -o_dir, --output_directory : Path to the directory where the output files will be saved. By default, outputs will be saved at ``$CWAS_WORKSPACE/random-mutation``.
- -s, --sample_info : Path to the txt file containing the sample information for each sample. This file must have three columns (``SAMPLE``, ``FAMILY``, ``PHENOTYPE``) with the exact name.
- -t, --out_tag : Prefix of output files. Each output file name will start with this tag. By default, ``rand_mut``.
- -n, --num_sim : Number of simulations to generate random variants. By default, 1. This number of sets, with each set containing the equal number of variants as the input VCF file, will be created.
- -r, --resume : Resumes the simulation from the last step. Assumes some generated output files are not truncated. By default, False.


.. code-block:: solidity

    cwas simulation -i INPUT.annotated.vcf -s SAMPLE_LIST.txt -n 10000 -p 8


