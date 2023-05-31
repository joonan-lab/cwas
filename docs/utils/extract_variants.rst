===================
Extract variants
===================

This process is to extract variants that are allocated to a specific category. When there is a category of interest, the users can extract the variants that belong to the category and look through.

The parameters of the command are as below:

- -i, --input_file: Path to the annotated VCF, resulted from annotation process. This file could be gzipped or not.
- -o, --output_directory: Path to the directory where the output files will be saved. By default, outputs will be saved at ``$CWAS_WORKSPACE``.
- -t, --tag: Tag used for the name of the output file. By default, None.
- -c, --category_set_path: Path to a txt file containing categories for extracting variants. By default, None. This file must contain ``Category`` column with the name of categories to be extracted.

  +-------------------------------------------------------+
  |Category                                               |
  +=======================================================+
  |All_CHD8Common_All_IntergenicRegion_EarlyCREMicro      |
  +-------------------------------------------------------+
  |All_CHD8Common_phastCons46way_PromoterRegion_EarlyCREL4|
  +-------------------------------------------------------+
  |All_DDD_All_PromoterRegion_EarlyCREOligo               |
  +-------------------------------------------------------+

- -ai, --annotation_info: Saves with annotation information attached (such as gene list, functional annotations, etc). By default, False.


**Without a category set**

.. code-block:: solidity

    cwas extract_variant -i INPUT.annotated.vcf -o_dir OUTPUT_DIR -t filtered -ai

**With a category set**

.. code-block:: solidity

    cwas extract_variant -i INPUT.annotated.vcf -o_dir OUTPUT_DIR -t filtered -ai -c CATEGORY_SET.txt -ai