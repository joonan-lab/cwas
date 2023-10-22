.. _tips:

**********************
Tips
**********************


How to configure ANNOTATION_KEY_CONFIG yaml file
######################################################

The ANNOTATION_KEY_CONFIG yaml file contains *functional_score* and *functional_annotation*. Each file is in bed format and should be written with their short name that represents them. This name will be further used in other analyses.

When setting the short name, avoid using underscores ('_'). Underscores are used to distinguish different domains within a single category. For example, a category 'A_B_C_D_E' will be recognized as five domains, but if the name of the category is 'A_B_C_D_E_F', it will cause error while association testing, as the category is divided into six domains.

An example for ANNOTATION_KEY_CONFIG yaml file looks like below.

.. code-block:: solidity

  functional_score:
    bed1.bed.gz: annot1
    bed2.bed.gz: annot2
  functional_annotation:
    bed3.bed.gz: annot_3 # Do not use underscores like this. Users can use 'annot3' instead.
    bed4.bed.gz: annot4

The files should be located inside **ANNOTATION_DATA_DIR**. For preparation step, these files should be indexed using tabix.

As an example, users can sort and index their bed file like below.

.. code-block:: solidity
  
  cat bed1.bed | sort -k1,1V -k2,2n -k3,3n -t$'\t' | bgzip -c > sorted.bed1.bed.gz
  tabix -p bed sorted.bed1.bed.gz




How to add or remove a gene set
######################################################


How to add or remove a functional annotation or score
######################################################

