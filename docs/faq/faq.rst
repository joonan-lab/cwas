.. _faq:

**********************
FAQ
**********************

Q. 
###########################

When preparing the ANNOTATION_KEY_CONFIG yaml file, please avoid using underscores ('_') in the annotation name. Underscores are used for distinguishing different domains within a single category.

For example, check below.

.. code-block:: solidity

  functional_score:
    bed1.bed.gz: annot1
    bed2.bed.gz: annot2
  functional_annotation:
    bed3.bed.gz: annot_3 # Do not use underscores like this. Users can use 'annot3' instead.
    bed4.bed.gz: annot4


After filling the configuration file, type the below command for configuration. This process will create a symlink to the annotation dataset directory, gene matrix and the annotation key file to the user's workspace. Also, based on the annotation key file, a category domain file and a redundant category file will be created. The category domain file contains all the inferior category groups that will be used for CWAS-Plus. The redundant category file contains the combination of categories that will be excluded in CWAS-Plus. This is for removing duplicated categories (for example, coding variants with all genes and coding variants with coding genes) and nonsense categories (for example, missense variants that are indels).

To force configuration (overwrite previous configurations), use ``-f`` option.

.. code-block:: solidity

    cwas configuration


After configuration, a file ``.cwas_env`` that contains environmental variables for CWAS-Plus will be created in the home directory.

