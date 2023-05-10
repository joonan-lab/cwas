====================
Configuration
====================

Configuration step is required for setting the path the required datasets and environmental variables.

-------------------------------------
Setting the environmental variables
-------------------------------------

Inside the CWAS-Plus workspace directory, there is a configuration file ``configuration.txt``. This file contains a few crucial environmental variables that will be used through CWAS-Plus. The description of each variable is as below:

- **ANNOTATION_DATA_DIR** : This is the path of the directory, which contains annotation datasets, such as bed files.
- **GENE_MATRIX** : This is the path of the gene matrix, which is a text file. The first column should be gene ID, and the second column should be gene name. The other columns will represent each gene list and show whether each row (=gene) are matched to the gene list or not by a binary code (0, 1). 1 if the gene is matched to a gene list, 0 if not.
- **ANNOTATION_KEY_CONFIG** : This is the path of the annotation key file, which is a yaml file. This file contains the name of the annotation datasets inside the annotation dataset directory and the key names that will be used to represent the dataset. All details should be written in yaml syntax. Also, to split the category group to functional score and functional annotation, the users should type each annotation dataset under the matched group dictionary. Below is an example of this file. The format should be (name) : (key) with a uniform indentation for each row.
- **VEP** : This is the path of VEP. If there is a pre-installed VEP, this line would be written in advance when the users typed the command ``cwas start``.
- **VEP_CONSERVATION_FILE** : This is the path of the conservation file (`loftee.sql`), which will be used for variant classification.
- **VEP_LOFTEE** : This is the path of the directory of loftee plugin, which will be used for variant classification.
- **VEP_HUMAN_ANCESTOR_FA** : This is the path of the human ancestor fasta file, which will be used for variant classification.
- **VEP_GERP_BIGWIG** : This is the path of the GERP bigwig file, which will be used for variant classification.
- **VEP_MPC** : This is the path of the MPC file, which is a text file. This will be used for variant classification.


After filling the configuration file, type the below command for configuration. This process will create a symlink to the annotation dataset directory, gene matrix and the annotation key file to the user's workspace. Also, based on the annotation key file, a category domain file and a redundant category file will be created. The category domain file contains all the inferior category groups that will be used for CWAS-Plus. The redundant category file contains the combination of categories that will be excluded in CWAS-Plus. This is for removing duplicated categories (for example, coding variants with all genes and coding variants with coding genes) and nonsense categories (for example, missense variants that are indels).

.. code-block:: solidity

    cwas configuration


After configuration, a file ``.cwas_env`` that contains environmental variables for CWAS-Plus will be created in the home directory.


------------------
Data preparation
------------------

For efficient annotation process, the users should merge all bed files by typing the below command. During this process, all bed files will be split into their intersected or non-intersected intervals with numbers that indicate which annotation datasets are matched to the interval in binary scale.

The parameters of the command are as below:

- p : The number of processors.

.. code-block:: solidity

    cwas preparation -p 8