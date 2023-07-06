================================
Advanced tutorial for CWAS-Plus
================================

This is an advanced tutorial for CWAS-Plus. Specific descriptions of arguments are described in the page of each step.



0. Data requirements

  1. Input vcf data (variant list)

  Prepare variants in vcf format. The order of the columns should follow the specification of `vcf <https://samtools.github.io/hts-specs/VCFv4.2.pdf>`_. The INFO field must contain a sample ID of each variant with this format ``SAMPLE={sample_id}``.

  .. code-block:: solidity

    #CHROM  POS ID  REF ALT QUAL    FILTER  INFO
    chr1    3747728 .        T       C       .       .       SAMPLE=11000.p1;BATCH=P231
    chr1    38338861        .       C       A       .       .       SAMPLE=11000.p1;BATCH=P231
    chr1    117942118       .      T       G       .       .       SAMPLE=11000.p1;BATCH=P231


  2. Sample information

  Prepare sample information in txt format. The file must be tab separated. It also must contain three columns, *SAMPLE*, *FAMILY*, and *PHENOTYPE*. A value in the *PHENOTYPE* muse be *case* or *ctrl*.
  The values in the SAMPLE column must be matched to the sample IDs of variants in the input vcf data.

  +----------+--------+-----------+
  |  SAMPLE  | FAMILY | PHENOTYPE |
  +==========+========+===========+
  | 11000.p1 | 11000  |   case    |
  +----------+--------+-----------+
  | 11000.s1 | 11000  |   ctrl    |
  +----------+--------+-----------+
  | 11002.p1 | 11002  |   case    |
  +----------+--------+-----------+
  | 11002.s1 | 11002  |   ctrl    |
  +----------+--------+-----------+

  3. Adjustment factors

  Adjustment factors are required if the users want to adjust the number of variants for each sample in CWAS-Plus. The file must be tab separated and must contain two columns, *SAMPLE* and *AdjustFactor*. A value in the *AdjustFactor* must be a float.
  The values in the SAMPLE column must be matched to the sample IDs of variants in the input vcf data.

  +----------+--------------+
  | SAMPLE   | AdjustFactor |
  +==========+==============+
  | 11000.p1 | 0.932        |
  +----------+--------------+
  | 11000.s1 | 1.082        |
  +----------+--------------+
  | 11002.p1 | 0.895        |
  +----------+--------------+
  | 11002.s1 | 1.113        |
  +----------+--------------+


  For example run, the above data are available at `joonan-lab/cwas-input-example <https://github.com/joonan-lab/cwas-input-example>`_.


  4. Annotation dataset

  CWAS-Plus requires annotation dataset to annotate and categorize variants. Users can customize their own annotation dataset based on ther interest.
  For example run, annotation datasets for investigating autism spectrum disorder are available at `joonan-lab/cwas-dataset <https://github.com/joonan-lab/cwas-dataset>`_.

  .. code-block:: solidity

    git clone https://github.com/joonan-lab/cwas-dataset.git
  



1. :ref:`Install CWAS-Plus <installation>`

  The users can install CWAS-Plus through Github.

  Create a conda environment with required installations through command below.

  .. code-block:: solidity
    
    git clone https://github.com/joonan-lab/cwas.git
    cd cwas
    conda env create -f environment.yml -n cwas


  After creating conda environment, activate the environment and install CWAS-Plus. If users want to force the installation, they can use ``-f`` option.

  .. code-block:: solidity
    
    conda activate cwas
    python setup.py install


  CWAS-Plus requires a working directory for efficiency. Users can create the working directory through command below.

  .. code-block:: solidity
    
    cwas start

  By default, the command creates a working directory (``.cwas``) in the home directory. However, if users want to set the working directory manually, they can use ``-w`` option to specify the path of the desired working directory.

  .. code-block:: solidity
    
    cwas start -w /path/to/the/working/directory

  The command ``cwas start``, also creates a configuration file inside the working directory. If there is a pre-installed VEP, the path of the VEP in the configuration file will be automatically set.


2. :ref:`Configuration <configuration>`

  Inside the CWAS working directory, there is a configuration file (``configuration.txt``). This file is needed for retrieving the path of specific files needed for CWAS-Plus run.
  With pre-installed VEP, the configuration file looks like below.


  .. code-block:: solidity
    
    ANNOTATION_DATA_DIR=
    GENE_MATRIX=
    ANNOTATION_KEY_CONFIG=
    VEP=/home/cwas_testing/miniconda3/envs/cwas/bin/vep
    VEP_CACHE_DIR=
    VEP_CONSERVATION_FILE=
    VEP_LOFTEE=
    VEP_HUMAN_ANCESTOR_FA=
    VEP_GERP_BIGWIG=
    VEP_MIS_DB=
    VEP_MIS_INFO_KEY=
    VEP_MIS_THRES=


  Fill in the file with paths of the required tools and data.

  - **ANNOTATION_DATA_DIR**: This is the path of the directory, which contains annotation datasets, such as bed files.
  - **GENE_MATRIX**: This is the path of the gene matrix, which is a text file. The first column should be gene ID, and the second column should be gene name. The other columns will represent each gene list and show whether each row (=gene) are matched to the gene list or not by a binary code (0, 1). 1 if the gene is matched to a gene list, 0 if not.
  - **ANNOTATION_KEY_CONFIG**: This is the path of the annotation key file, which is a yaml file. This file contains the name of the annotation datasets inside the annotation dataset directory and the key names that will be used to represent the dataset. All details should be written in yaml syntax. Also, to split the category group to functional score and functional annotation, the users should type each annotation dataset under the matched group dictionary. Below is an example of this file. The format should be (name): (key) with a uniform indentation for each row. Be aware that the name of the annotations should not contain '_'. As domains will combined with '_' as a delimiter, using '_' in the annotation name will cause errors.
  - **VEP**: This is the path of VEP. If there is a pre-installed VEP, this line would be written in advance when the users typed the command ``cwas start``.
  - **VEP_CONSERVATION_FILE**: This is the path of the conservation file (`loftee.sql`), which will be used for variant classification.
  - **VEP_LOFTEE**: This is the path of the directory of loftee plugin, which will be used for variant classification.
  - **VEP_HUMAN_ANCESTOR_FA**: This is the path of the human ancestor fasta file, which will be used for variant classification.
  - **VEP_GERP_BIGWIG**: This is the path of the GERP bigwig file, which will be used for variant classification.
  - **VEP_MIS_DB**: This is the path of the database in vcf format. This will be used for variant classification. Users can manually prepare this file to classify damaging missense variants.
  - **VEP_MIS_INFO_KEY**: The name that will be used for the missense classification database.
  - **VEP_MIS_THRES**: The cutoff that will be used for the missense classification. The missense variants scoring equal to or above *VEP_MIS_THRES* will be classified as damaging missense mutations.


  Each data is available as below.

  - By default, CWAS-Plus provides *ANNOTATION_DATA_DIR*, *GENE_MATRIX* and *ANNOTATION_KEY_CONFIG* through `joonan-lab/cwas-dataset <https://github.com/joonan-lab/cwas-dataset>`_.
  - *VEP* can be installed through github or conda. The command to install VEP through conda is as below.

  .. code-block:: solidity

    conda install -c bioconda ensembl-vep


  To download VEP resources in GRCh38 version in one step, run the command below.

  .. code-block:: solidity

    Run.

  Specific paths for downloading resources are as below.


  - *VEP_CONSERVATION_FILE* can be downloaded as below.

  .. code-block:: solidity
    
    wget https://personal.broadinstitute.org/konradk/loftee_data/GRCh38/loftee.sql.gz
    gunzip -k loftee.sql.gz

  - *VEP_LOFTEE* can be downloaded as below.

  .. code-block:: solidity
    
    cd /home/cwas_testing/.vep/Plugins
    git clone -b grch38 https://github.com/konradjk/loftee.git

  - *VEP_HUMAN_ANCESTOR_FA* can be downloaded as below.

  .. code-block:: solidity
    
    wget https://personal.broadinstitute.org/konradk/loftee_data/GRCh38/human_ancestor.fa.gz
    wget https://personal.broadinstitute.org/konradk/loftee_data/GRCh38/human_ancestor.fa.gz.fai
    wget https://personal.broadinstitute.org/konradk/loftee_data/GRCh38/human_ancestor.fa.gz.gzi

  - *VEP_GERP_BIGWIG* can be downloaded as below.

  .. code-block:: solidity

    wget https://personal.broadinstitute.org/konradk/loftee_data/GRCh38/gerp_conservation_scores.homo_sapiens.GRCh38.bw






   For example run, copy the ``configuration.txt`` in the ``cwas-dataset`` to the working directory. The file should be as below.

  .. code-block:: solidity
    
    ANNOTATION_DATA_DIR=$HOME/cwas-dataset
    GENE_MATRIX=$HOME/cwas-dataset/gene_matrix.txt
    ANNOTATION_KEY_CONFIG=$HOME/cwas-dataset/annotation_keys.yaml
    VEP=$HOME/miniconda3/envs/cwas/bin/vep
    VEP_CACHE_DIR=$HOME/.vep
    VEP_CONSERVATION_FILE=$HOME/.vep/loftee.sql
    VEP_LOFTEE=$HOME/.vep/Plugins/loftee
    VEP_HUMAN_ANCESTOR_FA=$HOME/.vep/human_ancestor.fa.gz
    VEP_GERP_BIGWIG=$HOME/.vep/gerp_conservation_scores.homo_sapiens.GRCh38.bw
    VEP_MIS_DB=$HOME/.vep/MPC_hg38.vcf.bgz
    VEP_MIS_INFO_KEY=MPC
    VEP_MIS_THRES=2



  After filling the configuration file, ``cwas configuration`` command will create symlinks of annotation datasets into the working directory and fill the ``.cwas_env`` file in the home directory for storing environmental variables.

  .. code-block:: solidity

    cwas configuration

3. :ref:`Prepare annotation datasets <data-prep-label>`

  Gather and merge functional annotations and scores into a single bed file.

  .. code-block:: solidity

    cwas preparation -p 8

4. :ref:`Annotation <annotation>`

  Annotate the input VCF file with VEP and bed custom annotation algorithm.

  .. code-block:: solidity

    cwas annotation -v INPUT.vcf -o_dir OUTPUT_DIR -n 8

5. :ref:`Categorization <categorization>`

  Categorize variants into groups based on the annotation datasets. A single category is a combination of five domains (i.e., variant type, gene biotype, gene list, functional annotation and functional score). Details are provided in the :ref:`Overview of annotation datasets <overview>`.

  .. code-block:: solidity

    cwas categorization -i INPUT.annotated.vcf -o_dir OUTPUT_DIR -p 8

6. :ref:`Burden test <burdentest>`

  Calculate the burden of each category by comparing the number of variants per case and control. Two types of tests are used for p-value calculation: binomial test and permutation test.
   
  - Binomial test

     .. code-block:: solidity

        cwas binomial_test -i INPUT.categorization_result.txt.gz -o_dir OUTPUT_DIR -s SAMPLE_LIST.txt -a ADJUST_FACTOR.txt

  - Permutation test
   
     .. code-block:: solidity

        cwas permutation_test -i INPUT.categorization_result.txt.gz -o_dir OUTPUT_DIR -s SAMPLE_LIST.txt -a ADJUST_FACTOR.txt -n 10000 -p 8 -b


7. :ref:`Caculate the correlation matrix <categorization>`

  Caculate the correlation matrix from the intersected number of variants (or samples) between every two categories.

  .. code-block:: solidity

    cwas categorization -i INPUT.annotated.vcf -o_dir OUTPUT_DIR -p 8 -m variant


8.  :ref:`Calculate the number of effective tests <effnumtest>`

  From correlation matrix, compute eigen values and vectors. Based on these outputs, users can calculate the number of effective tests.

  .. code-block:: solidity

    cwas effective_num_test -i INPUT.correlation_matrix.pkl -o_dir OUTPUT_DIR -t test -c CATEGORY_SET.txt -ef


9.  :ref:`Risk score analysis <riskscore>`

  Identify the best predictor of the phenotype by training Lasso regression model with the number of variants within each category across samples.

  .. code-block:: solidity

    cwas risk_score -i INPUT.categorization_result.txt.gz \
    -o_dir OUTPUT_DIR \
    -s SAMPLE_LIST.txt \
    -a ADJUST_FACTOR.txt \
    -c CATEGORY_SET.txt \
    -thr 3 \
    -tf 0.7 \
    -n_reg 10 \
    -f 5 \
    -n 1000 \
    -p 8


10.  :ref:`Burden shift analysis <riskscore>`

  Identify the overrepresented domains associated to the phenotype.

  .. code-block:: solidity


11.  :ref:`DAWN analysis <dawn>`

  Investigate the relationship between categories and identify the specific type of categories clustered within the network.

  .. code-block:: solidity
  
      cwas dawn -i_dir INPUT_DIR \
      -o_dir OUTPUT_DIR \
      -r 2,500 \
      -s 123 \
      -t test \
      -c CATEGORY_SET.txt \
      -c_count CATEGORY_COUNTS.txt \
      -CT 2 \
      -CR 0.7 \
      -S 20 \
      -p 8


