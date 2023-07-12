*********************************
Advanced tutorial for CWAS-Plus
*********************************

This is an advanced tutorial for CWAS-Plus. Specific descriptions of arguments are described in the page of each step.



0. Data requirements
#####################

  1. Input vcf file (variant list)

  Prepare sorted variants in vcf format. The order of the columns should follow the specification of `vcf <https://samtools.github.io/hts-specs/VCFv4.2.pdf>`_. The INFO field must contain a sample ID of each variant with this format ``SAMPLE={sample_id}``.

  .. code-block:: solidity

    #CHROM  POS ID  REF ALT QUAL    FILTER  INFO
    chr1    3747728 .        T       C       .       .       SAMPLE=11000.p1;BATCH=P231
    chr1    38338861        .       C       A       .       .       SAMPLE=11000.p1;BATCH=P231
    chr1    117942118       .      T       G       .       .       SAMPLE=11000.p1;BATCH=P231


  2. Sample information

  Prepare sample information in txt format. The file must be tab separated. It also must contain three columns, *SAMPLE*, *FAMILY*, and *PHENOTYPE*. A value in the *PHENOTYPE* muse be *case* or *ctrl*.
  The values in the SAMPLE column must be matched to the sample IDs of variants in the input vcf file.

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
  The values in the SAMPLE column must be matched to the sample IDs of variants in the input vcf file.

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

  .. code-block:: solidity

    git clone https://github.com/joonan-lab/cwas-input-example.git



  4. Annotation dataset

  CWAS-Plus requires annotation dataset to annotate and categorize variants. Users can customize their own annotation dataset based on ther interest.
  For example run, annotation datasets for investigating autism spectrum disorder are available at `joonan-lab/cwas-dataset <https://github.com/joonan-lab/cwas-dataset>`_.

  .. code-block:: solidity

    git clone https://github.com/joonan-lab/cwas-dataset.git
  



1. :ref:`Install CWAS-Plus <installation>`
############################################

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
############################################

  Inside the CWAS working directory, there is a configuration file (``configuration.txt``). This file is needed for retrieving the path of specific files needed for CWAS-Plus run.
  With pre-installed VEP, the configuration file looks like below.


  .. code-block:: solidity
    
    ANNOTATION_DATA_DIR=
    GENE_MATRIX=
    ANNOTATION_KEY_CONFIG=
    VEP=/path/to/VEP
    VEP_CACHE_DIR=
    VEP_CONSERVATION_FILE=
    VEP_LOFTEE=
    VEP_HUMAN_ANCESTOR_FA=
    VEP_GERP_BIGWIG=
    VEP_MIS_DB=
    VEP_MIS_INFO_KEY=
    VEP_MIS_THRES=


  The descriptions of each path are as follows.

  - **ANNOTATION_DATA_DIR**: This is the path of the directory, which contains annotation datasets, such as bed files.
  - **GENE_MATRIX**: This is the path of the gene matrix, which is a text file. The first column should be gene ID, and the second column should be gene name. The other columns will represent each gene list and show whether each row (=gene) are matched to the gene list or not by a binary code (0, 1). 1 if the gene is matched to a gene list, 0 if not.
  - **ANNOTATION_KEY_CONFIG**: This is the path of the annotation key file, which is a yaml file. This file contains the name of the annotation datasets inside the annotation dataset directory and the key names that will be used to represent the dataset. All details should be written in yaml syntax. Also, to split the category group to functional score and functional annotation, the users should type each annotation dataset under the matched group dictionary. Below is an example of this file. The format should be (name): (key) with a uniform indentation for each row. Be aware that the name of the annotations should not contain '_'. As domains will combined with '_' as a delimiter, using '_' in the annotation name will cause errors.
  - **VEP**: This is the path of VEP. If there is a pre-installed VEP, this line would be written in advance when the users typed the command ``cwas start``.
  - **VEP_CONSERVATION_FILE**: This is the path of the conservation file (`loftee.sql`), which will be used for variant classification.
  - **VEP_LOFTEE**: This is the path of the directory of loftee plugin, which will be used for variant classification.
  - **VEP_HUMAN_ANCESTOR_FA**: This is the path of the human ancestor fasta file, which will be used for variant classification.
  - **VEP_GERP_BIGWIG**: This is the path of the GERP bigwig file, which will be used for variant classification.
  - **VEP_MIS_DB**: This is the path of the database in vcf format. This will be used for variant classification. Users can manually prepare this file to classify damaging missense variants.
  - **VEP_MIS_INFO_KEY**: The name of the score in the missense classification database. It must be present in the INFO field of the database. The score must be specified by this name in the field. For example, if the user is using MPC score in the database, the database will look like below.
  
    +------+------+----+-----+-----+-----+--------+-----------+
    |#CHROM| POS  |  ID| REF |  ALT| QUAL| FILTER |INFO       |
    +------+------+----+-----+-----+-----+--------+-----------+
    |chr1  | 69094|  . | G   |  A  | .   | .      |MPC=2.73403|
    +------+------+----+-----+-----+-----+--------+-----------+
    |chr1  | 69094|  . | G   |  C  | .   | .      |MPC=2.29136|
    +------+------+----+-----+-----+-----+--------+-----------+
    |chr1  | 69094|  . | G   |  T  | .   | .      |MPC=2.29136|
    +------+------+----+-----+-----+-----+--------+-----------+
    |chr1  | 69095|  . | T   |  A  | .   | .      |MPC=4.31666|
    +------+------+----+-----+-----+-----+--------+-----------+

  - **VEP_MIS_THRES**: The cutoff that will be used for the missense classification. The missense variants scoring equal to or above *VEP_MIS_THRES* will be classified as damaging missense mutations.


  
  By default, CWAS-Plus provides all of the data above (except for VEP) and configuration file (``configuration.txt``) through `joonan-lab/cwas-dataset <https://github.com/joonan-lab/cwas-dataset>`_. Please note that the provided data serves as default examples, which users can customize to their specific needs.
  
  - *VEP* can be installed through github or conda. The command to install VEP through conda is as below.

  .. code-block:: solidity

    conda install -c bioconda ensembl-vep


  To download required resources and annotation datasets in GRCh38 version in one step, run the command below. It will create a directory (``$HOME/.vep``) and download resources in the directory.

  .. code-block:: solidity

    cd $HOME
    git clone https://github.com/joonan-lab/cwas-dataset.git
    cd cwas-dataset
    sh download_vep_resources.sh

  The downloading process might take a while.

  The descriptions of the files in the cwas-dataset are as below.

  - *annotation_keys.yaml*: List of annotation datasets with the exact file names and short names used for CWAS-Plus annotation.
  - *gene_matrix.txt*: List of genes with their functional annotations.
  - *download_vep_resources.sh*: Code to download VEP resources.
  - *configuration.txt*: Configuration file for CWAS-Plus specifying VEP path and required resources.
  - *functional_annotations.tar.gz*: BED files for annotating variants. After decompressing, **please move the files within the directory to the parent directory "cwas-dataset."**
  - *MPC_hg38.vcf.bgz*: Database for annotation damaging missense variants. For further information, please refer to the provided reference.
  - BED files for vertebrate conservation scores

    - PhyloP46way and PhastCons46Way
    - Due to the large file sizes, we provide an alternative download link for the original files.


  After preparing all resources, fill in the ``configuration.txt`` file with specific paths to the file.

  For example run, you can copy the ``configuration.txt`` in the ``cwas-dataset`` to the working directory. The file should be as below.
  
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
    VEP_MIS_DB=$HOME/cwas-dataset/MPC_hg38.vcf.bgz
    VEP_MIS_INFO_KEY=MPC
    VEP_MIS_THRES=2

  Please check the VEP path and modify *VEP* with the exact path.

  After filling the configuration file, ``cwas configuration`` command will create symlinks of annotation datasets into the working directory.
  The command will also add environment variables for CWAS-Plus in the ``.cwas_env`` file in the home directory. 

  .. code-block:: solidity

    cwas configuration

1. :ref:`Prepare annotation datasets <data-prep-label>`
############################################################

  Gather and merge functional annotations and scores into a single bed file. The annotation datasets in the *ANNOTATION_DATA_DIR* will be merged to a single bed file in the working directory.
  
  The parameters of the command are as below:

   - p: The number of processors.

  .. code-block:: solidity

    cwas preparation -p 8

4. :ref:`Annotation <annotation>`
############################################

  Annotate the input VCF file with VEP and bed custom annotation algorithm.
  When using more than one worker processes, CWAS-Plus automatically gzip and indexes non-gzipped input files for efficient multiprocessing.
  Output files are stored in the designated output directory (``-o_dir``) or, by default, in the working directory (``$CWAS_WORKSPACE``).

  The parameters of the command are as below:

   - -v, --vcf_file: Path to the input vcf file. This file could be gzipped or not.
   - -n, --num_cores: Number of worker processes that will be used for the annotation process. By default, 1.
   - -o_dir, --output_directory: Path to the directory where the output files will be saved. By default, outputs will be saved at ``$CWAS_WORKSPACE``.

  .. code-block:: solidity

    cwas annotation -v INPUT.vcf -o_dir OUTPUT_DIR -n 8

  The specific descriptions of the output files are as below. Each output file containing a specific pattern (i.e., ``.vep.vcf.gz``, ``.vep.vcf.gz.tbi``, ``.annotated.vcf``) in the file name as below will be found in the output directory.

  - OUTPUT.vep.vcf.gz: VEP annotated output file. This file is an intermediate output that has not been annotated with bed annotation files yet.
  - OUTPUT.vep.vcf.gz.tbi: Index file of the OUTPUT.vep.vcf.gz.
  - OUTPUT.annotated.vcf: The final output file. This file will be used as an input for categorization process.

5. :ref:`Categorization <categorization>`
############################################

  Categorize variants into groups based on the annotation datasets. A single category is a combination of five domains (i.e., variant type, gene biotype, gene list, functional annotation and functional score). Details are provided in the :ref:`Overview of annotation datasets <overview>`.

  The input file is the final output file resulted from annotation process. If users want to generate a matrix that contains correlation values between every two CWAS-Plus categories, they can use ``-m`` option. With this option, users must specify whether they want to calculate the correlation in variant-level (``-m variant``) or sample-level (``-m sample``). The generated correlation matrix will be used to calculate the number of effective tests for multiple comparisons.

  The parameters of the command are as below:

   - -i, --input_file: Path to the annotated VCF, resulted from annotation process. This file contains a specific pattern of ``.annotated.vcf`` in the file name. This file could be gzipped or not.
   - -o_dir, --output_directory: Path to the directory where the output files will be saved. By default, outputs will be saved at ``$CWAS_WORKSPACE``.
   - -p, --num_proc: Number of worker processes that will be used for the categorization process. To prevent crashes caused by insufficient RAM when processing large input VCF files (e.g., over 10 million variants) using multiple cores, using small number of cores and monitoring the memory usage are recommended. By default, 1.
   - -m, --matrix: Generate a correlation matrix and a matrix with intersected number of variants (or samples) between every two categories. Available options are ``variant`` or ``sample``. By default, False.

     - variant: Use the intersected number of variants between two categories.
     - sample: Use the intersected number of samples between two categories.

  .. code-block:: solidity

    cwas categorization -i INPUT.annotated.vcf -o_dir OUTPUT_DIR -p 8

  Caculate the correlation matrix from the intersected number of variants (or samples) between every two categories.

  .. code-block:: solidity

    cwas categorization -i INPUT.annotated.vcf -o_dir OUTPUT_DIR -p 8 -m variant


  The specific descriptions of the output files are as below. Each output file containing a specific pattern (i.e., ``.categorization_result.txt.gz``, ``.intersection_matrix.pkl``, ``.correlation_matrix.pkl``) in the file name as below will be found in the output directory.

  - OUTPUT.categorization_result.txt.gz: The final output file containing the number of variants in each category across samples. This file will be used as input in the burden test process.
  - OUTPUT.intersection_matrix.pkl: The matrix containing the number of intersected variants (or samples) between every two categories. This file will be generated only with ``-m`` option given.
  - OUTPUT.correlation_matrix.pkl: The matrix containing the correlation values between every two categories. This file will be generated only with ``-m`` option given. This file will be used for :ref:`calculating the number of effective tests <effnumtest>`. This file will be used as an input for :ref:`DAWN analysis <dawn>`.


6. :ref:`Burden test <burdentest>`
############################################

  Calculate the burden of each category by calculating the burden of each category by comparing the rate of variants per cases and the rate of variants per controls.
  
  For burden measurement, the package uses relative risk (RR), which is calculated by comparing the number of variants per phenotype group (RR>1, case burden; RR<1, control burden). The burden test in CWAS-Plus contains two types of p-value computation methods, binomial test and permutation test, to find more accurate p statistics.
   
  - Binomial test

    - -i, --input_file: Path to the categorized txt file, resulted from categorization process. This file could be gzipped or not.
    - -o_dir, --output_directory: Path to the directory where the output files will be saved. By default, outputs will be saved at ``$CWAS_WORKSPACE``.
    - -s, --sample_info: Path to the txt file containing the sample information for each sample. This file must have three columns (``SAMPLE``, ``FAMILY``, ``PHENOTYPE``) with the exact name.
    - -a, --adjustment_factor: Path to the txt file containing the adjust factors for each sample. This is optional. With this option, CWAS-Plus multiplies the number of variants (or carriers, in -u option) with the adjust factor per sample.
    - -u, --use_n_carrier: Enables the use of the number of samples with variants in each category for burden test instead of the number of variants. With this option, CWAS-Plus counts the number of samples that carry at least one variant of each category.

     .. code-block:: solidity

        cwas binomial_test -i INPUT.categorization_result.txt.gz -o_dir OUTPUT_DIR -s SAMPLE_LIST.txt -a ADJUST_FACTOR.txt

  - Permutation test
   
    - -i, --input_file: Path to the categorized txt file, resulted from categorization process. This file could be gzipped or not.
    - -o_dir, --output_directory: Path to the directory where the output files will be saved. By default, outputs will be saved at ``$CWAS_WORKSPACE``.
    - -s, --sample_info: Path to the txt file containing the sample information for each sample. This file must have three columns (``SAMPLE``, ``FAMILY``, ``PHENOTYPE``) with the exact name.
    - -a, --adjustment_factor: Path to the txt file containing the adjust factors for each sample. This is optional. With this option, CWAS-Plus multiplies the number of variants (or carriers, in -u option) with the adjust factor per sample.
    - -n, --num_perm: Number of permutations for label-swapping. By default, 10000.
    - -p, --num_proc: Number of worker processes that will be used for the permutation process. By default, 1.
    - -b, --burden_shift: Generates an output file containing binomial p-values for each label-swapped permutation. By default, False.
    - -rr, --perm_rr: Generates an output file containing relative risks for each label-swapped permutation. By default, False.
    - -u, --use_n_carrier: Enables the use of the number of samples with variants in each category for burden test instead of the number of variants. With this option, CWAS-Plus counts the number of samples that carry at least one variant of each category.

     .. code-block:: solidity

        cwas permutation_test -i INPUT.categorization_result.txt.gz -o_dir OUTPUT_DIR -s SAMPLE_LIST.txt -a ADJUST_FACTOR.txt -n 10000 -p 8 -b

  The specific descriptions of the output files are as below. Each output file containing a specific pattern (i.e., ``.burden_test.txt.gz``, ``.permutation_test.txt.gz``, ``.binom_pvals.txt.gz``) in the file name as below will be found in the output directory.

  - OUTPUT.burden_test.txt.gz: The final output file containing relative risk, two-sided binomial p-value and one-sided binomial p-value of each category.
  - OUTPUT.permutation_test.txt.gz: The final output file containing p-values calculated from permutations. This file will be used for :ref:`DAWN analysis <dawn>`.
  - OUTPUT.binom_pvals.txt.gz: The matrix containing binomial p-values generated from each permutation. This file will be generated only with ``-b`` option given.



7.  :ref:`Calculate the number of effective tests <effnumtest>`
####################################################################

  From correlation matrix, compute eigen values and vectors. Based on these outputs, users can calculate the number of effective tests.
  
  The parameters of the command are as below:

    - -i, --input_file: Path to the concatenated z-scores.
    - -if, --input_format: Specify the format of the input file. Available options are ``corr`` or ``inter``. By default, ``corr`` will be used. Each format refers to the following:

      - corr: A matrix with correlation values between categories.
      - inter: A matrix with intersected number of variants (or samples) between categories.

    - -o_dir, --output_directory: Path to the directory where the output files will be saved. By default, outputs will be saved at ``$CWAS_WORKSPACE``.
    - -n, --num_sim: Number of eigen values to use in calculating the number of effective tests. The maximum number is equivalent to the number of categories. By default, 10000.
    - -s, --sample_info: Path to the txt file containing the sample information for each sample. This file must have three columns (``SAMPLE``, ``FAMILY``, ``PHENOTYPE``) with the exact name. Required only when input format is set to ``inter``. By default, None.
    - -t, --tag: Tag used for the name of the output files. By default, None.
    - -c, --category_set_path: Path to a text file containing categories for eigen decomposition. If not specified, all of the categories in the z-score file will be used. This file must contain ``Category`` column with the name of categories to be used.
    - -ef, --eff_num_test: Calculate the effective number of tests. By default, False.

    +-------------------------------------------------------+
    |Category                                               |
    +=======================================================+
    |All_CHD8Common_All_IntergenicRegion_EarlyCREMicro      |
    +-------------------------------------------------------+
    |All_CHD8Common_phastCons46way_PromoterRegion_EarlyCREL4|
    +-------------------------------------------------------+
    |All_DDD_All_PromoterRegion_EarlyCREOligo               |
    +-------------------------------------------------------+


  .. code-block:: solidity

    cwas effective_num_test -i INPUT.correlation_matrix.pkl -o_dir OUTPUT_DIR -t test -ef -if corr -n 7918 -c CATEGORY_SET.txt


  The specific descriptions of the output files are as below. Each output file containing a specific pattern (i.e., ``.neg_lap.*.pickle``, ``.eig_vals.*.pickle``, ``.eig_vecs.*.txt.gz``) in the file name as below will be found in the output directory. If users set tag, the tag will be inserted in the file name like this: ``OUTPUT.eig_vecs.tag.txt.gz``.

  - OUTPUT.neg_lap.pickle: The negative laplacian matrix. This file is an intermediate output during eigen decomposition.
  - OUTPUT.eig_vals.pickle: The matrix containing eigen values. This file will be used to calculate the number of effective tests.
  - OUTPUT.eig_vecs.txt.gz: The matrix containing eigen vectors. This file will be used as an input for :ref:`DAWN analysis <dawn>`.

  In addition, the number of effective tests will be printed as below when ``-ef`` option is given. The number will also be written in ``.cwas_env`` as environment variable ``N_EFFECTIVE_TEST``.

  .. code-block:: solidity
    
    [RESULT] The number of effective tests is 1438.



8.  :ref:`Risk score analysis <riskscore>`
############################################

  Identify the best predictor of the phenotype by training Lasso regression model with the number of variants within each category across samples.
  
  CWAS-Plus utilizes categorized results to estimate the optimal predictor for the phenotype. It trains a Lasso regression model using the number of variants within each category across samples. After training the model with a subset of samples, the remaining test set is employed to calculate the |R2|. The significance of the |R2| value is determined by calculating it from samples with a randomly shuffled phenotype. The number of regressions (-n_reg) can be set to obtain the average |R2| value from all regressions.

  .. |R2| replace:: R\ :sup:`2`

  The parameters of the command are as below:
  
  - -i, --input_file: Path to the categorized txt file, resulted from categorization process. This file could be gzipped or not.
  - -o_dir, --output_directory: Path to the directory where the output files will be saved. By default, outputs will be saved at ``$CWAS_WORKSPACE``.
  - -s, --sample_info: Path to the txt file containing the sample information for each sample. This file must have three columns (``SAMPLE``, ``FAMILY``, ``PHENOTYPE``) with the exact name.
  - -a, --adjustment_factor: Path to the txt file containing the adjust factors for each sample. This is optional. With this option, CWAS-Plus multiplies the number of variants (or carriers, in -u option) with the adjust factor per sample.
  - -c, --category_set_path: Path to a text file containing categories for training. If not specified, all of the categories categorization file will be used. This file must contain ``Category`` column with the name of categories to be used.
  - -t, --tag: Tag used for the name of the output files. By default, None.
  - -u, --use_n_carrier: Enables the use of the number of samples with variants in each category for burden test instead of the number of variants. With this option, CWAS-Plus counts the number of samples that carry at least one variant of each category.
  - -thr, --threshold: The number of variants in controls (or the number of control carriers) used to select rare categories. For example, if set to 3, categories with less than 3 variants in controls will be used for training. By default, 3.
  - -tf, --train_set_fraction: The fraction of the training set. For example, if set to 0.7, 70% of the samples will be used as training set and 30% will be used as test set. By default, 0.7.
  - -n_reg, --num_regression: Number of regression trials to calculate a mean of R squares. By default, 10.
  - -f, --fold: Number of folds for cross-validation.
  - -l, --logistic: (hold) Make a logistic model with L1 penalty. By default, False.
  - -n, --n_permute: The number of permutations used to calculate the p-value. By default, 1,000.
  - --predict_only: If set, only predict the risk score and skip the permutation process. By default, False.
  - -p, --num_proc: Number of worker processes that will be used for the permutation process. By default, 1.


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


  The specific descriptions of the output files are as below. Each output file containing a specific pattern (i.e., ``.lasso_results_thres_*.txt``, ``.lasso_null_models_thres_*.txt``, ``.lasso_histogram_thres_*.pdf``, ``lasso_coef_thres_*.txt``) in the file name as below will be found in the output directory. If users set tag, the tag will be inserted in the file name like this: ``OUTPUT.eig_vecs.tag.txt.gz``.

  - OUTPUT.lasso_results_thres_*.txt: 
  - OUTPUT.lasso_null_models_thres_*.txt: 
  - OUTPUT.lasso_histogram_thres_*.pdf: Histogram plot for the observed predictive |R2| and random distribution. The random distribution is obtained from samples with a randomly shuffled phenotype. The x axis refers to the observed |R2| and the y axis refers to the frequency of |R2| s.
  - OUTPUT.lasso_coef_thres_*.txt: 


9.  :ref:`Burden shift analysis <burdenshift>`
################################################

  Identify the overrepresented domains associated to the phenotype.

  *In progress.*

  .. code-block:: solidity


10.  :ref:`DAWN analysis <dawn>`
##################################

  Investigate the relationship between categories and identify the specific type of categories clustered within the network.

  The parameters of the command are as below:

  - -i_dir, --input_directory: Path to the directory where the input files are stored. This directory must include three required files.

    - Eigen vector file: This is the output file from :ref:`calculation of effective number of tests <effnumtest>`. The file name must have pattern ``*eig_vecs*.txt.gz``.
    - Category correlation matrix file: This is the output file from :ref:`categorization <categorization>`. The file name must have pattern ``*correlation_matrix*.pkl``.
    - Permutation test file: This is the output file from :ref:`burden test <permtest>`. The file name must have pattern ``*permutation_test*.txt.gz``.

  - -o_dir, --output_directory: Path to the directory where the output files will be saved. By default, outputs will be saved at ``$CWAS_WORKSPACE``.
  - -r, --range: Range (i.e., (start,end)) to find optimal K for k-means clustering. It must contain two integers that are comma-separated. The first integer refers to the start number and must be above 1. The second integer refers to the end.
  - -k, --k_val: K for K-means clustering. With this argument, users can determine K manually. ``-r`` and ``-k`` arguments are mutually exclusive. If ``-k`` is given, ``-r`` will be ignored.
  - -s, --seed: Seed value for t-SNE. Same seed will generate same results for the same inputs.
  - -t, --tag: Tag used for the name of the output files. By default, None.
  - -c, --category_set_path: Path to a text file containing categories for training. If not specified, all of the categories categorization file will be used. This file must contain ``Category`` column with the name of categories to be used.
  - -c_count, --cat_count
  - -CT, --count_threshold: The treshold of variant (or sample) counts. The least amount of variants a category should have.
  - -CR, --corr_threshold: The threshold of correlation values between clusters. Computed by the mean value of correlation values of categories within a cluster.
  - -S, --size_threshold: The threshold of the number of categories per cluster. The least amount of categories a cluster should have.
  - -p, --num_proc: Number of worker processes that will be used for the DAWN analysis. By default, 1.


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


