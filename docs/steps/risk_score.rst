.. _riskscore:

*************************
Risk score analysis
*************************

CWAS-Plus utilizes categorized results to estimate the optimal predictor for the phenotype. It trains a Lasso regression model using the number of variants within each category across samples. After training the model with a subset of samples, the remaining test set is employed to calculate the |R2|. The significance of the |R2| value is determined by calculating it from samples with a randomly shuffled phenotype. The number of regressions (-n_reg) can be set to obtain the average |R2| value from all regressions.

.. |R2| replace:: R\ :sup:`2`


- -i, --input_file: Path to the categorized zarr directory, resulted from categorization process.
- -o_dir, --output_directory: Path to the directory where the output files will be saved. By default, outputs will be saved at ``$CWAS_WORKSPACE``.
- -s, --sample_info: Path to the txt file containing the sample information for each sample. This file must have three columns (``SAMPLE``, ``FAMILY``, ``PHENOTYPE``) with the exact name.

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

- -a, --adjustment_factor: Path to the txt file containing the adjust factors for each sample. This is optional. With this option, CWAS-Plus multiplies the number of variants (or carriers, in -u option) with the adjust factor per sample.

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

- -c_info, --category_info: Path to a text file category information (`*.category_info.txt`).
- -d, --domain_list: Domain list to filter categories based on GENCODE domain. If 'run_all' is given, all available options will be tested. Available options are `run_all,all,coding,noncoding,ptv,missense,damaging_missense,promoter,noncoding_wo_promoter,intron,intergenic,utr,lincRNA`. By default, all.
- -t, --tag: Tag used for the name of the output files. By default, None.
- --do_each_one: Use each annotation from functional annotation to calculate risk score. By default, False.
- --leave_one_out: Calculate risk score while excluding one annotation from functional annotation. This option is not used when the `--do_each_one` flag is enabled. By default, False.
- -u, --use_n_carrier: Enables the sample-level analysis (the use of the number of samples with variants in each category for burden test instead of the number of variants). With this option, CWAS-Plus counts the number of samples that carry at least one variant of each category.
- -thr, --threshold: The number of variants in controls (or the number of control carriers) used to select rare categories. For example, if set to 3, categories with less than 3 variants in controls will be used for training. By default, 3.
- -tf, --train_set_fraction: The fraction of the training set. For example, if set to 0.7, 70% of the samples will be used as training set and 30% will be used as test set. By default, 0.7.
- -n_reg, --num_regression: Number of regression trials to calculate a mean of R squares. By default, 10.
- -f, --fold: Number of folds for cross-validation.
- -n, --n_permute: The number of permutations used to calculate the p-value. By default, 1,000.
- --predict_only: If set, only predict the risk score and skip the permutation process. By default, False.
- -S, --seed: Seed of random state. By default, 42.
- -p, --num_proc: Number of worker processes that will be used for the permutation process. By default, 1.

.. code-block:: solidity
  
  cwas risk_score -i INPUT.categorization_result.txt.gz \
  -o_dir OUTPUT_DIR \
  -s SAMPLE_LIST.txt \
  -a ADJUST_FACTOR.txt \
  -c_info CATEGORY_SET.txt \
  -thr 3 \
  -tf 0.7 \
  -n_reg 10 \
  -f 5 \
  -n 1000 \
  -p 8


Users can perform two types of risk score analyses in a loop to identify annotations with the best predictive performance and composition within the annotation set.

1. Risk score analysis for categories containing a single annotation within a specific domain

    .. code-block:: solidity
    
        cwas risk_score -i INPUT.categorization_result.txt.gz \
        -o_dir OUTPUT_DIR \
        -s SAMPLE_LIST.txt \
        -a ADJUST_FACTOR.txt \
        -c_info CATEGORY_SET.txt \
        -thr 3 \
        -tf 0.7 \
        -n_reg 10 \
        -f 5 \
        -n 1000 \
        -p 8 \
        --do_each_one

2. Risk score analysis for categories with one annotation excluded from the total annotations

    .. code-block:: solidity
        
        cwas risk_score -i INPUT.categorization_result.txt.gz \
        -o_dir OUTPUT_DIR \
        -s SAMPLE_LIST.txt \
        -a ADJUST_FACTOR.txt \
        -c_info CATEGORY_SET.txt \
        -thr 3 \
        -tf 0.7 \
        -n_reg 10 \
        -f 5 \
        -n 1000 \
        -p 8 \
        --leave_one_out


