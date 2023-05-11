.. _burdentest:

###############################
Burden test
###############################

With categorized results, CWAS-Plus calculate the burden of each category by comparing the rate of variants per cases and the rate of variants per controls.

For burden measurement, the package uses relative risk (RR), which is calculated by comparing the number of variants per phenotype group (RR>1, case burden; RR<1, control burden). The burden test in CWAS-Plus contains two types of p-value computation methods, binomial test and permutation test, to find more accurate p statistics.

--------------------------------
Binomial test
--------------------------------

For binomial tests, two types of tests are used: two-sided binomial tests and one-sided binomial tests with an alternative hypothesis of "greater".

- -i, --input_file : Path to the categorized txt file, resulted from categorization process. This file could be gzipped or not.
- -o_dir, --output_directory : Path to the directory where the output files will be saved. By default, outputs will be saved at ``$CWAS_WORKSPACE``.
- -s, --sample_info : Path to the txt file containing the sample information for each sample. This file must have three columns (``SAMPLE``, ``FAMILY``, ``PHENOTYPE``) with the exact name.

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

- -a, --adjustment_factor : Path to the txt file containing the adjust factors for each sample. This is optional. With this option, CWAS-Plus multiplies the number of variants (or carriers, in -u option) with the adjust factor per sample.

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

- -u, --use_n_carrier : Enables the use of the number of samples with variants in each category for burden test instead of the number of variants. With this option, CWAS-Plus counts the number of samples that carry at least one variant of each category.

.. code-block:: solidity

    cwas binomial_test -i INPUT.categorization_result.txt.gz -o_dir OUTPUT_DIR -s SAMPLE_LIST.txt -a ADJUST_FACTOR.txt

--------------------------------
Permutation test
--------------------------------

In permutation tests, phenotype labels (case, control) are randomly swapped by randomly selecting cases. Here, the number of selected cases are same as the original number of cases. After swapping, relative risks of each category are calculated. This swap processes are conducted 10,000 times by default (set the number of permutations by ``-n`` option).

Permutation p-values are calculated by comparing the relative risk of the original results with the relative risk from swapped results. The direction of burden is also considered. For example, when permutated 10,000 times, category **A** with relative risk=1.3 and p-value=0.01 refers to the following:

- Out of the 10,000 relative risks of category **A** obtained from random swapping, 100 were larger than the original relative risk.

In another example, when permutated 10,000 times, category **B** with relative risk=0.7 and p-value=0.01 refers to the following:

- Out of the 10,000 relative risks of category **B** obtained from random swapping, 100 were smaller than the original relative risk.

The parameters of the command are as below:

- -i, --input_file : Path to the categorized txt file, resulted from categorization process. This file could be gzipped or not.
- -o_dir, --output_directory : Path to the directory where the output files will be saved. By default, outputs will be saved at ``$CWAS_WORKSPACE``.
- -s, --sample_info : Path to the txt file containing the sample information for each sample. This file must have three columns (``SAMPLE``, ``FAMILY``, ``PHENOTYPE``) with the exact name.
- -a, --adjustment_factor : Path to the txt file containing the adjust factors for each sample. This is optional. With this option, CWAS-Plus multiplies the number of variants (or carriers, in -u option) with the adjust factor per sample.
- -n, --num_perm : Number of permutations for label-swapping. By default, 10000.
- -p, --num_proc : Number of worker processes that will be used for the permutation process. By default, 1.
- -b, --burden_shift : Generates an output file containing binomial p-values for each label-swapped permutation. By default, False.
- -rr, --perm_rr : Generates an output file containing relative risks for each label-swapped permutation. By default, False.
- -u, --use_n_carrier : Enables the use of the number of samples with variants in each category for burden test instead of the number of variants. With this option, CWAS-Plus counts the number of samples that carry at least one variant of each category.

.. code-block:: solidity

    cwas permutation_test -i INPUT.categorization_result.txt.gz -o_dir OUTPUT_DIR -s SAMPLE_LIST.txt -a ADJUST_FACTOR.txt -n 10000 -p 8 -b



