#!/usr/bin/env python
"""
Script for mutation burden analysis for each CWAS category using binomial test

For more detailed information, please refer to An et al., 2018 (PMID 30545852).

"""
import argparse
import os
import re
from datetime import datetime

import numpy as np
import pandas as pd
from scipy.stats import binom_test


def main(cat_result_path, adj_file_path, outfile_path):
    # Print the description and run settings
    print(__doc__)
    print(f'[Setting] The input CWAS categorization result: {cat_result_path}')
    print(f'[Setting] The file with adjustment factors of No. DNVs for each sample: {adj_file_path}')
    print(f'[Setting] The output path: {outfile_path}')
    print()

    # Check the validity of the settings
    assert os.path.isfile(cat_result_path), f'The input file "{cat_result_path}" cannot be found.'
    assert adj_file_path is None or os.path.isfile(adj_file_path), f'The input file "{adj_file_path}" cannot be found.'
    outfile_dir = os.path.dirname(outfile_path)
    assert outfile_dir == '' or os.path.isdir(outfile_dir), f'The outfile directory "{outfile_dir}" cannot be found.'

    # Load and parse the input data
    print(f'[{get_curr_time()}, Progress] Load the categorization result into DataFrame')
    cwas_cat_df = pd.read_table(cat_result_path, index_col='SampleID')
    check_sample_id_format(cwas_cat_df.index.values)  # It can raise AssertionError.

    if adj_file_path is not None:
        print(f'[{get_curr_time()}, Progress] Adjust No. DNVs of each sample')
        adj_factor_df = pd.read_table(adj_file_path)

        # Match the format of sample IDs in the adjustment file with that of the previous categorization result
        adj_factor_df['SampleID'] = np.vectorize(lambda x: x.replace('_', '.'))(adj_factor_df['SampleID'].values)
        are_same_samples = cwas_cat_df.index.values == adj_factor_df['SampleID'].values
        assert np.all(are_same_samples), "The lists of sample IDs from the two input files are not consistent."

        # Adjust the number of de novo variants of each sample by adjustment factors
        cwas_cat_df = cwas_cat_df.multiply(adj_factor_df['AdjustFactor'].values, axis='index')
        cwas_cat_df.index.name = 'SampleID'
        cwas_cat_df = cwas_cat_df.astype('int64')

    # Run burden analysis
    print(f'[{get_curr_time()}, Progress] Burden analysis via binomial tests')
    burden_df = run_burden_binom(cwas_cat_df)

    # Write the burden analysis result
    print(f'[{get_curr_time()}, Progress] Write the result of the burden analysis')
    burden_df.to_csv(outfile_path, sep='\t')

    print(f'[{get_curr_time()}, Progress] Done')


def check_sample_id_format(sample_ids: np.ndarray):
    """ Function to check whether formats of the sample IDs are matched with the asserted format """
    sample_id_f = re.compile(r'^\d+[.]([ps])\d$')

    for sample_id in sample_ids:
        assert sample_id_f.match(sample_id), \
            f'Wrong sample ID format "{sample_id}". The regular expression of a sample ID is "^\d+[.]([ps])\d$".'


def run_burden_binom(cwas_cat_df: pd.DataFrame) -> pd.DataFrame:
    """ Function for the mutation burden analysis using binomial test

    :param cwas_cat_df: A DataFrame that contains the result of CWAS categorization
    :return: A DataFrame that contains binomial p-values and other statistics for each CWAS category
    """
    # Count the number of de novo variants (DNV) for probands and siblings
    is_prob_func = lambda sample_id: 'p' in sample_id
    sample_ids = cwas_cat_df.index.values
    are_prob = np.vectorize(is_prob_func)(sample_ids)
    cat_df_vals = cwas_cat_df.values
    prob_dnv_cnt = cat_df_vals[are_prob, :].sum(axis=0)
    sib_dnv_cnt = cat_df_vals[~are_prob, :].sum(axis=0)
    dnv_cnt_arr = np.concatenate([prob_dnv_cnt[:, np.newaxis], sib_dnv_cnt[:, np.newaxis]], axis=1)

    # Make a DataFrame for burden analysis results
    burden_df = \
        pd.DataFrame(dnv_cnt_arr, index=cwas_cat_df.columns.values, columns=['Case_DNV_Count', 'Ctrl_DNV_Count'])
    burden_df.index.name = 'Category'
    burden_df['Relative_Risk'] = burden_df['Case_DNV_Count'].values / burden_df['Ctrl_DNV_Count'].values

    # Burden analysis via binomial test
    binom_two_tail = lambda n1, n2: binom_test(x=n1, n=n1+n2, p=0.5, alternative='two-sided')
    binom_one_tail = lambda n1, n2: binom_test(x=n1, n=n1+n2, p=0.5, alternative='greater') if n1 > n2 \
        else binom_test(x=n2, n=n1+n2, p=0.5, alternative='greater')
    burden_df['P_Two-tail'] = \
        np.vectorize(binom_two_tail)(burden_df['Case_DNV_Count'].values, burden_df['Ctrl_DNV_Count'].values)
    burden_df['P_One-tail'] = \
        np.vectorize(binom_one_tail)(burden_df['Case_DNV_Count'].values, burden_df['Ctrl_DNV_Count'].values)

    return burden_df


def get_curr_time() -> str:
    now = datetime.now()
    curr_time = now.strftime('%H:%M:%S %m/%d/%y')
    return curr_time


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-i', '--infile', dest='cat_result_path', required=True, type=str,
                        help='The path of a result of the CWAS categorization')
    parser.add_argument('-a', '--adj', dest='adj_file_path', required=False, type=str,
                        help='The file that contains adjustment factors for No. DNVs of each sample',
                        default=None)
    parser.add_argument('-o', '--outfile', dest='outfile_path', required=False, type=str,
                        help='The path of a result of this script', default='cwas_burden_binom_result.txt')
    args = parser.parse_args()
    main(args.cat_result_path, args.adj_file_path, args.outfile_path)
