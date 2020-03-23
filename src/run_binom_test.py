#!/usr/bin/env python
__version__ = 0.2
__author__ = 'Joon An'
__date__ = 'October 5th, 2018'

description = '''
				Script for burden test.
			'''

import argparse

import numpy as np
import pandas as pd
from scipy.stats import binom_test


def main(cat_result_path, adj_file_path, outfile_path):
    if adj_file_path is None:
        cwas_cat_df = pd.read_table(cat_result_path, index_col='SampleID')
        burden_df = run_burden_binom(cwas_cat_df)
        burden_df.to_csv(outfile_path, sep='\t')
    else:
        # Load the input data
        cwas_cat_df = pd.read_table(cat_result_path, index_col='SampleID')
        adj_factor_df = pd.read_table(adj_file_path)

        # Match the format of sample IDs in the adjustment file with that of the previous categorization result
        adj_factor_df['SampleID'] = np.vectorize(lambda x: x.replace('_', '.'))(adj_factor_df['SampleID'].values)
        are_same_samples = cwas_cat_df.index.values == adj_factor_df['SampleID'].values
        assert np.all(are_same_samples)

        # Adjust the rate of de novo variants by adjustment factors for every sample
        adj_cwas_cat_df = cwas_cat_df.multiply(adj_factor_df['AdjustFactor'].values, axis='index')
        adj_cwas_cat_df.index.name = 'SampleID'
        adj_cwas_cat_df = adj_cwas_cat_df.astype('int64')

        # Run burden analysis using binomial test
        burden_df = run_burden_binom(cwas_cat_df)
        adj_burden_df = run_burden_binom(adj_cwas_cat_df)

        # Write the result
        adj_burden_df = adj_burden_df.rename(columns=lambda colname: 'Adj_' + colname)
        output_df = pd.concat([burden_df, adj_burden_df], axis='columns')
        output_df.to_csv(outfile_path, sep='\t')


def run_burden_binom(cwas_cat_df):
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
    burden_df['Binom_P'] = \
        np.vectorize(binom_two_tail)(burden_df['Case_DNV_Count'].values, burden_df['Ctrl_DNV_Count'].values)
    burden_df['Binom_P_1tail'] = \
        np.vectorize(binom_one_tail)(burden_df['Case_DNV_Count'].values, burden_df['Ctrl_DNV_Count'].values)

    return burden_df


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('-i', '--infile', dest='cat_result_path', required=True, type=str,
                        help='The path of a result of the CWAS categorization')
    parser.add_argument('-a', '--adj', dest='adj_file_path', required=False, type=str,
                        help='The file that contains adjustment factors for No. DNVs of each sample',
                        default=None)
    parser.add_argument('-o', '--outfile', dest='outfile_path', required=False, type=str,
                        help='The path of a result of this script', default='cwas_burden_binom_result.txt')
    args = parser.parse_args()
    main(args.cat_result_path, args.adj_file_path, args.outfile_path)
