#!/usr/bin/env python
"""
Script to run burden tests for de novo variants (DNVs) of cases and controls in each CWAS category.
Binomial tests and permutation tests are available in this script.

For more detailed information, please refer to An et al., 2018 (PMID 30545852).

"""
import argparse
import multiprocessing as mp
import os
import re
import sys
from datetime import datetime
from functools import partial

import numpy as np
import pandas as pd
from scipy.stats import binom_test


def main(cat_result_path, adj_file_path, outfile_path):
    # Print the description and run settings
    print(__doc__)
    print(f'[Setting] The input CWAS categorization result: {cat_result_path}')
    print(f'[Setting] The list of adjustment factors for No. DNVs of each sample: {adj_file_path}')
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

    # Run burden tests
    print(f'[{get_curr_time()}, Progress] Burden tests via binomial tests')
    burden_df = run_burden_binom(cwas_cat_df)

    # Write results of the burden tests
    print(f'[{get_curr_time()}, Progress] Write the result of the burden tests')
    burden_df.to_csv(outfile_path, sep='\t')

    print(f'[{get_curr_time()}, Progress] Done')


def check_sample_id_format(sample_ids: np.ndarray):
    """ Function to check whether formats of the sample IDs are matched with the asserted format """
    sample_id_f = re.compile(r'^\d+[.]([ps])\d$')

    for sample_id in sample_ids:
        assert sample_id_f.match(sample_id), \
            f'Wrong sample ID format "{sample_id}". The regular expression of a sample ID is "^\d+[.]([ps])\d$".'


def run_burden_binom(cwas_cat_df: pd.DataFrame) -> pd.DataFrame:
    """ Function for burden tests via binomial tests

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

    # Make a DataFrame for the results of binomial tests
    burden_df = \
        pd.DataFrame(dnv_cnt_arr, index=cwas_cat_df.columns.values, columns=['Case_DNV_Count', 'Ctrl_DNV_Count'])
    burden_df.index.name = 'Category'
    burden_df['Relative_Risk'] = burden_df['Case_DNV_Count'].values / burden_df['Ctrl_DNV_Count'].values

    # Binomial tests
    binom_two_tail = lambda n1, n2: binom_test(x=n1, n=n1+n2, p=0.5, alternative='two-sided')
    binom_one_tail = lambda n1, n2: binom_test(x=n1, n=n1+n2, p=0.5, alternative='greater') if n1 > n2 \
        else binom_test(x=n2, n=n1+n2, p=0.5, alternative='greater')
    burden_df['P_Two-tail'] = \
        np.vectorize(binom_two_tail)(burden_df['Case_DNV_Count'].values, burden_df['Ctrl_DNV_Count'].values)
    burden_df['P_One-tail'] = \
        np.vectorize(binom_one_tail)(burden_df['Case_DNV_Count'].values, burden_df['Ctrl_DNV_Count'].values)

    return burden_df


def run_burden_perm(cwas_cat_df: pd.DataFrame, num_perm: int, num_proc: int) -> (pd.DataFrame, pd.DataFrame):
    """ Function for burden tests via permutation tests

    :param cwas_cat_df: A DataFrame that contains the result of CWAS categorization
    :param num_perm: The number of permutation trials
    :param num_proc: The number of processes used in this function (for multiprocessing)
    :returns:
        1. A DataFrame that contains permutation p-values and other statistics for each CWAS category
        2. A DataFrame that contains relative risks for each category after each permutation trial
    """
    # Calculate relative risks after permutation
    if num_proc == 1:
        perm_rr_list = get_perm_rr(num_perm, cwas_cat_df)
    else:
        num_perms = div_dist_num(num_perm, num_proc)
        pool = mp.Pool(num_proc)
        proc_outputs = pool.map(partial(get_perm_rr, cwas_cat_df=cwas_cat_df), num_perms)
        pool.close()
        pool.join()

        perm_rr_list = []

        for proc_output in proc_outputs:
            perm_rr_list += proc_output

    perm_rrs = np.concatenate(perm_rr_list, axis=0)

    # Calculate an original relative risk
    is_prob_func = lambda sample_id: 'p' in sample_id
    sample_ids = cwas_cat_df.index.values
    cwas_cats = cwas_cat_df.columns.values
    cat_df_vals = cwas_cat_df.values

    are_prob = np.vectorize(is_prob_func)(sample_ids)
    prob_dnv_cnt = cat_df_vals[are_prob, :].sum(axis=0)
    sib_dnv_cnt = cat_df_vals[~are_prob, :].sum(axis=0)
    rr = prob_dnv_cnt / sib_dnv_cnt

    # Permutation tests
    # Check whether a permutation RR is more extreme than the original RR
    are_ext_rr = ((rr >= 1) & (perm_rrs >= rr)) | ((rr < 1) & (perm_rrs <= rr))
    ext_rr_cnt = are_ext_rr.sum(axis=0)
    perm_p = ext_rr_cnt / num_perm

    # Return the results as a DataFrame
    burden_df = pd.concat([
        pd.Series(prob_dnv_cnt, index=cwas_cats, name='Case_DNV_Count'),
        pd.Series(sib_dnv_cnt, index=cwas_cats, name='Ctrl_DNV_Count'),
        pd.Series(rr, index=cwas_cats, name='Relative_Risk'),
        pd.Series(perm_p, index=cwas_cats, name='P_Perm'),
    ], axis=1)
    burden_df.index.name = 'Category'

    perm_rr_df = pd.DataFrame(perm_rrs, columns=cwas_cats)
    perm_rr_df.index += 1
    perm_rr_df.index.name = 'Trial'

    return burden_df, perm_rr_df


def get_perm_rr(num_perm: int, cwas_cat_df: pd.DataFrame) -> list:
    """ Calculate relative risks of each category after each permutation trial.
    The length of the returned list equals to the number of permutations.
    """
    is_prob_func = lambda sample_id: 'p' in sample_id
    sample_ids = cwas_cat_df.index.values
    cat_df_vals = cwas_cat_df.values
    perm_rr_list = []

    for _ in range(num_perm):
        swap_sample_ids = swap_label(sample_ids)
        are_prob = np.vectorize(is_prob_func)(swap_sample_ids)
        prob_dnv_cnt = cat_df_vals[are_prob, :].sum(axis=0)
        sib_dnv_cnt = cat_df_vals[~are_prob, :].sum(axis=0)
        perm_rr = prob_dnv_cnt / sib_dnv_cnt
        perm_rr_list.append(perm_rr[np.newaxis, :])

    return perm_rr_list


def swap_label(sample_ids: np.ndarray) -> np.ndarray:
    """ Randomly swap labels (proband and sibling) of each family in the samples
    and return a list of the sample IDs after swapped.
    """
    fam_to_hit_cnt = {}  # Key: A family, Value: The number of times the same family is referred
    fam_to_idx = {}  # Key: A family, Value: The index of a sample ID firstly matched with the family
    swap_sample_ids = np.copy(sample_ids)

    for i, sample_id in enumerate(sample_ids):
        fam = sample_id.split('.')[0]

        if fam_to_hit_cnt.get(fam, 0) == 0:
            fam_to_hit_cnt[fam] = 1
            fam_to_idx[fam] = i
        elif fam_to_hit_cnt.get(fam, 0) == 1:
            fam_to_hit_cnt[fam] += 1
            prev_idx = fam_to_idx[fam]
            # Swap
            do_swap = np.random.binomial(1, 0.5)
            if do_swap:
                swap_sample_ids[i] = sample_ids[prev_idx]
                swap_sample_ids[prev_idx] = sample_id
        else:
            print(f'[ERROR] The fam {fam} has more than 2 samples.', file=sys.stderr)
            raise AssertionError('One family must have two samples (one proband and one sibling).')

    return swap_sample_ids


def div_dist_num(num: int, num_group: int) -> list:
    """ Divide and distribute the number to each group almost equally. """
    num_per_groups = []
    num_per_group = num // num_group
    remain_num = num % num_group

    for _ in range(num_group):
        if remain_num == 0:
            num_per_groups.append(num_per_group)
        else:
            num_per_groups.append(num_per_group + 1)
            remain_num -= 1

    return num_per_groups


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
