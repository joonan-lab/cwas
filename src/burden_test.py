#!/usr/bin/env python
"""
Script to run burden tests for de novo variants (DNVs) of cases and controls in each CWAS category.
Binomial tests and permutation tests are available in this script.

For more detailed information, please refer to An et al., 2018 (PMID 30545852).

"""
import argparse
import multiprocessing as mp
import os
import sys
from datetime import datetime
from functools import partial

import numpy as np
import pandas as pd
from scipy.stats import binom_test


def main():
    # Cteate the top-level argument parser
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(description='Types of burden tests', dest='test_type')

    # Create the parser for binomial tests
    parser_binom = subparsers.add_parser('binom', description='Burden tests via binomial tests',
                                         help='Binomial tests (arg "binom -h" for usage)')
    parser_binom.add_argument('-i', '--infile', dest='cat_result_path', required=True, type=str,
                              help='Path of a result of the CWAS categorization')
    parser_binom.add_argument('-s', '--sample_file', dest='sample_file_path', required=True, type=str,
                              help='File listing sample IDs with their families and phenotypes (case or ctrl)')
    parser_binom.add_argument('-a', '--adj_file', dest='adj_file_path', required=False, type=str,
                              help='File that contains adjustment factors for No. DNVs of each sample',
                              default='')
    parser_binom.add_argument('-o', '--outfile', dest='outfile_path', required=False, type=str,
                              help='Path of results of burden tests', default='cwas_burden_binom_result.txt')

    # Create the parser for permutation tests
    parser_perm = subparsers.add_parser('perm', description='Burden tests via permutation tests',
                                        help='Permutation tests (arg "perm -h" for usage)')
    parser_perm.add_argument('-i', '--infile', dest='cat_result_path', required=True, type=str,
                              help='Path of a result of the CWAS categorization')
    parser_perm.add_argument('-s', '--sample_file', dest='sample_file_path', required=True, type=str,
                              help='File listing sample IDs with their families and phenotypes (case or ctrl)')
    parser_perm.add_argument('-a', '--adj_file', dest='adj_file_path', required=False, type=str,
                              help='File that contains adjustment factors for No. DNVs of each sample',
                              default='')
    parser_perm.add_argument('-o', '--outfile', dest='outfile_path', required=False, type=str,
                              help='Path of results of burden tests', default='cwas_burden_perm_result.txt')
    parser_perm.add_argument('-n', '--num_perm', dest='num_perm', required=False, type=int,
                             help='Number of label-swapping permutations',
                             default=10000)
    parser_perm.add_argument('-p', '--num_proc', dest='num_proc', required=False, type=int,
                             help='Number of processes used in permutation tests',
                             default=1)
    parser_perm.add_argument('-po', '--perm_outfile', dest='perm_rr_path', required=False, type=str,
                             help='Path of relative risk (RR) outputs from permutations',
                             default='')

    # Parse the arguments
    args = parser.parse_args()

    # Print the description
    print(__doc__)

    # Print and check the validity of the settings
    print(f'[Setting] Types of burden tests: {"Binomial test" if args.test_type == "binom" else "Permutation test"}')
    print(f'[Setting] The input CWAS categorization result: {args.cat_result_path}')
    print(f'[Setting] The list of sample IDs: {args.sample_file_path}')
    print(f'[Setting] The list of adjustment factors for No. DNVs of each sample: '
          f'{args.adj_file_path if args.adj_file_path else "None"}')
    print(f'[Setting] The output path: {args.outfile_path}')
    assert os.path.isfile(args.cat_result_path), f'The input file "{args.cat_result_path}" cannot be found.'
    assert os.path.isfile(args.sample_file_path), f'The input file "{args.sample_file_path}" cannot be found.'
    assert args.adj_file_path == '' or os.path.isfile(args.adj_file_path), \
        f'The input file "{args.adj_file_path}" cannot be found.'
    outfile_dir = os.path.dirname(args.outfile_path)
    assert outfile_dir == '' or os.path.isdir(outfile_dir), f'The outfile directory "{outfile_dir}" cannot be found.'

    if args.test_type == 'perm':
        print(f'[Setting] No. label-swapping permutations: {args.num_perm:,d}')
        print(f'[Setting] No. processes for the tests: {args.num_proc:,d}')
        print(f'[Setting] The permutation RR output path: {args.perm_rr_path if args.perm_rr_path else "None"}')
        assert 1 <= args.num_proc <= mp.cpu_count(), \
            f'Invalid number of processes "{args.num_proc:,d}". It must be in the range [1, {mp.cpu_count()}].'
        assert args.num_perm >= args.num_proc, f'No. processes must be equal or less than No. permutations.'
        perm_rr_dir = os.path.dirname(args.perm_rr_path)
        assert perm_rr_dir == '' or os.path.isdir(perm_rr_dir), \
            f'The outfile directory for RRs from permutations "{perm_rr_dir}" cannot be found.'
    print()

    # Load and parse the input data
    print(f'[{get_curr_time()}, Progress] Parse the categorization result into DataFrame')
    cwas_cat_df = pd.read_table(args.cat_result_path, index_col='SAMPLE')

    # Load and parse the file listing sample IDs
    print(f'[{get_curr_time()}, Progress] Parse the file listing sample IDs')
    sample_df = pd.read_table(args.sample_file_path, index_col='SAMPLE')
    assert cmp_two_arr(cwas_cat_df.index.values, sample_df.index.values), \
        f'The samples IDs of the categorization result are not the same ' \
        f'with the sample IDs of the file listing samples.'

    # Adjust No. DNVs of each sample in the categorization result
    if args.adj_file_path:
        print(f'[{get_curr_time()}, Progress] Adjust No. DNVs of each sample in the categorization result')
        adj_factor_df = pd.read_table(args.adj_file_path, index_col='SAMPLE')
        assert cmp_two_arr(adj_factor_df.index.values, sample_df.index.values), \
            f'The samples IDs of the categorization result are not the same ' \
            f'with the sample IDs of the file listing samples.'
        cwas_cat_df = adjust_cat_df(cwas_cat_df, adj_factor_df)

    # Run burden tests
    if args.test_type == 'binom':
        print(f'[{get_curr_time()}, Progress] Run burden tests via binomial tests')
        burden_df = run_burden_binom(cwas_cat_df, sample_df)
    else:  # args.test_type == 'perm'
        print(f'[{get_curr_time()}, Progress] Run burden tests via permutation tests')
        burden_df, perm_rr_df = run_burden_perm(cwas_cat_df, sample_df, args.num_perm, args.num_proc)

        if args.perm_rr_path:
            print(f'[{get_curr_time()}, Progress] Write lists of relative risks from label-swapping permutations')
            perm_rr_df.to_csv(args.perm_rr_path, sep='\t')

    # Write results of the burden tests
    print(f'[{get_curr_time()}, Progress] Write the result of the burden tests')
    burden_df.to_csv(args.outfile_path, sep='\t')

    print(f'[{get_curr_time()}, Progress] Done')


def cmp_two_arr(array1: np.ndarray, array2: np.ndarray) -> bool:
    """ Return True if two arrays have the same items regardless of the order, else return False """
    if len(array1) != len(array2):
        return False

    array1_item_set = set(array1)

    for item in array2:
        if item not in array1_item_set:
            return False

    return True


def adjust_cat_df(cwas_cat_df: pd.DataFrame, adj_factor_df: pd.DataFrame) -> pd.DataFrame:
    """ Adjust No. DNVs of each sample in a DataFrame for the result of CWAS categorization
    by adjustment factors for each sample

    :param cwas_cat_df: A DataFrame that contains the CWAS categorization result
    :param adj_factor_df: A DataFrame that contains a list of adjustment factors for each sample
    :return: A DataFrame that contains the CWAS categorization result with adjusted No. DNVs
    """
    # Reorder the adjustment factors
    sample_ids = cwas_cat_df.index.values
    adj_factor_dict = adj_factor_df.to_dict()['AdjustFactor']
    adj_factors = [adj_factor_dict[sample_id] for sample_id in sample_ids]

    # Adjust the number of de novo variants of each sample by adjustment factors
    cwas_cat_df = cwas_cat_df.multiply(adj_factors, axis='index')
    cwas_cat_df.index.name = 'SampleID'

    return cwas_cat_df


def run_burden_binom(cwas_cat_df: pd.DataFrame, sample_df: pd.DataFrame) -> pd.DataFrame:
    """ Function for burden tests via binomial tests

    :param cwas_cat_df: A DataFrame that contains the result of CWAS categorization
    :param sample_df: A DataFrame listing sample IDs with their families and phenotypes
    :return: A DataFrame that contains binomial p-values and other statistics for each CWAS category
    """
    # Count the number of de novo variants (DNV) for cases and controls
    sample_info_dict = sample_df.to_dict()
    sample_ids = cwas_cat_df.index.values
    sample_types = np.asarray([sample_info_dict['PHENOTYPE'][sample_id] for sample_id in sample_ids])
    are_case = sample_types == 'case'

    cat_df_vals = cwas_cat_df.values
    case_dnv_cnt = cat_df_vals[are_case, :].sum(axis=0)
    ctrl_dnv_cnt = cat_df_vals[~are_case, :].sum(axis=0)
    dnv_cnt_arr = np.concatenate([case_dnv_cnt[:, np.newaxis], ctrl_dnv_cnt[:, np.newaxis]], axis=1)

    # Make a DataFrame for the results of binomial tests
    burden_df = \
        pd.DataFrame(dnv_cnt_arr, index=cwas_cat_df.columns.values, columns=['Case_DNV_Count', 'Ctrl_DNV_Count'])
    burden_df.index.name = 'Category'
    burden_df['Relative_Risk'] = case_dnv_cnt / ctrl_dnv_cnt

    # Binomial tests
    binom_two_tail = lambda n1, n2: binom_test(x=n1, n=n1+n2, p=0.5, alternative='two-sided')
    binom_one_tail = lambda n1, n2: binom_test(x=n1, n=n1+n2, p=0.5, alternative='greater') if n1 > n2 \
        else binom_test(x=n2, n=n1+n2, p=0.5, alternative='greater')
    burden_df['P'] = \
        np.vectorize(binom_two_tail)(case_dnv_cnt.round(), ctrl_dnv_cnt.round())
    burden_df['P_1side'] = \
        np.vectorize(binom_one_tail)(case_dnv_cnt.round(), ctrl_dnv_cnt.round())

    return burden_df


def run_burden_perm(cwas_cat_df: pd.DataFrame, sample_df: pd.DataFrame, num_perm: int, num_proc: int) \
        -> (pd.DataFrame, pd.DataFrame):
    """ Function for burden tests via permutation tests

    :param cwas_cat_df: A DataFrame that contains the result of CWAS categorization
    :param sample_df: A DataFrame listing sample IDs with their families and phenotypes
    :param num_perm: The number of label-swapping permutation trials
    :param num_proc: The number of processes used in this function (for multiprocessing)
    :returns:
        1. A DataFrame that contains permutation p-values and other statistics for each CWAS category
        2. A DataFrame that contains relative risks for each category from each permutation trial
    """
    # Calculate relative risks from label-swapping permutations
    if num_proc == 1:
        perm_rr_list = cal_perm_rr(num_perm, cwas_cat_df, sample_df)
    else:
        num_perms = div_dist_num(num_perm, num_proc)
        pool = mp.Pool(num_proc)
        proc_outputs = pool.map(partial(cal_perm_rr, cwas_cat_df=cwas_cat_df, sample_df=sample_df), num_perms)
        pool.close()
        pool.join()

        perm_rr_list = []

        for proc_output in proc_outputs:
            perm_rr_list += proc_output

    perm_rrs = np.concatenate(perm_rr_list, axis=0)

    # Calculate an original relative risk
    sample_info_dict = sample_df.to_dict()
    sample_ids = cwas_cat_df.index.values
    sample_types = np.asarray([sample_info_dict['PHENOTYPE'][sample_id] for sample_id in sample_ids])
    are_case = sample_types == 'case'

    cat_df_vals = cwas_cat_df.values
    case_dnv_cnt = cat_df_vals[are_case, :].sum(axis=0)
    ctrl_dnv_cnt = cat_df_vals[~are_case, :].sum(axis=0)
    rr = case_dnv_cnt / ctrl_dnv_cnt

    # Permutation tests
    # Check whether a permutation RR is more extreme than the original RR
    are_ext_rr = ((rr >= 1) & (perm_rrs >= rr)) | ((rr < 1) & (perm_rrs <= rr))
    ext_rr_cnt = are_ext_rr.sum(axis=0)
    perm_p = ext_rr_cnt / num_perm

    # Return the results as a DataFrame
    cwas_cats = cwas_cat_df.columns.values

    burden_df = pd.concat([
        pd.Series(case_dnv_cnt, index=cwas_cats, name='Case_DNV_Count'),
        pd.Series(ctrl_dnv_cnt, index=cwas_cats, name='Ctrl_DNV_Count'),
        pd.Series(rr, index=cwas_cats, name='Relative_Risk'),
        pd.Series(perm_p, index=cwas_cats, name='P'),
    ], axis=1)
    burden_df.index.name = 'Category'

    perm_rr_df = pd.DataFrame(perm_rrs, columns=cwas_cats)
    perm_rr_df.index += 1
    perm_rr_df.index.name = 'Trial'

    return burden_df, perm_rr_df


def cal_perm_rr(num_perm: int, cwas_cat_df: pd.DataFrame, sample_df: pd.DataFrame) -> list:
    """ Calculate relative risks of each category in each permutation trial.
    The length of the returned list equals to the number of the permutations.
    """
    sample_info_dict = sample_df.to_dict()
    sample_ids = cwas_cat_df.index.values
    family_ids = np.asarray([sample_info_dict['FAMILY'][sample_id] for sample_id in sample_ids])
    sample_types = np.asarray([sample_info_dict['PHENOTYPE'][sample_id] for sample_id in sample_ids])
    cat_df_vals = cwas_cat_df.values
    perm_rr_list = []

    for _ in range(num_perm):
        swap_sample_types = swap_label(sample_types, family_ids)
        are_case = swap_sample_types == 'case'
        case_dnv_cnt = cat_df_vals[are_case, :].sum(axis=0)
        ctrl_dnv_cnt = cat_df_vals[~are_case, :].sum(axis=0)
        perm_rr = case_dnv_cnt / ctrl_dnv_cnt
        perm_rr_list.append(perm_rr[np.newaxis, :])

    return perm_rr_list


def swap_label(labels: np.ndarray, group_ids: np.ndarray) -> np.ndarray:
    """ Randomly swap labels (case or control) in each group and return a list of swapped labels.

    :param labels: Array of labels
    :param group_ids: Array of group IDs corresponding to each label
    :return: Swapped labels
    """
    group_to_hit_cnt = {group_id: 0 for group_id in group_ids}  # Key: A group, Value: No. times the key is referred
    group_to_idx = {}  # Key: A group, Value: The index of a label firstly matched with the group
    swap_labels = np.copy(labels)

    # Make an array for random swapping
    num_group = len(group_to_hit_cnt.keys())
    do_swaps = np.random.binomial(1, 0.5, size=num_group)
    group_idx = 0

    for i, label in enumerate(labels):
        group_id = group_ids[i]

        if group_to_hit_cnt.get(group_id, 0) == 0:
            group_to_hit_cnt[group_id] = 1
            group_to_idx[group_id] = i
        elif group_to_hit_cnt.get(group_id, 0) == 1:
            group_to_hit_cnt[group_id] += 1
            prev_idx = group_to_idx[group_id]

            # Swap
            do_swap = do_swaps[group_idx]
            if do_swap:
                swap_labels[i] = labels[prev_idx]
                swap_labels[prev_idx] = label

            group_idx += 1
        else:
            print(f'[ERROR] The group {group_id} has more than 2 labels.', file=sys.stderr)
            raise AssertionError('One group must have at most two labels.')

    return swap_labels


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
    main()
