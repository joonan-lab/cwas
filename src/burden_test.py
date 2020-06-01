#!/usr/bin/env python
"""
Script to run burden tests for de novo variants (DNVs) of cases and controls in each CWAS category.
Binomial tests and permutation tests are available in this script.

For more detailed information, please refer to An et al., 2018 (PMID 30545852).

"""
import argparse
import multiprocessing as mp
import os
from functools import partial

import numpy as np
import pandas as pd
from scipy.stats import binom_test

from utils import cmp_two_arr, div_dist_num, get_curr_time, swap_label


def main():
    # Parse the arguments
    parser = create_arg_parser()
    args = parser.parse_args()

    # Print the script description
    print(__doc__)

    # Print and check the validity of the settings
    print_args(args)
    check_args_validity(args)
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


def create_arg_parser() -> argparse.ArgumentParser:
    """ Create an argument parser for this script and return it """
    # Create the top-level argument parser
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(description='Types of burden tests', dest='test_type')

    def add_common_args(subparser: argparse.ArgumentParser):
        """ Add common arguments to the subparser """
        subparser.add_argument('-i', '--infile', dest='cat_result_path', required=True, type=str,
                               help='Path of a result of the CWAS categorization')
        subparser.add_argument('-s', '--sample_file', dest='sample_file_path', required=True, type=str,
                               help='File listing sample IDs with their families and sample_types (case or ctrl)')
        subparser.add_argument('-a', '--adj_file', dest='adj_file_path', required=False, type=str,
                               help='File that contains adjustment factors for No. DNVs of each sample', default='')

    # Create the parser for binomial tests
    parser_binom = subparsers.add_parser('binom', description='Burden tests via binomial tests',
                                         help='Binomial tests (arg "binom -h" for usage)')
    add_common_args(parser_binom)
    parser_binom.add_argument('-o', '--outfile', dest='outfile_path', required=False, type=str,
                              help='Path of results of burden tests', default='cwas_burden_binom_result.txt')

    # Create the parser for permutation tests
    parser_perm = subparsers.add_parser('perm', description='Burden tests via permutation tests',
                                        help='Permutation tests (arg "perm -h" for usage)')
    add_common_args(parser_perm)
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

    return parser


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


def print_args(args: argparse.Namespace):
    print(f'[Setting] Types of burden tests: {"Binomial test" if args.test_type == "binom" else "Permutation test"}')
    print(f'[Setting] Input CWAS categorization result: {args.cat_result_path}')
    print(f'[Setting] List of sample IDs: {args.sample_file_path}')
    print(f'[Setting] List of adjustment factors for No. DNVs of each sample: '
          f'{args.adj_file_path if args.adj_file_path else "None"}')
    print(f'[Setting] Output path: {args.outfile_path}')

    if args.test_type == 'perm':
        print(f'[Setting] No. label-swapping permutations: {args.num_perm:,d}')
        print(f'[Setting] No. processes for the tests: {args.num_proc:,d}')
        print(f'[Setting] Permutation RR output path: {args.perm_rr_path if args.perm_rr_path else "None"}')


def check_args_validity(args: argparse.Namespace):
    assert os.path.isfile(args.cat_result_path), f'The input file "{args.cat_result_path}" cannot be found.'
    assert os.path.isfile(args.sample_file_path), f'The input file "{args.sample_file_path}" cannot be found.'
    assert args.adj_file_path == '' or os.path.isfile(args.adj_file_path), \
        f'The input file "{args.adj_file_path}" cannot be found.'
    outfile_dir = os.path.dirname(args.outfile_path)
    assert outfile_dir == '' or os.path.isdir(outfile_dir), f'The outfile directory "{outfile_dir}" cannot be found.'

    if args.test_type == 'perm':
        assert 1 <= args.num_proc <= mp.cpu_count(), \
            f'Invalid number of processes "{args.num_proc:,d}". It must be in the range [1, {mp.cpu_count()}].'
        assert args.num_perm >= args.num_proc, f'No. processes must be equal or less than No. permutations.'
        perm_rr_dir = os.path.dirname(args.perm_rr_path)
        assert perm_rr_dir == '' or os.path.isdir(perm_rr_dir), \
            f'The outfile directory for RRs from permutations "{perm_rr_dir}" cannot be found.'


def run_burden_binom(cwas_cat_df: pd.DataFrame, sample_df: pd.DataFrame) -> pd.DataFrame:
    """ Function for burden tests via binomial tests

    :param cwas_cat_df: A DataFrame that contains the result of CWAS categorization
    :param sample_df: A DataFrame listing sample IDs with their families and sample_types
    :return: A DataFrame that contains binomial p-values and other statistics for each CWAS category
    """
    # Count the number of de novo variants (DNV) for cases and controls
    cwas_cat_vals = cwas_cat_df.values
    sample_info_dict = sample_df.to_dict()
    sample_ids = cwas_cat_df.index.values
    sample_types = np.vectorize(lambda sample_id: sample_info_dict['PHENOTYPE'][sample_id])(sample_ids)
    case_dnv_cnt, ctrl_dnv_cnt = cnt_case_ctrl_dnv(cwas_cat_vals, sample_types)
    dnv_cnt_arr = np.concatenate([case_dnv_cnt[:, np.newaxis], ctrl_dnv_cnt[:, np.newaxis]], axis=1)

    # Make a DataFrame for the results of binomial tests
    burden_df = \
        pd.DataFrame(dnv_cnt_arr, index=cwas_cat_df.columns.values, columns=['Case_DNV_Count', 'Ctrl_DNV_Count'])
    burden_df.index.name = 'Category'
    burden_df['Relative_Risk'] = case_dnv_cnt / ctrl_dnv_cnt

    # Binomial tests
    def binom_two_tail(n1, n2):
        return binom_test(x=n1, n=n1 + n2, p=0.5, alternative='two-sided')

    def binom_one_tail(n1, n2):
        return binom_test(x=n1, n=n1 + n2, p=0.5, alternative='greater') if n1 > n2 \
            else binom_test(x=n2, n=n1 + n2, p=0.5, alternative='greater')

    burden_df['P'] = \
        np.vectorize(binom_two_tail)(case_dnv_cnt.round(), ctrl_dnv_cnt.round())
    burden_df['P_1side'] = \
        np.vectorize(binom_one_tail)(case_dnv_cnt.round(), ctrl_dnv_cnt.round())

    return burden_df


def run_burden_perm(cwas_cat_df: pd.DataFrame, sample_df: pd.DataFrame, num_perm: int, num_proc: int) \
        -> (pd.DataFrame, pd.DataFrame):
    """ Function for burden tests via permutation tests

    :param cwas_cat_df: A DataFrame that contains the result of CWAS categorization
    :param sample_df: A DataFrame listing sample IDs with their families and sample_types
    :param num_perm: The number of label-swapping permutation trials
    :param num_proc: The number of processes used in this function (for multiprocessing)
    :returns:
        1. A DataFrame that contains permutation p-values and other statistics for each CWAS category
        2. A DataFrame that contains relative risks for each category from each permutation trial
    """
    # Arrays and a dictionary from the input DataFrames
    cwas_cat_vals = cwas_cat_df.values
    sample_info_dict = sample_df.to_dict()
    sample_ids = cwas_cat_df.index.values
    sample_types = np.vectorize(lambda sample_id: sample_info_dict['PHENOTYPE'][sample_id])(sample_ids)
    family_ids = np.vectorize(lambda sample_id: sample_info_dict['FAMILY'][sample_id])(sample_ids)

    # Calculate relative risks from label-swapping permutations
    if num_proc == 1:
        perm_rr_list = cal_perm_rr(num_perm, cwas_cat_vals, sample_types, family_ids)
    else:
        num_perms = div_dist_num(num_perm, num_proc)
        pool = mp.Pool(num_proc)
        proc_outputs = \
            pool.map(
                partial(cal_perm_rr,
                        sample_cat_vals=cwas_cat_vals,
                        sample_types=sample_types,
                        family_ids=family_ids
                        ),
                num_perms
            )
        pool.close()
        pool.join()

        perm_rr_list = []

        for proc_output in proc_outputs:
            perm_rr_list += proc_output

    perm_rrs = np.concatenate(perm_rr_list, axis=0)
    case_dnv_cnt, ctrl_dnv_cnt = cnt_case_ctrl_dnv(cwas_cat_vals, sample_types)
    rr = case_dnv_cnt / ctrl_dnv_cnt  # Original relative risks

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


def cnt_case_ctrl_dnv(sample_cat_vals: np.ndarray, sample_types: np.ndarray) -> (float, float):
    """ Count the number of the de novo variants for each phenotype, case and control.
    """
    are_case = sample_types == 'case'
    case_dnv_cnt = sample_cat_vals[are_case, :].sum(axis=0)
    ctrl_dnv_cnt = sample_cat_vals[~are_case, :].sum(axis=0)

    return case_dnv_cnt, ctrl_dnv_cnt


def cal_perm_rr(num_perm: int, sample_cat_vals: np.ndarray, sample_types: np.ndarray, family_ids: np.ndarray) -> list:
    """ Calculate relative risks of each category in each permutation trial.
    The length of the returned list equals to the number of the permutations.
    """
    perm_rr_list = []

    for _ in range(num_perm):
        swap_sample_types = swap_label(sample_types, family_ids)
        case_dnv_cnt, ctrl_dnv_cnt = cnt_case_ctrl_dnv(sample_cat_vals, swap_sample_types)
        perm_rr = case_dnv_cnt / ctrl_dnv_cnt
        perm_rr_list.append(perm_rr[np.newaxis, :])

    return perm_rr_list


if __name__ == "__main__":
    main()
