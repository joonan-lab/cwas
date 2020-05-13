#!/usr/bin/env python
"""
# TODO: Add explanation
"""
import argparse
import multiprocessing as mp
import os
import sys
from datetime import datetime
from functools import partial

import numpy as np
import pandas as pd
from glmnet_python.cvglmnet import cvglmnet
from glmnet_python.cvglmnetPredict import cvglmnetPredict
from scipy import stats


def main():
    # Parse the arguments
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-i', '--infile', dest='cat_result_path', required=True, type=str,
                        help='Path of a result of the CWAS categorization')
    parser.add_argument('-s', '--sample_file', dest='sample_file_path', required=True, type=str,
                        help='File listing sample IDs with their families and phenotypes (case or ctrl)')
    parser.add_argument('-a', '--adj_file', dest='adj_file_path', required=False, type=str,
                        help='File that contains adjustment factors for No. DNVs of each sample', default='')
    parser.add_argument('-p', '--num_proc', dest='num_proc', required=False, type=int,
                        help='Number of processes used in permutation tests',
                        default=1)
    args = parser.parse_args()

    # Print the script description
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

    # Arrays and dictionary from the DataFrame (in order to improving performance)
    cwas_cat_vals = cwas_cat_df.values
    sample_ids = cwas_cat_df.index.values
    sample_info_dict = sample_df.to_dict()
    sample_types = np.vectorize(lambda sample_id: sample_info_dict['PHENOTYPE'][sample_id])(sample_ids)
    sample_responses = np.vectorize(lambda sample_type: 1.0 if sample_type == 'case' else -1.0)(sample_types)
    sample_families = np.vectorize(lambda sample_id: sample_info_dict['FAMILY'][sample_id])(sample_ids)

    # Filter categories and leave only rare categories (few variants in controls)
    print(f'[{get_curr_time()}, Progress] Filter categories and leave only rare categories')
    var_cnt_cutoff = 3
    case_dnv_cnt, ctrl_dnv_cnt = cnt_case_ctrl_dnv(cwas_cat_vals, sample_types)
    is_rare_cat = ctrl_dnv_cnt < var_cnt_cutoff
    rare_cat_vals = cwas_cat_vals[:, is_rare_cat]

    # Determine a training set
    print(f'[{get_curr_time()}, Progress] Divide the samples into training and test set')
    # FIXME: This part is for testing reproducibility, so it is not generalized version.
    project_dir = os.path.abspath('..')
    werling_sample_path = os.path.join(project_dir, 'data', 'SuppTable01_Sample_information.xlsx')
    werling_sample_df = pd.read_excel(werling_sample_path, sheet_name='qcMetricsAndCounts_transform.tx')
    werling_sample_ids = werling_sample_df['sampleID'].values
    werling_sample_set = set(werling_sample_ids)
    is_train_set = np.vectorize(lambda sample_id: sample_id in werling_sample_set)(sample_ids)

    # Train and test a lasso model multiple times
    print(f'[{get_curr_time()}, Progress] Train and test a lasso model to generate de novo risk scores')
    num_trial = 10
    num_fold = 5  # For cross-validation
    num_parallel = min(mp.cpu_count(), num_fold)
    coeffs = []
    rsqs = []

    for _ in range(num_trial):
        coeff, rsq = \
            lasso_regression(rare_cat_vals, sample_responses, sample_families, is_train_set, num_fold, num_parallel)
        coeffs.append(coeff)
        rsqs.append(rsq)

    # Permutation tests to get null distribution of R squares
    print(f'[{get_curr_time()}, Progress] Permutation tests')
    num_perm = 1000

    if args.num_proc == 1:
        perm_rsqs = get_perm_rsq(num_perm, cwas_cat_vals, sample_types, sample_families, is_train_set,
                                 var_cnt_cutoff, num_fold, num_parallel)
    else:
        num_perms = div_dist_num(num_perm, args.num_proc)
        pool = mp.Pool(args.num_proc)
        proc_outputs = \
            pool.map(
                partial(get_perm_rsq,
                        sample_cat_vals=cwas_cat_vals,
                        sample_types=sample_types,
                        sample_groups=sample_families,
                        is_train_set=is_train_set,
                        var_cnt_cutoff=var_cnt_cutoff,
                        num_fold=num_fold,
                        num_parallel=1
                        ),
                num_perms
            )
        pool.close()
        pool.join()

        perm_rsqs = []

        for proc_output in proc_outputs:
            perm_rsqs += proc_output

    # Statistical test (Z test)
    m_rsq = np.mean(rsqs)
    m_perm_rsq = np.mean(perm_rsqs)
    s_perm_rsq = np.std(perm_rsqs)
    z = (m_rsq - m_perm_rsq) / s_perm_rsq
    p = stats.norm.sf(abs(z)) * 2  # Two-sided

    print(f'[{get_curr_time()}, Progress] Print the results')
    print(f'Mean R square\t{m_rsq * 100}%')
    print(f'P-value\t{p:.2e}')


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


def cnt_case_ctrl_dnv(sample_cat_vals: np.ndarray, sample_types: np.ndarray) -> (float, float):
    """ Count the number of the de novo variants for each phenotype, case and control.
    """
    are_case = sample_types == 'case'
    case_dnv_cnt = sample_cat_vals[are_case, :].sum(axis=0)
    ctrl_dnv_cnt = sample_cat_vals[~are_case, :].sum(axis=0)

    return case_dnv_cnt, ctrl_dnv_cnt


def lasso_regression(sample_covariates: np.ndarray, sample_responses: np.ndarray, sample_groups: np.ndarray,
                     is_train_set: np.ndarray, num_fold: int, num_parallel: int) -> (np.ndarray, float):
    """ Lasso regression to generate a de novo risk score """
    # Divide the input data into a train and a test set
    train_covariates = sample_covariates[is_train_set]
    train_responses = sample_responses[is_train_set]
    train_groups = sample_groups[is_train_set]
    test_covariates = sample_covariates[~is_train_set]
    test_responses = sample_responses[~is_train_set]

    # Allocate fold IDs for cross validation
    group_to_fold_id = allocate_fold_ids(train_groups, num_fold)
    train_fold_ids = np.vectorize(lambda family_id: group_to_fold_id[family_id])(train_groups)

    # Train Lasso model
    lasso_model = cvglmnet(x=train_covariates, y=train_responses, ptype='deviance', foldid=train_fold_ids, alpha=1,
                           standardize=True, nlambda=20, parallel=num_parallel)
    opt_model_idx = np.argmin(lasso_model['cvm'])
    coeff = lasso_model['glmnet_fit']['beta'][:, opt_model_idx]

    # Prediction using the Lasso model
    pred_responses = cvglmnetPredict(lasso_model, newx=test_covariates, s='lambda_min').flatten()
    mse = np.mean((pred_responses - test_responses) ** 2)
    rsq = 1 - mse

    return coeff, rsq


def get_perm_rsq(num_perm: int, sample_cat_vals: np.ndarray, sample_types: np.ndarray, sample_groups: np.ndarray,
                 is_train_set: np.ndarray, var_cnt_cutoff: int, num_fold: int, num_parallel: int) -> list:
    """ Get R squares of the lasso regression after each label swapping trial.
    The length of the returned list equals to the number of the permutations.
    """
    perm_rsqs = []

    for _ in range(num_perm):
        # Label swapping
        swap_sample_types = swap_label(sample_types, sample_groups)
        swap_responses = np.vectorize(lambda sample_type: 1.0 if sample_type == 'case' else -1.0)(swap_sample_types)

        # Filter categories and leave only rare categories (few variants in controls)
        case_dnv_cnt, ctrl_dnv_cnt = cnt_case_ctrl_dnv(sample_cat_vals, swap_sample_types)
        is_rare_cat = ctrl_dnv_cnt < var_cnt_cutoff
        rare_cat_vals = sample_cat_vals[:, is_rare_cat]

        _, rsq = lasso_regression(rare_cat_vals, swap_responses, sample_groups, is_train_set, num_fold, num_parallel)
        perm_rsqs.append(rsq)

    return perm_rsqs


def allocate_fold_ids(samples: np.ndarray, num_fold: int) -> dict:
    """ Randomly allocate fold IDs to each sample and
    return a dictionary which key and value are sample and its ID, respectively.
    """
    unique_samples = np.unique(samples)
    sample_size = len(unique_samples)
    sample_to_fold_id = {}
    sizes_per_fold = div_dist_num(sample_size, num_fold)
    perm_samples = np.random.permutation(unique_samples)
    start_idx = 0

    for fold_id in range(num_fold):
        end_idx = start_idx + sizes_per_fold[fold_id]

        for j in range(start_idx, end_idx):
            sample_to_fold_id[perm_samples[j]] = fold_id

        start_idx = end_idx

    return sample_to_fold_id


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


def cmp_two_arr(array1: np.ndarray, array2: np.ndarray) -> bool:
    """ Return True if two arrays have the same items regardless of the order, else return False """
    if len(array1) != len(array2):
        return False

    array1_item_set = set(array1)

    for item in array2:
        if item not in array1_item_set:
            return False

    return True


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


if __name__ == '__main__':
    main()
