#!/usr/bin/env python
"""
# TODO: Add explanation
"""
import argparse
import os
import sys
from datetime import datetime

import numpy as np
import pandas as pd
from glmnet_python.cvglmnet import cvglmnet
from glmnet_python.cvglmnetPredict import cvglmnetPredict


def main():
    # Parse the arguments
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-i', '--infile', dest='cat_result_path', required=True, type=str,
                        help='Path of a result of the CWAS categorization')
    parser.add_argument('-s', '--sample_file', dest='sample_file_path', required=True, type=str,
                        help='File listing sample IDs with their families and phenotypes (case or ctrl)')
    parser.add_argument('-a', '--adj_file', dest='adj_file_path', required=False, type=str,
                        help='File that contains adjustment factors for No. DNVs of each sample', default='')
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
    cwas_cat_indice = cwas_cat_df.index.values  # Sample IDs
    cwas_cat_cols = cwas_cat_df.columns.values
    cwas_cat_vals = cwas_cat_df.values
    sample_info_dict = sample_df.to_dict()

    # Filter categories and leave only rare categories (few variants in controls)
    ctrl_var_cnt_cutoff = 3
    case_dnv_cnt, ctrl_dnv_cnt = cnt_case_ctrl_dnv(cwas_cat_df, sample_df)
    is_rare_cat = ctrl_dnv_cnt < ctrl_var_cnt_cutoff
    cwas_cat_cols = cwas_cat_cols[is_rare_cat]
    cwas_cat_vals = cwas_cat_vals[:, is_rare_cat]

    # Divide the cwas categorization result into training and test set
    # FIXME: This part is for testing reproducibility, so it is not generalized version.
    project_dir = os.path.abspath('..')
    werling_sample_path = os.path.join(project_dir, 'data', 'SuppTable01_Sample_information.xlsx')
    werling_sample_df = pd.read_excel(werling_sample_path, sheet_name='qcMetricsAndCounts_transform.tx')
    werling_sample_ids = werling_sample_df['sampleID'].values
    werling_sample_set = set(werling_sample_ids)

    is_train_set = np.vectorize(lambda sample_id: sample_id in werling_sample_set)(cwas_cat_indice)
    train_sample_ids = cwas_cat_indice[is_train_set]
    train_family_ids = np.unique([sample_info_dict['FAMILY'][sample_id] for sample_id in train_sample_ids])
    train_covariates = cwas_cat_vals[is_train_set]
    train_responses = np.vectorize(
        lambda sample_id: 1.0 if sample_info_dict['PHENOTYPE'][sample_id] == 'case' else -1.0)(train_sample_ids)
    test_sample_ids = cwas_cat_indice[~is_train_set]
    test_covariates = cwas_cat_vals[~is_train_set]
    test_responses = np.vectorize(
        lambda sample_id: 1.0 if sample_info_dict['PHENOTYPE'][sample_id] == 'case' else -1.0)(test_sample_ids)

    # Train and test a lasso model to generate de novo risk scores
    num_trial = 10
    num_fold = 5  # For cross-validation
    num_family_per_fold = div_dist_num(len(train_family_ids), num_fold)
    coeffs = []
    r_sqs = []

    for _ in range(num_trial):
        # Randomly allocate fold IDs to each sample for cross validation
        train_family_ids = np.random.permutation(train_family_ids)
        family_to_fold_id = {}
        start_family_idx = 0

        for i in range(num_fold):
            num_family = num_family_per_fold[i]
            end_family_idx = start_family_idx + num_family

            for j in range(start_family_idx, end_family_idx):
                family_to_fold_id[train_family_ids[j]] = i

            start_family_idx = end_family_idx

        train_fold_ids = np.vectorize(lambda sample_id: family_to_fold_id[sample_info_dict['FAMILY'][sample_id]])(
            train_sample_ids)

        # Train Lasso model
        lasso_model = cvglmnet(x=train_covariates, y=train_responses, ptype='deviance', foldid=train_fold_ids, alpha=1,
                               standardize=True, nlambda=20)
        opt_model_idx = np.argmin(lasso_model['cvm'])
        opt_coeff = lasso_model['glmnet_fit']['beta'][:, opt_model_idx]
        coeffs.append(opt_coeff)

        # Prediction using the Lasso model
        pred_responses = cvglmnetPredict(lasso_model, newx=test_covariates, s='lambda_min').flatten()
        mse = np.mean((pred_responses - test_responses) ** 2)
        r_sq = 1 - mse
        r_sqs.append(r_sq)

    print(np.mean(r_sqs), file=sys.stdout)


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


def cnt_case_ctrl_dnv(cwas_cat_df: pd.DataFrame, sample_df: pd.DataFrame) -> (float, float):
    """ Count the number of the de novo variants for each phenotype, case and control.
    """
    sample_info_dict = sample_df.to_dict()
    sample_ids = cwas_cat_df.index.values
    sample_types = np.asarray([sample_info_dict['PHENOTYPE'][sample_id] for sample_id in sample_ids])
    are_case = sample_types == 'case'

    cat_df_vals = cwas_cat_df.values
    case_dnv_cnt = cat_df_vals[are_case, :].sum(axis=0)
    ctrl_dnv_cnt = cat_df_vals[~are_case, :].sum(axis=0)

    return case_dnv_cnt, ctrl_dnv_cnt


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
