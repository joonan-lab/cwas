#!/usr/bin/env python
"""
Script for an analysis of de novo risk scores using the number of variants for each category
Lasso regression is used for category selection.

For more detailed information, please refer to An et al., 2018 (PMID 30545852).

"""
import argparse
import multiprocessing as mp
import os

import numpy as np
import pandas as pd
import yaml
from glmnet import ElasticNet
from scipy import stats

from utils import cmp_two_arr, get_curr_time, swap_label


def main():
    # Paths to essential configuration files
    curr_dir = os.path.dirname(os.path.abspath(__file__))
    project_dir = os.path.dirname(curr_dir)
    cat_filt_path = os.path.join(project_dir, 'conf', 'category_group.yaml')

    # Parse the arguments
    parser = create_arg_parser()
    args = parser.parse_args()

    # Print the script description
    print(__doc__)

    # Print and check the validity of the settings
    print_args(args)
    check_args_validity(args)
    print()

    # Print the script description
    # Load and parse the input data
    print(f'[{get_curr_time()}, Progress] Parse the categorization result into DataFrame')
    cwas_cat_df = pd.read_table(args.cat_result_path, index_col='SAMPLE')

    # Filter the columns (categories)
    if args.cat_group != 'all':
        print(f'[{get_curr_time()}, Progress] Filter categories and leave \'{args.cat_group}\' categories only')
        with open(cat_filt_path, 'r') as cat_filt_file:
            cat_filt_conf = yaml.safe_load(cat_filt_file)

        cat_filt_dict = cat_filt_conf.get(args.cat_group)
        assert cat_filt_dict is not None, f'{args.cat_group} is not an available category group.'

        for annot_group in cat_filt_dict:
            cat_filt_dict[annot_group] = set(cat_filt_dict[annot_group])

        cwas_cat_df = filter_categories(cwas_cat_df, cat_filt_dict)

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
    case_dnv_cnt, ctrl_dnv_cnt = cnt_case_ctrl_dnv(cwas_cat_vals, sample_types)
    is_rare_cat = ctrl_dnv_cnt < args.rare_cat_cutoff
    rare_cat_vals = cwas_cat_vals[:, is_rare_cat]

    # Determine a training set
    print(f'[{get_curr_time()}, Progress] Divide the samples into training and test set')
    random_seed = 0
    is_train_set = determine_train_set(sample_ids, sample_families, args.train_set_f, random_seed)

    # Train and test a lasso model multiple times
    print(f'[{get_curr_time()}, Progress] Train and test a lasso model to generate de novo risk scores')
    num_parallel = min(mp.cpu_count(), args.num_cv_fold) if args.use_parallel else 1
    coeffs = []
    rsqs = []

    for seed in range(args.num_reg):
        coeff, rsq = \
            lasso_regression(rare_cat_vals, sample_responses, is_train_set, args.num_cv_fold, num_parallel, seed)
        coeffs.append(coeff)
        rsqs.append(rsq)

    m_rsq = np.mean(rsqs)

    if args.do_test:
        # Permutation tests to get null distribution of R squares
        print(f'[{get_curr_time()}, Progress] Permutation tests')
        perm_rsqs = get_perm_rsq(args.num_perm, cwas_cat_vals, sample_types, sample_families, is_train_set,
                                 args.rare_cat_cutoff, args.num_cv_fold, num_parallel)

        # Statistical test (Z test)
        print(f'[{get_curr_time()}, Progress] Statistical test (Z-test)')
        m_perm_rsq = np.mean(perm_rsqs)
        s_perm_rsq = np.std(perm_rsqs, ddof=1)
        z = (m_rsq - m_perm_rsq) / s_perm_rsq
        p = stats.norm.sf(abs(z)) * 2  # Two-sided
    else:
        p = None

    # Make a result DataFrame
    print(f'[{get_curr_time()}, Progress] Make a DataFrame for the de novo risk score analysis')
    num_select = np.sum(np.vectorize(lambda w: w != 0)(coeffs), axis=0)
    m_coeff = np.sum(coeffs, axis=0) / num_select
    is_select_coeff = num_select >= int(args.num_reg * 0.8)
    select_cats = cwas_cat_df.columns.values[is_rare_cat][is_select_coeff]
    cat_case_dnv_cnt = case_dnv_cnt[is_rare_cat][is_select_coeff]
    cat_ctrl_dnv_cnt = ctrl_dnv_cnt[is_rare_cat][is_select_coeff]
    select_m_coeff = m_coeff[is_select_coeff]

    order = select_m_coeff.argsort()[::-1]  # Sort by lasso coefficients in descending order
    result_mat = np.concatenate([
        cat_case_dnv_cnt[:, np.newaxis],
        cat_ctrl_dnv_cnt[:, np.newaxis],
        select_m_coeff[:, np.newaxis]
    ], axis=1)
    result_mat = result_mat[order]
    select_cats = select_cats[order]

    result_df = \
        pd.DataFrame(result_mat, index=select_cats, columns=['Case_DNV_Count', 'Ctrl_DNV_Count', 'Lasso_Coeff'])
    result_df.index.name = 'Category'

    # Write the result
    print(f'[{get_curr_time()}, Progress] Write the result')
    with open(args.outfile_path, 'w') as outfile:
        print(f'#De novo risk score analysis result for {args.cat_group} regions', file=outfile)
        print(f'#Mean R square: {m_rsq * 100:.2f}%', file=outfile)
        if args.do_test:
            print(f'#P-value: {p:.2e}', file=outfile)
        result_df.to_csv(outfile, sep='\t')

    print(f'[{get_curr_time()}, Progress] Done')


def create_arg_parser() -> argparse.ArgumentParser:
    """ Create an argument parser for this script and return it """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-i', '--infile', dest='cat_result_path', required=True, type=str,
                        help='Path of a result of the CWAS categorization')
    parser.add_argument('-s', '--sample_file', dest='sample_file_path', required=True, type=str,
                        help='File listing sample IDs with their families and phenotypes (case or ctrl)')
    parser.add_argument('-a', '--adj_file', dest='adj_file_path', required=False, type=str,
                        help='File that contains adjustment factors for No. DNVs of each sample', default='')
    parser.add_argument('-o', '--outfile', dest='outfile_path', required=False, type=str,
                        help='Path of results of burden tests', default='cwas_denovo_risk_score_result.txt')
    parser.add_argument('--category_group', dest='cat_group', required=False, type=str,
                        choices=['all', 'coding', 'coding-no-ptv', 'noncoding', 'promoter', 'noncoding-wo-promoter'],
                        help='Group of CWAS categories for this analysis (Default: all)', default='all')
    parser.add_argument('--rare_category_cutoff', dest='rare_cat_cutoff', required=False, type=int,
                        help='Rare category cutoff for No. variants of a control (Default: 3)', default=3)
    parser.add_argument('--train_set_fraction', dest='train_set_f', required=False, type=float,
                        help='Fraction of the training set (Default: 0.7)', default=0.7)
    parser.add_argument('--num_regression', dest='num_reg', required=False, type=int,
                        help='No. regression trials to calculate a mean of R squares (Default: 10)', default=10)
    parser.add_argument('--num_cv_fold', dest='num_cv_fold', required=False, type=int,
                        help='No. cross-validation folds (Default: 5)', default=5)
    parser.add_argument('--use_parallel', dest='use_parallel', required=False, type=int, choices={0, 1},
                        help='Use multiprocessing for cross-validation (Default: 1 (True))', default=1)
    parser.add_argument('--do_test', dest='do_test', required=False, type=int, choices={0, 1},
                        help='Do permutation tests to get a p-value for the R square (Default: 0 (False))', default=0)
    parser.add_argument('--num_perm', dest='num_perm', required=False, type=int,
                        help='No. label-swapping permutations for permutation tests (Default: 1,000)', default=1000)
    return parser


def print_args(args: argparse.Namespace):
    print(f'[Setting] Input CWAS categorization result: {args.cat_result_path}')
    print(f'[Setting] List of sample IDs: {args.sample_file_path}')
    print(f'[Setting] List of adjustment factors for No. DNVs of each sample: '
          f'{args.adj_file_path if args.adj_file_path else "None"}')
    print(f'[Setting] Output path: {args.outfile_path}')
    print(f'[Setting] Category group: {args.cat_group}')
    print(f'[Setting] Rare category cutoff for No. variants of a control: {args.rare_cat_cutoff:,d}')
    print(f'[Setting] Fraction of the training set: {args.train_set_f}')
    print(f'[Setting] No. lasso regression trials: {args.num_reg:,d}')
    print(f'[Setting] No. cross-validation folds: {args.num_cv_fold:,d}')
    print(f'[Setting] Use multiprocessing for the cross-validation: {bool(args.use_parallel)}')
    if args.do_test:
        print(f'[Setting] Do the significance test (permutaion test)')
        print(f'[Setting] No. label-swapping permutations: {args.num_perm:,d}')
    else:
        print(f'[Setting] Skip the significance test')


def check_args_validity(args: argparse.Namespace):
    assert os.path.isfile(args.cat_result_path), f'The input file "{args.cat_result_path}" cannot be found.'
    assert os.path.isfile(args.sample_file_path), f'The input file "{args.sample_file_path}" cannot be found.'
    assert args.adj_file_path == '' or os.path.isfile(args.adj_file_path), \
        f'The input file "{args.adj_file_path}" cannot be found.'
    outfile_dir = os.path.dirname(args.outfile_path)
    assert outfile_dir == '' or os.path.isdir(outfile_dir), f'The outfile directory "{outfile_dir}" cannot be found.'
    assert 0 < args.train_set_f < 1, f'The fraction of the training set must be in a range (0, 1).'


def filter_categories(cwas_cat_df: pd.DataFrame, cat_filt_dict: dict) -> pd.DataFrame:
    """ Filter columns (categories) of a DataFrame for the CWAS categorization result
    and leave specific columns according to the input dictionary

    :param cwas_cat_df: A DataFrame that contains the CWAS categorization result
    :param cat_filt_dict: A dictionary which key and value are a group name of annotation terms and
                          the group's annotation terms that must be included in categories, respectively
    :return: A filtered DataFrame
    """
    categories = cwas_cat_df.columns.values
    var_type_annot_set = cat_filt_dict.get('var_type')
    gene_list_annot_set = cat_filt_dict.get('gene_list')
    cons_annot_set = cat_filt_dict.get('cons')
    effect_annot_set = cat_filt_dict.get('effect')
    region_annot_set = cat_filt_dict.get('region')

    def is_passed_category(category):
        var_type, gene_list, cons, effect, region = category.split('_')
        return (var_type_annot_set is None or var_type in var_type_annot_set) and \
               (gene_list_annot_set is None or gene_list in gene_list_annot_set) and \
               (cons_annot_set is None or cons in cons_annot_set) and \
               (effect_annot_set is None or effect in effect_annot_set) and \
               (region_annot_set is None or region in region_annot_set)

    is_passed = np.vectorize(is_passed_category)(categories)
    return cwas_cat_df.loc[:, is_passed]


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


def determine_train_set(samples, sample_groups, frac, random_seed):
    """ Randomly determine a train set from the samples
    """
    uniq_groups = np.unique(sample_groups)
    n_train_group = int(np.rint(len(uniq_groups) * frac))
    np.random.seed(random_seed)
    train_group_set = set(np.random.choice(uniq_groups, size=n_train_group, replace=False))
    is_train_set = np.full(len(samples), False)

    for i in range(len(samples)):
        sample_group = sample_groups[i]
        if sample_group in train_group_set:
            is_train_set[i] = True

    return is_train_set


def lasso_regression(sample_covariates: np.ndarray, sample_responses: np.ndarray, is_train_set: np.ndarray,
                     num_cv_fold: int, num_parallel: int, random_state: int = None) -> (np.ndarray, float):
    """ Lasso regression to generate a de novo risk score """
    # Divide the input data into a train and a test set
    train_covariates = sample_covariates[is_train_set]
    train_responses = sample_responses[is_train_set]
    test_covariates = sample_covariates[~is_train_set]
    test_responses = sample_responses[~is_train_set]

    # Train Lasso model
    lasso_model = ElasticNet(alpha=1, n_lambda=100, standardize=True, n_splits=num_cv_fold, n_jobs=num_parallel,
                             scoring='mean_squared_error', random_state=random_state)
    lasso_model.fit(train_covariates, train_responses)
    opt_model_idx = np.argmax(getattr(lasso_model, 'cv_mean_score_'))
    coeffs = getattr(lasso_model, 'coef_path_')
    opt_coeff = coeffs[:, opt_model_idx]
    opt_lambda = getattr(lasso_model, 'lambda_max_')

    # Prediction using the Lasso model
    pred_responses = lasso_model.predict(test_covariates, lamb=opt_lambda)
    mean_response = np.mean(test_responses)
    rsq = 1 - np.sum((test_responses - pred_responses) ** 2) / np.sum((test_responses - mean_response) ** 2)

    return opt_coeff, rsq


def get_perm_rsq(num_perm: int, sample_cat_vals: np.ndarray, sample_types: np.ndarray, sample_groups: np.ndarray,
                 is_train_set: np.ndarray, rare_cat_cutoff: int, num_cv_fold: int, num_parallel: int) -> list:
    """ Get R squares of the lasso regression after each label swapping trial.
    The length of the returned list equals to the number of the permutations.
    """
    perm_rsqs = []

    for seed in range(num_perm):
        # Label swapping
        swap_sample_types = swap_label(sample_types, sample_groups)
        swap_responses = np.vectorize(lambda sample_type: 1.0 if sample_type == 'case' else -1.0)(swap_sample_types)

        # Filter categories and leave only rare categories (few variants in controls)
        case_dnv_cnt, ctrl_dnv_cnt = cnt_case_ctrl_dnv(sample_cat_vals, swap_sample_types)
        is_rare_cat = ctrl_dnv_cnt < rare_cat_cutoff
        rare_cat_vals = sample_cat_vals[:, is_rare_cat]

        _, rsq = lasso_regression(rare_cat_vals, swap_responses, is_train_set, num_cv_fold, num_parallel, seed)
        perm_rsqs.append(rsq)

    return perm_rsqs


if __name__ == '__main__':
    main()
