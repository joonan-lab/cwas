import argparse
import os, sys
import pandas as pd
import numpy as np
from pathlib import Path
import rpy2.robjects as ro
from rpy2.robjects import numpy2ri
from rpy2.robjects.packages import importr

from sklearn.preprocessing import StandardScaler
from sklearn.metrics import r2_score

import cwas.utils.log as log
from cwas.core.common import cmp_two_arr
from cwas.utils.check import check_is_file, check_num_proc, check_is_dir
from cwas.runnable import Runnable
from typing import Optional, Tuple
from contextlib import contextmanager
from collections import defaultdict
import matplotlib.pyplot as plt
import polars as pl
import re
import parmap
from tqdm import tqdm
from functools import partial
import zarr
from concurrent.futures import ProcessPoolExecutor
import gc


class RiskScore(Runnable):
    def __init__(self, args: argparse.Namespace):
        super().__init__(args)
        self._sample_info = None
        self._categorization_root = None
        self._categorization_result = None
        self._sample_ids = None
        self._categories = None
        self._adj_factor = None
        self._category_set_path = None
        self._category_set = None
        self._datasets = None
        self._covariates = None
        self._test_covariates = None
        self._response = None
        self._test_response = None
        self._result_dict = defaultdict(dict)
        self._permutation_dict = defaultdict(dict)
        self._filtered_combs = None
        self.cv_glmnet = importr("glmnet").cv_glmnet
        self._annotation_list = None
        self._result_for_loop = defaultdict(dict)
        self._result_for_n_of_one_leave = defaultdict(dict)
    
    @staticmethod
    def _print_args(args: argparse.Namespace):
        log.print_arg(
            "Categorization result file", 
            args.categorization_result_path
                if args.categorization_result_path
                else "Not specified: $CATEGORIZATION_RESULT will be used")
        log.print_arg("Sample information file", args.sample_info_path)
        log.print_arg("Adjustment factor list", args.adj_factor_path)
        log.print_arg(
            "Category set file", 
            args.category_set_path)
        log.print_arg(
            "Domain list", 
            args.domain_list)        
        if args.tag:
            log.print_arg("Output tag (prefix of output files)", args.tag)
        log.print_arg("If the number of carriers is used for calculating risk score or not", args.use_n_carrier)
        if args.do_loop:
            log.print_arg("Use each annotation from functional annotation to calculate risk score", args.do_loop)
        if args.n_of_one_leave:
            log.print_arg("Exclude one annotation from functional annotation and functional score and calculate risk score", args.n_of_one_leave)
        log.print_arg(
            "Threshold for selecting rare categories",
            f"{args.ctrl_thres: ,d}",
        )
        log.print_arg("Fraction of the training set", f"{args.train_set_f: ,f}")
        log.print_arg(
            "No. regression trials to calculate a mean of R squared values",
            f"{args.num_reg: ,d}",
        )
        log.print_arg(
            "No. folds for Cross-Vadidation",
            f"{args.fold: ,d}",
        )
        #log.print_arg("Use Logistic regression", args.logistic)
        log.print_arg(
            "No. permutation used to calculate the p-value",
            f"{args.n_permute: ,d}",
        )
        log.print_arg("Skip the permutation test", args.predict_only)
        log.print_arg(
            "No. worker processes",
            f"{args.num_proc: ,d}",
        )

    @staticmethod
    def _check_args_validity(args: argparse.Namespace):
        check_is_file(args.sample_info_path)
        check_num_proc(args.num_proc)
        if args.categorization_result_path:
            check_is_dir(args.categorization_result_path)
        if args.adj_factor_path is not None:
            check_is_file(args.adj_factor_path)
        if args.category_set_path:
            check_is_file(args.category_set_path)
    
    @property
    def categorization_result_path(self) -> Path:
        return (
            self.args.categorization_result_path.resolve()
            if self.args.categorization_result_path 
            else Path(self.get_env("CATEGORIZATION_RESULT"))
        )
    
    @property
    def category_set_path(self) -> Path:
        return self.args.category_set_path.resolve()

    @property
    def sample_info_path(self) -> Path:
        return self.args.sample_info_path.resolve()

    @property
    def out_dir(self) -> Path:
        return(self.args.output_dir_path.resolve())

    @property
    def adj_factor_path(self) -> Optional[Path]:
        return (
            self.args.adj_factor_path.resolve()
            if self.args.adj_factor_path
            else None
        )

    @property
    def domain_list(self) -> str:
        if self.args.domain_list == 'all':
            return ['all']
        elif self.args.domain_list=='run_all':
            all_domains = ['all'] + [col[3:] for col in self.category_set.columns if col.startswith('is_')]
            return all_domains
        else:
            if 'all' in self.args.domain_list:
                all_domains = [col[3:] for col in self.category_set.columns if col.startswith('is_')]
                matching_values = ['all']+[self._check_domain_list(str.lower(d.strip()), all_domains) for d in self.args.domain_list.split(',')]
                return matching_values
            else:
                all_domains = [col[3:] for col in self.category_set.columns if col.startswith('is_')]
                matching_values = [self._check_domain_list(str.lower(d.strip()), all_domains) for d in self.args.domain_list.split(',')]
                return matching_values

    def _check_domain_list(self, d, all_domain_list):
        if not d in map(str.lower, all_domain_list):
            raise ValueError(
                "Invalid domain name: "
                "{}".format(d)
            )
        else:
            idx = list(map(str.lower, all_domain_list)).index(d)
            return all_domain_list[idx]

    @property
    def tag(self) -> str:
        return self.args.tag
    
    @property
    def do_loop(self) -> bool:
        return self.args.do_loop

    @property
    def n_of_one_leave(self) -> bool:
        return False if self.do_loop else self.args.n_of_one_leave

    @property
    def use_n_carrier(self) -> bool:
        return self.args.use_n_carrier

    @property
    def ctrl_thres(self) -> int:
        return self.args.ctrl_thres

    @property
    def train_set_f(self) -> float:
        return self.args.train_set_f

    @property
    def num_reg(self) -> int:
        return self.args.num_reg

    @property
    def fold(self) -> int:
        return self.args.fold

    @property
    def n_permute(self) -> int:
        return self.args.n_permute
    
    @property
    def predict_only(self) -> bool:
        return self.args.predict_only

    @property
    def num_proc(self) -> int:
        return self.args.num_proc

    @property
    def seed(self) -> int:
        return self.args.seed
    
    @property
    def categorization_root(self):
        if self._categorization_root is None:
            self._categorization_root = zarr.open(self.categorization_result_path, mode='r')
        return self._categorization_root
    
    @property
    def categorization_result(self):
        if self._categorization_result is None:
            self._categorization_result = self.categorization_root['data']
        return self._categorization_result

    @property
    def sample_ids(self):
        if self._sample_ids is None:
            self._sample_ids = self.categorization_root['metadata'].attrs['sample_id']
        return self._sample_ids

    @property
    def categories(self):
        if self._categories is None:
            self._categories = self.categorization_root['metadata'].attrs['category']
        return self._categories

    @property
    def annotation_list(self):
        if self._annotation_list is None:
            self._annotation_list = sorted([value for value in np.unique(self.category_set['functional_annotation']) if value != 'Any'])
        return self._annotation_list
            
    @property
    def category_set(self) -> pd.DataFrame:
        if self._category_set is None:
            self._category_set = pd.read_csv(self.category_set_path, sep='\t')
            if sorted(self._category_set['Category'].tolist()) != sorted(self.categories):
                raise ValueError("The categories in the 'input_file' and 'category_set' do not match.")
            self._category_set = self._category_set.loc[self._category_set["Category"].isin(self.categories)]
            self._category_set['Category'] = pd.Categorical(self._category_set['Category'],
                                                            categories=self.categories,
                                                            ordered=True)
            self._category_set.sort_values('Category', ignore_index=True, inplace=True)
        return self._category_set
    
    @property
    def sample_info(self) -> pd.DataFrame:
        if self._sample_info is None:
            self._sample_info = pd.read_table(self.sample_info_path,
                                              index_col="SAMPLE",
                                              dtype={"SAMPLE": str},
                                              sep="\t").reindex(self.sample_ids)
            if ("SET" not in self._sample_info.columns):
                log.print_log("LOG", 
                              "No 'SET' column in sample information file. "
                              "Training and test sets will be assigned randomly.")
                
                case_count = sum(self._sample_info['PHENOTYPE'] == 'case')
                ctrl_count = sum(self._sample_info['PHENOTYPE'] == 'ctrl')
                
                self._case_f = round(case_count*self.train_set_f)
                self._ctrl_f = round(ctrl_count*self.train_set_f)
                
                case_train_idx = self._sample_info.loc[self._sample_info.PHENOTYPE=='case'].sample(n=self._case_f, random_state=self.seed).index
                ctrl_train_idx = self._sample_info.loc[self._sample_info.PHENOTYPE=='ctrl'].sample(n=self._ctrl_f, random_state=self.seed).index
                self._sample_info['SET'] = ''
                self._sample_info.loc[case_train_idx, 'SET'] = 'training'
                self._sample_info.loc[ctrl_train_idx, 'SET'] = 'training'
                self._sample_info.loc[self._sample_info['SET'] == '', 'SET'] = 'test'
            else:
                self._case_f = sum((self._sample_info['PHENOTYPE'] == 'case') & (self._sample_info['SET'] == 'training'))
                self._ctrl_f = sum((self._sample_info['PHENOTYPE'] == 'ctrl') & (self._sample_info['SET'] == 'training'))
              
            log.print_log("LOG",
                        "Use {} cases and {} controls for training set".format(self._case_f, self._ctrl_f)
                        )
        return self._sample_info

    @property
    def adj_factor(self) -> pd.DataFrame:
        if (self._adj_factor is None) and (self.adj_factor_path is not None):
            self._adj_factor = pd.read_table(self.adj_factor_path,
                                             index_col="SAMPLE",
                                             dtype={"SAMPLE": str},
                                             sep="\t")
            if not cmp_two_arr(self.sample_ids, self._adj_factor.index.values):
                raise ValueError(
                    "The sample IDs from the adjustment factor list are"
                    "not the same with the sample IDs "
                    "from the categorization result."
                )
            self._adj_factor = self._adj_factor.reindex(self.sample_ids) 
        return self._adj_factor
    
    @property
    def coef_path(self) -> Path:
        tag = '' if self.tag is None else ''.join([self.tag, '_'])
        f_name = re.sub(r'.categorization_result\.zarr\.gz|.categorization_result\.zarr', f'.lasso_coef_{tag}thres_{self.ctrl_thres}.txt', str(self.categorization_result_path.name))
        return Path(
            f"{self.out_dir}/" +
            f"{f_name}"
        )
    
    @property
    def result_path(self) -> Path:
        tag = '' if self.tag is None else ''.join([self.tag, '_'])
        f_name = re.sub(r'.categorization_result\.zarr\.gz|.categorization_result\.zarr', f'.lasso_results_{tag}thres_{self.ctrl_thres}.txt', str(self.categorization_result_path.name))
        return Path(
            f"{self.out_dir}/" +
            f"{f_name}"
        )

    @property
    def null_model_path(self) -> Path:
        tag = '' if self.tag is None else ''.join([self.tag, '_'])
        f_name = re.sub(r'.categorization_result\.zarr\.gz|.categorization_result\.zarr', f'.lasso_null_models_{tag}thres_{self.ctrl_thres}.txt', str(self.categorization_result_path.name))
        return Path(
            f"{self.out_dir}/" +
            f"{f_name}"
        )
        
    @property
    def plot_path(self) -> Optional[Path]:
        tag = '' if self.tag is None else ''.join([self.tag, '_'])
        f_name = re.sub(r'.categorization_result\.zarr\.gz|.categorization_result\.zarr', f'.lasso_histogram_{tag}thres_{self.ctrl_thres}.pdf', str(self.categorization_result_path.name))
        return Path(
            f"{self.out_dir}/" +
            f"{f_name}"
        )

    def run(self):
        self.prepare()
        if self.do_loop:
            for i in self.annotation_list:
                log.print_progress("Start loop for each annotation")
                log.print_progress(f"Generate risk scores for annotation: {i}")
                self.filtered_category_set = self.category_set[self.category_set['functional_annotation'] == i]
                self.risk_scores()
                if not self.predict_only:
                    self.permute_pvalues()
                self.gather_results_for_loop(i)
            self.save_results_for_loop()
        if self.n_of_one_leave:
            for i in self.annotation_list:
                log.print_progress("Start N of one leave for each annotation")
                log.print_progress(f"Generate risk scores excluding annotation: {i}")
                self.filtered_category_set = self.category_set[self.category_set['functional_annotation'] != i]
                self.risk_scores()
                if not self.predict_only:
                    self.permute_pvalues()
                self.gather_results_for_loop(i)
                self.save_results_for_loop()
        if not (self.do_loop or self.n_of_one_leave):
            self.filtered_category_set = self.category_set
            self.risk_scores()
            if not self.predict_only:
                self.permute_pvalues()
            self.save_results()
        self.update_env()
        log.print_progress("Done.")

    def prepare(self):
        if not cmp_two_arr(self.sample_ids, self.sample_info.index.values):
            raise ValueError(
                "The sample IDs from the sample information are "
                "not the same with the sample IDs "
                "from the categorization result."
            )

    def risk_scores(self):
        """Generate risk scores for various seeds """
        log.print_progress(self.risk_scores.__doc__)

        #domain_values = [domain.strip() for domain in self.domain_list.split(',')]
        seeds = np.arange(self.seed, self.seed + self.num_reg * 10, 10)
        num_proc2 = len(seeds) if self.num_proc > len(seeds) else self.num_proc
        pool = ProcessPoolExecutor(max_workers=num_proc2)
        
        for domain in self.domain_list:
            log.print_progress(f"Generate risk score for the domain: {domain}")
            filtered_combs = self.filtered_category_set.loc[self.filtered_category_set['is_'+domain]==1]["Category"] if domain != 'all' else pd.Series(self.categories)
            
            rare_categories, cov, test_cov, response, test_response = self._create_covariates(seed=self.seed,
                                                                                              swap_label=False,
                                                                                              filtered_combs=filtered_combs)
            _risk_score_per_category_ = partial(self._risk_score_per_category,
                                                swap_label = False,
                                                rare_categories = rare_categories,
                                                cov = cov,
                                                test_cov = test_cov,
                                                response = response,
                                                test_response = test_response,
                                                filtered_combs = filtered_combs)
            map_result = pool.map(_risk_score_per_category_, seeds)
            self._result_dict[domain] = {key: value for x in map_result for key, value in x.items()}
            gc.collect()
            
    def permute_pvalues(self):
        """Run LassoCV to get permutated pvalues"""
        log.print_progress(self.permute_pvalues.__doc__)
                    
        seeds = np.arange(self.seed, self.seed+self.n_permute)
        pool = ProcessPoolExecutor(max_workers=self.num_proc)

        for domain in self.domain_list:
            log.print_progress(f"Generate permutation p-values for the domain: {domain}")
            filtered_combs = self.filtered_category_set.loc[self.filtered_category_set['is_'+domain]==1]['Category'] if domain != 'all' else pd.Series(self.categories)
            
            _risk_score_per_category_ = partial(self._risk_score_per_category,
                                                swap_label = True,
                                                rare_categories = None,
                                                cov = None,
                                                test_cov = None,
                                                response = None,
                                                test_response = None,
                                                filtered_combs = filtered_combs)
            #map_result = parmap.map(_risk_score_per_category_, seeds, pm_pbar=True, pm_processes=self.num_proc)
            map_result = list(tqdm(pool.map(_risk_score_per_category_, seeds), total=len(seeds), desc="Permutation p-values"))
            self._permutation_dict[domain] = {key: value for x in map_result for key, value in x.items()}
            gc.collect()
    
    def _create_covariates(self, seed: int, swap_label: bool = False, filtered_combs = None):
        """Create covariates for risk score model"""
        if swap_label:
            np.random.seed(seed)
            
            perm_sample_info = self.sample_info.copy()
            perm_sample_info['Perm_PHENOTYPE'] = np.random.permutation(self.sample_info['PHENOTYPE'])
            perm_sample_info.drop(columns=['PHENOTYPE', 'SET'], inplace=True)
                        
            perm_case_train_idx = perm_sample_info.loc[perm_sample_info.Perm_PHENOTYPE=='case'].sample(n=self._case_f, random_state=seed).index
            perm_ctrl_train_idx = perm_sample_info.loc[perm_sample_info.Perm_PHENOTYPE=='ctrl'].sample(n=self._ctrl_f, random_state=seed).index
                        
            perm_sample_info['Perm_SET'] = ''
            perm_sample_info.loc[perm_case_train_idx, 'Perm_SET'] = 'training'
            perm_sample_info.loc[perm_ctrl_train_idx, 'Perm_SET'] = 'training'
            perm_sample_info.loc[perm_sample_info['Perm_SET'] == '', 'Perm_SET'] = 'test'
                
            ctrls_coords = list(np.where(perm_sample_info["Perm_PHENOTYPE"]=='ctrl')[0])
            
            train_samples_coords = list(np.where(perm_sample_info['Perm_SET']=='training')[0])
            test_samples_coords = list(np.where(perm_sample_info['Perm_SET']=='test')[0])
            
            response = (perm_sample_info.iloc[train_samples_coords]['Perm_PHENOTYPE'] == 'case').values
            test_response = (perm_sample_info.iloc[test_samples_coords]['Perm_PHENOTYPE'] == 'case').values
        else:
            ctrls_coords = list(np.where(self.sample_info.PHENOTYPE=='ctrl')[0]) # 4348
        
            train_samples_coords = list(np.where(self.sample_info['SET']=='training')[0])
            test_samples_coords = list(np.where(self.sample_info['SET']=='test')[0])
            
            response = (self.sample_info.iloc[train_samples_coords]['PHENOTYPE'] == 'case').values
            test_response = (self.sample_info.iloc[test_samples_coords]['PHENOTYPE'] == 'case').values
            
        filtered_combs_coords = filtered_combs.index.tolist()
        
        if not self.adj_factor is None:
            if self.use_n_carrier:
                ctrl_categorization_result = (self.categorization_result.get_orthogonal_selection((ctrls_coords, filtered_combs_coords)) > 0).astype('uint8')
                ctrl_adj_factors = self.adj_factor.iloc[ctrls_coords]['AdjustFactor'].values[:,np.newaxis]
                ctrl_var_counts = pd.Series(data=(np.where(self.categorization_result.get_orthogonal_selection((ctrls_coords, filtered_combs_coords)) > 1, 1, 0) * ctrl_adj_factors).sum(axis=0),
                                            index=filtered_combs)
                rare_categories = ctrl_var_counts[ctrl_var_counts < self.ctrl_thres].index.tolist()
                rare_category_coords = list(np.where(pd.Series(self.categories).isin(rare_categories))[0])
                
                cov = ((self.categorization_result.get_orthogonal_selection((train_samples_coords, rare_category_coords)) > 0).astype('uint8') * self.adj_factor.iloc[train_samples_coords]['AdjustFactor'].values[:,np.newaxis]).astype('float32')
                test_cov = ((self.categorization_result.get_orthogonal_selection((test_samples_coords, rare_category_coords)) > 0).astype('uint8') * self.adj_factor.iloc[test_samples_coords]['AdjustFactor'].values[:,np.newaxis]).astype('float32')
            else:
                ctrl_adj_factors = self.adj_factor.iloc[ctrls_coords]['AdjustFactor'].values[:,np.newaxis]
                ctrl_var_counts = pd.Series(data=(self.categorization_result.get_orthogonal_selection((ctrls_coords, filtered_combs_coords)) * ctrl_adj_factors).sum(axis=0),
                                            index=filtered_combs)
                rare_categories = ctrl_var_counts[ctrl_var_counts < self.ctrl_thres].index.tolist()
                rare_category_coords = list(np.where(pd.Series(self.categories).isin(rare_categories))[0])
                
                cov = (self.categorization_result.get_orthogonal_selection((train_samples_coords, rare_category_coords)) * self.adj_factor.iloc[train_samples_coords]['AdjustFactor'].values[:,np.newaxis]).astype('float32')
                test_cov = (self.categorization_result.get_orthogonal_selection((test_samples_coords, rare_category_coords)) * self.adj_factor.iloc[test_samples_coords]['AdjustFactor'].values[:,np.newaxis]).astype('float32')
        else:
            if self.use_n_carrier:
                ctrl_categorization_result = (self.categorization_result.get_orthogonal_selection((ctrls_coords, filtered_combs_coords)) > 0).astype('uint8')
                ctrl_var_counts = pd.Series(data=ctrl_categorization_result.sum(axis=0),
                                            index=filtered_combs)
                rare_categories = ctrl_var_counts[ctrl_var_counts < self.ctrl_thres].index.tolist()
                rare_category_coords = list(np.where(pd.Series(self.categories).isin(rare_categories))[0])

                cov = (self.categorization_result.get_orthogonal_selection((train_samples_coords, rare_category_coords)) > 0).astype('uint8')
                test_cov = (self.categorization_result.get_orthogonal_selection((test_samples_coords, rare_category_coords)) > 0).astype('uint8')
            else:
                ctrl_var_counts = pd.Series(data=self.categorization_result.get_orthogonal_selection((ctrls_coords, filtered_combs_coords)).sum(axis=0),
                                            index=filtered_combs)
                rare_categories = ctrl_var_counts[ctrl_var_counts < self.ctrl_thres].index.tolist()
                rare_category_coords = list(np.where(pd.Series(self.categories).isin(rare_categories))[0])

                cov = self.categorization_result.get_orthogonal_selection((train_samples_coords, rare_category_coords))
                test_cov = self.categorization_result.get_orthogonal_selection((test_samples_coords, rare_category_coords))
            
        if not swap_label:
            log.print_progress(f"# of rare categories (Seed: {seed}): {len(rare_categories)}")
        
        gc.collect()
        return rare_categories, numpy2ri.py2rpy(cov), numpy2ri.py2rpy(test_cov), response, test_response
    
    def _risk_score_per_category(self, seed: int, swap_label: bool = False, rare_categories = None, cov = None, test_cov = None, response = None, test_response = None, filtered_combs = None):
        """Lasso model selection"""
        output_dict = defaultdict(dict)
        
        if swap_label:
            rare_categories, cov, test_cov, response, test_response = self._create_covariates(seed, swap_label, filtered_combs)
            
        if len(rare_categories) == 0:
            if not swap_label:
                log.print_warn(f"There are no rare categories (Seed: {seed}).")
            return

        y = np.where(response, 1., 0.)
        test_y = np.where(test_response, 1., 0.)
        
        if not (swap_label or self.do_loop or self.n_of_one_leave):
            log.print_progress(f"Running LassoCV (Seed: {seed})")
        
        # Create a glmnet model
        cvfit = self.cv_glmnet(x = cov,
                               y = ro.FloatVector(y),
                               alpha = 1,
                               foldid = numpy2ri.py2rpy(self._custom_cv_folds(len(response), seed)),
                               type_measure='deviance',
                               nlambda=100)
        
        # Get the lambda with the minimum mean cross-validated error
        opt_lambda = cvfit.rx2("lambda.min")[0]
        # The index is 1-based.
        lambda_values = np.array(cvfit.rx2("lambda"))
        i_choose = np.where(lambda_values == opt_lambda)[0][0]
        
        # Get lasso coefficients
        full_glmnet = cvfit.rx2("glmnet.fit")
        beta_matrix = np.array(ro.r["as.matrix"](ro.r["coef"](full_glmnet, s=opt_lambda)))

        rare_idx = filtered_combs.isin(rare_categories)
        opt_coeff = np.zeros(len(rare_idx))
        
        # beta_matrix includes the intercept term, resulting in one additional column compared to your original input data.
        opt_coeff[rare_idx] = beta_matrix[1:, 0]
                
        n_select = np.sum(np.abs(opt_coeff) != 0.0)

        # Compute the predictive R-squared using the test set
        predict_function = ro.r["predict"]
        predicted_values = predict_function(full_glmnet, newx=test_cov)
        predict_y = np.array(predicted_values)[:, i_choose]
        rsq = r2_score(test_y, predict_y)

        output_dict[seed] = [opt_lambda, rsq, n_select, opt_coeff]
                
        if not (swap_label or self.do_loop or self.n_of_one_leave):
            log.print_progress(f"Done (Seed: {seed})")
                
        gc.collect()
        return output_dict    
        
    def _custom_cv_folds(self, nobs: int, seed: int = 42) -> Tuple[np.ndarray, np.ndarray]:
        """Customize k-fold cross-validation"""
        np.random.seed(seed)
        rand_idx = np.random.permutation(nobs)
        foldid = np.repeat(0, nobs)
        i=0
        
        while i<=self.fold:
            idx_val = rand_idx[np.arange(nobs * (i - 1) / self.fold, nobs * i / self.fold, dtype=int)]
            foldid[idx_val] = i
            i+=1
            
        return foldid
    
    def save_results(self):
        """Save the results to a file """
        log.print_progress(self.save_results.__doc__)

        domain_list = list(self._result_dict.keys())
        null_models = []
        fin_res = pd.DataFrame()

        for domain in domain_list:
            filtered_combs = self.filtered_category_set.loc[self.filtered_category_set['is_'+domain]==1]['Category'] if domain != 'all' else pd.Series(self.categories)

            result_table = []

            choose_idx = np.all([self._result_dict[domain][seed][3] != 0 
                                for seed in self._result_dict[domain].keys()], axis=0)
            
            ## Get the categories which are selected by the LassoCV for all seeds
            coef_df = pd.DataFrame.from_dict(
                {seed: self._result_dict[domain][seed][3][choose_idx]
                for seed in self._result_dict[domain].keys()},
                orient="index",
                columns = filtered_combs[choose_idx]
            )

            coef_df.to_csv(str(self.coef_path).replace('.txt', f'.{domain}.txt'), sep="\t")

            for seed in self._result_dict[domain].keys():
                result_table += [[domain] + [str(seed)] + self._result_dict[domain][seed][:-1]]

            result_df = pd.DataFrame(result_table, columns=["Domain", "Seed", "Parameter", "R2", "N_select"])
            new_df = pd.DataFrame([domain, 'average', result_df['Parameter'].mean(), result_df['R2'].mean(), sum(choose_idx)]).T
            new_df.columns = result_df.columns
            result_df = pd.concat([new_df, result_df], ignore_index=True)

            if not self.predict_only:
                r2_scores = np.array([self._permutation_dict[domain][seed][1]
                                      for seed in self._permutation_dict[domain].keys()])
                null_models.append([domain, 'average', r2_scores.mean(), r2_scores.std()])
                null_models.extend([[domain] + [i+1] + [str(r2_scores[i])] + [''] for i in range(len(r2_scores))])

                new_values = []
                for row in result_df['R2']:
                    new_value = (np.sum(r2_scores >= row) + 1) / (len(r2_scores) + 1)
                    new_values.append(new_value)

                result_df['Perm_P'] = new_values
                
                self.draw_histogram_plot(domain=domain,
                                         r2 = float(result_df.loc[result_df['Seed'] == 'average']['R2'].values),
                                         perm_r2 = r2_scores)
            
            fin_res = pd.concat([fin_res, result_df], ignore_index=True)

        fin_res.to_csv(self.result_path, sep="\t", index=False)
        
        if not self.predict_only:
            null_models = pd.DataFrame(null_models,
                                       columns=["Domain", "N_perm", "R2", "std"])
            null_models.to_csv(self.null_model_path, sep="\t", index=False)
  
    def draw_histogram_plot(self, domain: str, r2: float, perm_r2: np.ndarray):
        log.print_progress("Save histogram plot")
        
        # Set the font size
        plt.rcParams.update({'font.size': 8})
        
        # Set the figure size
        plt.figure(figsize=(7, 7))

        # Create the histogram plot
        plt.hist(perm_r2, bins=20, color='lightgrey', edgecolor='black')
        
        text_label1 = 'P={:.2f}'.format((sum(perm_r2>=r2)+1)/(len(perm_r2)+1))
        text_label2 = '$R^2$={:.2f}%'.format(r2*100)

        # Add labels and title
        plt.title(f'Histogram Plot (Domain: {domain})', fontsize = 8)
        plt.xlabel('$R^2$')
        plt.ylabel('Frequency')
        plt.axvline(x=r2, color='red')
        plt.text(0.05, 0.95, text_label1, transform=plt.gca().transAxes, ha='left', va='top', fontsize=8, color='black')
        plt.text(0.05, 0.85, text_label2, transform=plt.gca().transAxes, ha='left', va='top', fontsize=8, color='red')
        plt.locator_params(axis='x', nbins=5)
        plt.tight_layout()
        plt.savefig(str(self.plot_path).replace('.pdf', f'.{domain}.pdf'), bbox_inches='tight')
        
    def update_env(self):
        """Update the environment variables """
        log.print_progress(self.update_env.__doc__)
        
        self.set_env("LASSO_RESULTS", self.result_path)
        if not self.predict_only:
            self.set_env("LASSO_NULL_MODELS", self.null_model_path)
        self.save_env()

    def gather_results_for_loop(self, annotation):
        """Gather the results"""
        log.print_progress(self.gather_results_for_loop.__doc__)

        domain_list = list(self._result_dict.keys())
        fin_res = pd.DataFrame()

        for domain in domain_list:
            filtered_combs = self.filtered_category_set.loc[self.filtered_category_set['is_'+domain]==1]['Category'] if domain != 'all' else pd.Series(self.categories)

            result_table = []

            choose_idx = np.all([self._result_dict[domain][seed][3] != 0 
                                for seed in self._result_dict[domain].keys()], axis=0)
            
            ## Get the categories which are selected by the LassoCV for all seeds
            coef_df = pd.DataFrame.from_dict(
                {seed: self._result_dict[domain][seed][3][choose_idx]
                for seed in self._result_dict[domain].keys()},
                orient="index",
                columns = filtered_combs[choose_idx]
            )

            file_suffix = '' if self.do_loop else '.excluded'
            file_name = f".{domain}.{annotation}{file_suffix}.txt"
            coef_df.to_csv(str(self.coef_path).replace('.txt', file_name), sep="\t")

            for seed in self._result_dict[domain].keys():
                result_table += [[domain] + [str(seed)] + self._result_dict[domain][seed][:-1]]

            result_df = pd.DataFrame(result_table, columns=["Domain", "Seed", "Parameter", "R2", "N_select"])
            new_df = pd.DataFrame([domain, 'average', result_df['Parameter'].mean(), result_df['R2'].mean(), sum(choose_idx)]).T
            new_df.columns = ["Domain", "Seed", "Parameter", "R2", "N_select"]
            #result_df = pd.concat([new_df, result_df], ignore_index=True)

            if not self.predict_only:
                r2_scores = np.array([self._permutation_dict[domain][seed][1]
                                      for seed in self._permutation_dict[domain].keys()])
                new_df['Perm_P'] = (np.sum(r2_scores >= new_df['R2'].values[0]) + 1) / (len(r2_scores) + 1)

            fin_res = pd.concat([fin_res, new_df], ignore_index=True)

        #fin_res.to_csv(self.result_path, sep="\t", index=False)
        self._result_for_loop[annotation] = fin_res

    def save_results_for_loop(self):
        # Initialize an empty list to store the DataFrames
        dataframe_list = []
        key_name = 'Annotation' if self.do_loop else 'Annotation_excluded'
        # Loop through the dictionary and append DataFrames to the list
        for key, df in self._result_for_loop.items():
            # Add a new column 'Domain' with the key from the dictionary
            df[key_name] = key
            dataframe_list.append(df)
        fin_res = pd.concat(dataframe_list, ignore_index=True)
        # Move 'key_name' column to the first column
        column_order = [key_name] + [col for col in fin_res.columns if col != key_name]
        fin_res = fin_res[column_order]

        file_suffix = 'each_annot_loop' if self.do_loop else 'n_of_one_leave'
        file_name = f".{file_suffix}.txt"
        fin_res.to_csv(str(self.result_path).replace('.txt', file_name), sep="\t", index=False)
