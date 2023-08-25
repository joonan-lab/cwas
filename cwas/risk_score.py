import os, sys, argparse
import pandas as pd
import numpy as np
import cwas.utils.log as log
from pathlib import Path
from tqdm import tqdm
#from glmnet import ElasticNet, LogitNet
from sklearn.linear_model import ElasticNetCV
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import r2_score
from cwas.core.common import cmp_two_arr
from cwas.utils.check import check_is_file, check_num_proc
from cwas.runnable import Runnable
from typing import Optional, Tuple
from contextlib import contextmanager
from collections import defaultdict
import matplotlib.pyplot as plt
import polars as pl
import re

class RiskScore(Runnable):
    def __init__(self, args: argparse.Namespace):
        super().__init__(args)
        self._sample_info = None
        self._categorization_result = None
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
        if args.tag:
            log.print_arg("Output tag (prefix of output files)", args.tag)
        log.print_arg("If the number of carriers is used for calculating R2 or not", args.use_n_carrier)
        log.print_arg(
            "Threshold for selecting rare categories",
            f"{args.ctrl_thres: ,d}",
        )
        log.print_arg("Fraction of the training set", f"{args.train_set_f: ,f}")
        log.print_arg(
            "No. regression trials to calculate a mean of R squares",
            f"{args.num_reg: ,d}",
        )
        log.print_arg(
            "No. folds for CV",
            f"{args.fold: ,d}",
        )
        #log.print_arg("Use Logistic regression", args.logistic)
        log.print_arg(
            "No. permutation used to calculate the p-value",
            f"{args.n_permute: ,d}",
        )
        log.print_arg("Skip the permutation test", args.predict_only)
        log.print_arg(
            "No. worker processes for permutation",
            f"{args.num_proc: ,d}",
        )

    @staticmethod
    def _check_args_validity(args: argparse.Namespace):
        check_is_file(args.sample_info_path)
        check_num_proc(args.num_proc)
        if args.categorization_result_path:
            check_is_file(args.categorization_result_path)
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
    def adj_factor_path(self) -> Optional[Path]:
        return (
            self.args.adj_factor_path.resolve()
            if self.args.adj_factor_path
            else None
        )

    @property
    def plot_path(self) -> Optional[Path]:
        tag = '' if self.tag is None else ''.join([self.tag, '_'])
        f_name = re.sub(r'.categorization_result\.txt\.gz|.categorization_result\.txt', f'.lasso_histogram_{tag}thres_{self.ctrl_thres}.pdf', str(self.categorization_result_path.name))
        return Path(
            f"{self.out_dir}/" +
            f"{f_name}"
        )

    @property
    def category_set_path(self) -> Optional[Path]:
        return (
            self.args.category_set_path.resolve()
            if self.args.category_set_path
            else None
        )

    #@property
    #def logistic(self) -> bool:
    #    return self.args.logistic

    @property
    def out_dir(self) -> Path:
        return(self.args.output_dir_path.resolve())

    @property
    def category_set(self) -> pd.DataFrame:
        if self._category_set is None and self.category_set_path:
            self._category_set = pd.read_csv(self.category_set_path, sep='\t')
        return self._category_set

    @property
    def sample_info_path(self) -> Path:
        return self.args.sample_info_path.resolve()

    @property
    def use_n_carrier(self) -> bool:
        return self.args.use_n_carrier

    @property
    def domain_list(self) -> str:
        if self.args.domain_list=='run_all':
            return 'all,coding,noncoding,ptv,missense,damaging_missense,promoter,noncoding_wo_promoter,intron,intergenic,utr,lincRNA'
        else:
            return self.args.domain_list
        

    @property
    def tag(self) -> str:
        return self.args.tag

    @property
    def num_proc(self) -> int:
        return self.args.num_proc

    @property
    def num_reg(self) -> int:
        return self.args.num_reg

    @property
    def train_set_f(self) -> float:
        return self.args.train_set_f

    @property
    def fold(self) -> int:
        return self.args.fold

    @property
    def n_permute(self) -> int:
        return self.args.n_permute
    
    @property
    def ctrl_thres(self) -> int:
        return self.args.ctrl_thres
    
    @property
    def predict_only(self) -> bool:
        return self.args.predict_only
    
    @property
    def seed(self) -> int:
        return self.args.seed
    
    @property
    def sample_info(self) -> pd.DataFrame:
        if self._sample_info is None:
            self._sample_info = pd.read_table(
                self.sample_info_path, index_col="SAMPLE", dtype={"SAMPLE": str}, sep="\t"
            )
            if ("SET" not in self._sample_info.columns):
                log.print_log("LOG", 
                              "No 'SET' column in sample information file. "
                              "Training and test sets will be assigned randomly.")
                
                case_count = sum(self._sample_info['PHENOTYPE'] == 'case')
                ctrl_count = sum(self._sample_info['PHENOTYPE'] == 'ctrl')
                
                self.case_f = round(case_count*self.train_set_f)
                self.ctrl_f = round(ctrl_count*self.train_set_f)
                
                
                log.print_log("LOG",
                              "Use {} cases and {} controls for training set".format(self.case_f, self.ctrl_f)
                              )
                              #"Use {} samples from each phenotype as training set."
                              #.format(self.min_size))
                
                
                case_test_idx = self._sample_info.loc[self._sample_info.PHENOTYPE=='case'].sample(n=self.case_f, random_state=42).index
                ctrl_test_idx = self._sample_info.loc[self._sample_info.PHENOTYPE=='ctrl'].sample(n=self.ctrl_f, random_state=42).index
                self._sample_info.loc[case_test_idx, 'SET'] = 'training'
                self._sample_info.loc[ctrl_test_idx, 'SET'] = 'training'
                self._sample_info["SET"] = self._sample_info['SET'].fillna('test')
                
                #self.min_size = int(np.rint(min(case_count, ctrl_count) * self.train_set_f))
                #test_idx = self._sample_info.groupby('PHENOTYPE').sample(n=self.min_size, random_state=42).index
                #self._sample_info["SET"] = np.where(self._sample_info.index.isin(test_idx), "test", "training")
        return self._sample_info

    @property
    def adj_factor(self) -> pd.DataFrame:
        if self._adj_factor is None and self.adj_factor_path:
            self._adj_factor = pd.read_table(
                self.adj_factor_path, index_col="SAMPLE", dtype={"SAMPLE": str}, sep="\t"
            )
        return self._adj_factor
    
    @property
    def categorization_result(self) -> pd.DataFrame:
        if self._categorization_result is None:
            log.print_progress("Load the categorization result")
            self._categorization_result = pl.read_csv(
                self.categorization_result_path, dtypes={"SAMPLE": str}, separator = '\t'
            )
            self._categorization_result = self._categorization_result.to_pandas().set_index("SAMPLE")
            if self.adj_factor is not None:
                self._adjust_categorization_result()
            if self.use_n_carrier:
                self._categorization_result = self._categorization_result.applymap(lambda x: 1 if x > 0 else 0)
            log.print_log("LOG",
                          "Categorization result has {} samples and {} categories."
                          .format(self._categorization_result.shape[0], self._categorization_result.shape[1]))
        return self._categorization_result
    
    def _adjust_categorization_result(self):
        if not self._contain_same_index(
           self._categorization_result, self._adj_factor
        ):
            raise ValueError(
                "The sample IDs from the adjustment factor list are "
                "not the same with the sample IDs "
                "from the categorization result."
            )
        '''
        adj_factors = [
            self.adj_factor.to_dict()["AdjustFactor"][sample_id]
            for sample_id in self._categorization_result.index.values
        ]
        '''
        adj_factors = self.adj_factor.loc[self._categorization_result.index.values, "AdjustFactor"].tolist()
        self._categorization_result = self._categorization_result.multiply(
            adj_factors, axis="index"
        )

    @staticmethod
    def _contain_same_index(table1: pd.DataFrame, table2: pd.DataFrame) -> bool:
        return cmp_two_arr(table1.index.values, table2.index.values)
    
    @property
    def datasets(self) -> np.ndarray:
        if self._datasets is None:
            sample_ids = self.categorization_result.index.values
            set_dict = self.sample_info.to_dict()["SET"]
            self._datasets = np.array([set_dict[sample_id] for sample_id in sample_ids])
        return self._datasets

    @property
    def covariates(self) -> pd.DataFrame:
        if self._covariates is None:
            self._covariates = self.categorization_result[self.datasets == "training"]
        return self._covariates
    
    @property
    def test_covariates(self) -> pd.DataFrame:
        if self._test_covariates is None:
            self._test_covariates = self.categorization_result[self.datasets == "test"]
        return self._test_covariates
    
    @property
    def response(self) -> np.ndarray:
        if self._response is None:
            sample_ids = self.categorization_result[self.datasets == "training"].index.values
            phenotype_dict = self.sample_info.to_dict()["PHENOTYPE"]
            self._response = np.array([phenotype_dict[sample_id] == "case" for sample_id in sample_ids])
        return self._response

    @property
    def test_response(self) -> np.ndarray:
        if self._test_response is None:
            sample_ids = self.categorization_result[self.datasets == "test"].index.values
            phenotype_dict = self.sample_info.to_dict()["PHENOTYPE"]
            self._test_response = np.array([phenotype_dict[sample_id] == "case" for sample_id in sample_ids])
        return self._test_response

    @property
    def coef_path(self) -> Path:
        tag = '' if self.tag is None else ''.join([self.tag, '_'])
        f_name = re.sub(r'.categorization_result\.txt\.gz|.categorization_result\.txt', f'.lasso_coef_{tag}thres_{self.ctrl_thres}.txt', str(self.categorization_result_path.name))
        return Path(
            f"{self.out_dir}/" +
            f"{f_name}"
        )
    
    @property
    def result_path(self) -> Path:
        tag = '' if self.tag is None else ''.join([self.tag, '_'])
        f_name = re.sub(r'.categorization_result\.txt\.gz|.categorization_result\.txt', f'.lasso_results_{tag}thres_{self.ctrl_thres}.txt', str(self.categorization_result_path.name))
        return Path(
            f"{self.out_dir}/" +
            f"{f_name}"
        )

    @property
    def null_model_path(self) -> Path:
        tag = '' if self.tag is None else ''.join([self.tag, '_'])
        f_name = re.sub(r'.categorization_result\.txt\.gz|.categorization_result\.txt', f'.lasso_null_models_{tag}thres_{self.ctrl_thres}.txt', str(self.categorization_result_path.name))
        return Path(
            f"{self.out_dir}/" +
            f"{f_name}"
        )

    def run(self):
        self.prepare()
        self.risk_scores()
        if not self.predict_only:
            self.permute_pvalues()
        self.save_results()
        self.update_env()
        log.print_progress("Done")

    def prepare(self):
        if not self._contain_same_index(
            self.categorization_result, self.sample_info
        ):
            raise ValueError(
                "The sample IDs from the sample information are "
                "not the same with the sample IDs "
                "from the categorization result."
            )
            
    def risk_scores(self):
        """Generate risk scores for various seeds """
        log.print_progress(self.risk_scores.__doc__)

        domain_values = [domain.strip() for domain in self.domain_list.split(',')]
        for domain in domain_values:
            if domain == 'all':
                filtered_combs = pd.Series(self.categorization_result.columns)
            else:
                filtered_combs = self.category_set.loc[self.category_set['is_'+domain]==1]['Category']
            
            seeds = np.arange(self.seed, self.seed + self.num_reg * 10, 10)
            for seed in seeds:
                self.risk_score_per_category(domain = domain, result_dict=self._result_dict, seed=seed, filtered_combs=filtered_combs)
            
    def risk_score_per_category(self, domain: str, result_dict: defaultdict, seed: int = 42, swap_label: bool = False, filtered_combs = None):
        """Lasso model selection """
        if swap_label:
            np.random.seed(seed)
            
            perm_sample_info = self.sample_info
            perm_sample_info['Perm_PHENOTYPE'] = np.random.permutation(self.sample_info['PHENOTYPE'])
            perm_sample_info = perm_sample_info.drop(columns=['PHENOTYPE', 'SET'])
            
            perm_case_test_idx = perm_sample_info.loc[perm_sample_info.Perm_PHENOTYPE=='case'].sample(n=self.case_f, random_state=seed).index
            perm_ctrl_test_idx = perm_sample_info.loc[perm_sample_info.Perm_PHENOTYPE=='ctrl'].sample(n=self.ctrl_f, random_state=seed).index
            perm_sample_info.loc[perm_case_test_idx, 'Perm_SET'] = 'training'
            perm_sample_info.loc[perm_ctrl_test_idx, 'Perm_SET'] = 'training'
            perm_sample_info["Perm_SET"] = perm_sample_info['Perm_SET'].fillna('test')

            #test_idx = perm_sample_info.groupby('Perm_PHENOTYPE').sample(n=self.min_size, random_state=seed).index
            #perm_sample_info["Perm_SET"] = np.where(perm_sample_info.index.isin(test_idx), "test", "training")
            
            sample_ids1 = self.categorization_result.index.values
            set_dict = perm_sample_info.to_dict()["Perm_SET"]
            datasets = np.array([set_dict[sample_id] for sample_id in sample_ids1])
            
            sample_ids2 = self.categorization_result[datasets == "test"].index.values
            phenotype_dict = perm_sample_info.to_dict()["Perm_PHENOTYPE"]
            test_response = np.array([phenotype_dict[sample_id] == "case" for sample_id in sample_ids2])
            
            sample_ids3 = self.categorization_result[datasets == "training"].index.values
            response = np.array([phenotype_dict[sample_id] == "case" for sample_id in sample_ids3])
            
            covariates = self.categorization_result[datasets == "training"]
            test_covariates = self.categorization_result[datasets == "test"]
            
            ctrl_var_counts = pd.concat(
                [covariates[filtered_combs][~response],
                 test_covariates[filtered_combs][~test_response]],
            ).sum()
            rare_idx = (ctrl_var_counts < self.ctrl_thres).values
            log.print_progress(f"# of rare categories (Seed: {seed}): {rare_idx.sum()}")
            cov = covariates[filtered_combs].iloc[:, rare_idx]
            test_cov = test_covariates[filtered_combs].iloc[:, rare_idx]
        else:
            response, test_response = self.response, self.test_response
            
            ctrl_var_counts = pd.concat(
                [self.covariates[filtered_combs][~response],
                 self.test_covariates[filtered_combs][~test_response]],
            ).sum()
            rare_idx = (ctrl_var_counts < self.ctrl_thres).values
            log.print_progress(f"# of rare categories (Seed: {seed}): {rare_idx.sum()}")
            cov = self.covariates[filtered_combs].iloc[:, rare_idx]
            test_cov = self.test_covariates[filtered_combs].iloc[:, rare_idx]

        if cov.shape[1] == 0:
            log.print_warn(f"There are no rare categories (Seed: {seed}).")
            return
        
        scaler = StandardScaler().fit(cov)
        cov2 = scaler.transform(cov)
        test_cov2 = scaler.transform(test_cov)

        y = np.where(response, 1.0, 0.0)
        test_y = np.where(test_response, 1.0, 0.0)
        log.print_progress(f"Running LassoCV (Seed: {seed})")
        
        #lasso_model = ElasticNet(alpha=1, n_lambda=100, standardize=True, n_splits=self.fold, n_jobs=self.num_proc,
        #                         scoring='mean_squared_error', random_state=seed)
        
        lasso_model = ElasticNetCV(l1_ratio=1, cv = self.custom_cv_folds(seed=seed), n_jobs = self.num_proc,
                                   random_state = seed, verbose = False, n_alphas=100, selection='random')

        lasso_model.fit(cov2, y)
        #opt_model_idx = np.argmax(getattr(lasso_model, 'cv_mean_score_'))
        #coeffs = getattr(lasso_model, 'coef_path_')
        coeffs = getattr(lasso_model, 'coef_')
        opt_coeff = np.zeros(len(rare_idx))
        opt_coeff[rare_idx] = coeffs
        
        #opt_lambda = getattr(lasso_model, 'lambda_max_')
        opt_lambda = getattr(lasso_model, 'alpha_')
        n_select = np.sum(np.abs(opt_coeff) != 0.0)
        y_pred = lasso_model.predict(test_cov2)
        rsq = r2_score(test_y, y_pred)
        result_dict[domain][seed] = [opt_lambda, rsq, n_select, opt_coeff]
        
        log.print_progress("Done")
            
    def custom_cv_folds(self, seed: int = 42) -> Tuple[np.ndarray, np.ndarray]:
        np.random.seed(seed)
        nobs = np.sum(self.datasets == "training")
        rand_idx = np.random.permutation(nobs)
        foldid = np.repeat(0, nobs)
        i=0
        while i<=self.fold:
            idx_val = rand_idx[np.arange(nobs * (i - 1) / self.fold, nobs * i / self.fold, dtype=int)]
            foldid[idx_val] = i
            i+=1
    
    def permute_pvalues(self):
        """Run LassoCV to get permutated pvalues"""
        log.print_progress(self.permute_pvalues.__doc__)
        
        @contextmanager
        def nullify_output(suppress_stdout: bool=True, suppress_stderr: bool=True):
            stdout = sys.stdout
            stderr = sys.stderr
            devnull = open(os.devnull, "w")
            try:
                if suppress_stdout:
                    sys.stdout = devnull
                if suppress_stderr:
                    sys.stderr = devnull
                yield
            finally:
                if suppress_stdout:
                    sys.stdout = stdout
                if suppress_stderr:
                    sys.stderr = stderr
                    
        seeds = np.arange(self.seed, self.seed+self.n_permute)

        domain_values = [domain.strip() for domain in self.domain_list.split(',')]
        for domain in domain_values:
            if domain == 'all':
                filtered_combs = pd.Series(self.categorization_result.columns)
            else:
                filtered_combs = self.category_set.loc[self.category_set['is_'+domain]==1]['Category']

            for seed in tqdm(seeds):
                with nullify_output():
                    self.risk_score_per_category(domain = domain, result_dict=self._permutation_dict, seed=seed, swap_label=True, filtered_combs=filtered_combs)
                

    def save_results(self):
        """Save the results to a file """
        log.print_progress(self.save_results.__doc__)

        domain_list = list(self._result_dict.keys())
        null_models = []
        fin_res = pd.DataFrame()

        for domain in domain_list:
            if domain == 'all':
                filtered_combs = pd.Series(self.categorization_result.columns)
            else:
                filtered_combs = self.category_set.loc[self.category_set['is_'+domain]==1]['Category']

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

            # Create the initial DataFrame with the 'Domain' column
            #domain_column_data = {'Domain': ['noncoding'] * coef_df.shape[0]}
            #domain_df = pd.DataFrame(domain_column_data)
            #domain_df

            coef_df.to_csv(
                str(self.coef_path).replace('.csv', f'.{domain}.csv'),
                sep="\t"
            )

            for seed in self._result_dict[domain].keys():
                result_table += [[domain] + [str(seed)] + self._result_dict[domain][seed][:-1]]
                #[opt_lambda, rsq, n_select, opt_coeff]
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
            null_models = pd.DataFrame(
                null_models,
                columns=["Domain", "N_perm", "R2", "std"]
            )
            null_models.to_csv(self.null_model_path, sep="\t", index=False)


    def update_env(self):
        """Update the environment variables """
        log.print_progress(self.update_env.__doc__)
        
        self.set_env("LASSO_RESULTS", self.result_path)
        if not self.predict_only:
            self.set_env("LASSO_NULL_MODELS", self.null_model_path)
        self.save_env()
        
    def draw_histogram_plot(self, domain: str, r2: float, perm_r2: np.ndarray):
        log.print_progress("Save histogram plot")
        
        # Set the font size
        plt.rcParams.update({'font.size': 8})
        
        # Set the figure size
        plt.figure(figsize=(4, 4))

        # Create the histogram plot
        plt.hist(perm_r2, bins=20, color='lightgrey', edgecolor='black')
        
        text_label1 = 'P={:.2f}'.format((sum(perm_r2>=r2)+1)/(len(perm_r2)+1))
        text_label2 = '$R^2$={:.2f}%'.format(r2*100)

        # Add labels and title
        plt.xlabel('$R^2$')
        plt.ylabel('Frequency')
        plt.title('Histogram Plot', fontsize = 8)
        plt.axvline(x=r2, color='red')
        plt.text(0.05, 0.95, text_label1, transform=plt.gca().transAxes, ha='left', va='top', fontsize=8, color='black')
        plt.text(0.05, 0.85, text_label2, transform=plt.gca().transAxes, ha='left', va='top', fontsize=8, color='red')
        plt.locator_params(axis='x', nbins=5)

        plt.savefig(str(self.plot_path).replace('.pdf', f'.{domain}.pdf'))
        



