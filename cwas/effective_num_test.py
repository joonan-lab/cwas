import argparse, os, sys, pickle
from pathlib import Path
import numpy as np
import pandas as pd
from tqdm import tqdm
from typing import Optional

from cwas.utils.log import print_progress, print_arg, print_warn, print_log
from cwas.runnable import Runnable
from scipy.stats import norm
from cwas.utils.check import check_is_file, check_is_dir

class EffectiveNumTest(Runnable):
    def __init__(self, args: argparse.Namespace):
        super().__init__(args)
        self._zscore_df = None
        self._intersection_matrix = None
        self._correlation_matrix = None
        self._category_set_path = None
        self._category_set = None
        self._tag = None
        self._num_sim = None
        self.eff_num_test_value = None
        self._replace_term = None
        self._binom_p = None
        self._sample_info = None
        
    @staticmethod
    def _print_args(args: argparse.Namespace):
        print_arg("Input file", args.input_path)
        print_arg("Input format", args.input_format)
        print_arg("Output directory", args.output_dir_path)
        if args.sample_info_path:
            print_arg("Sample information file", args.sample_info_path)
        if args.tag:
            print_arg("Output tag (prefix of output files)", args.tag)
        if args.category_set_path :
            print_arg("Category set file", args.category_set_path)

    @staticmethod
    def _check_args_validity(args: argparse.Namespace):
        check_is_file(args.input_path)
        check_is_dir(args.output_dir_path)
        if args.sample_info_path:
            check_is_file(args.sample_info_path)
        if args.category_set_path :
            check_is_file(args.category_set_path)

    @property
    def input_path(self) -> Path:
        return self.args.input_path.resolve()

    @property
    def input_format(self) -> str:
        return self.args.input_format
    
    @property
    def zscore_df(self) -> pd.DataFrame:
        if self._zscore_df is None and self.input_format == 'zscores':
            self._zscore_df = pd.read_table(self.input_path, index_col='Simulation')
        return self._zscore_df

    @property
    def intersection_matrix(self) -> pd.DataFrame:
        if self._intersection_matrix is None and self.input_format == 'inter':
            with self.input_path.open('rb') as f:
                self._intersection_matrix = pickle.load(f)
            self._intersection_matrix.index = list(map(str, self._intersection_matrix.index.tolist()))
            self._intersection_matrix.columns = list(map(str, self._intersection_matrix.columns.tolist()))
        return self._intersection_matrix

    @property
    def correlation_matrix(self) -> pd.DataFrame:
        if self._correlation_matrix is None and self.input_format == 'corr':
            with self.input_path.open('rb') as f:
                self._correlation_matrix = pickle.load(f)
            self._correlation_matrix.index = list(map(str, self._correlation_matrix.index.tolist()))
            self._correlation_matrix.columns = list(map(str, self._correlation_matrix.columns.tolist()))
        return self._correlation_matrix

    @property
    def sample_info_path(self) -> Optional[Path]:
        return self.args.sample_info_path.resolve()

    @property
    def sample_info(self) -> pd.DataFrame:
        if self._sample_info is None:
            self._sample_info = pd.read_table(
                self.sample_info_path, index_col="SAMPLE"
            )
        return self._sample_info

    @property
    def binom_p(self) -> float:
        if self._binom_p is None:
            self._binom_p = (self.sample_info["PHENOTYPE"] == "case").sum() / np.isin(self.sample_info["PHENOTYPE"], ["case", "ctrl"]).sum()
        return self._binom_p

    @property
    def num_sim(self) -> int:
        if self._num_sim is None:
            self._num_sim = self.args.num_sim
            return self._num_sim

    @property
    def category_set_path(self) -> Optional[Path]:
        return (
            self.args.category_set_path.resolve()
            if self.args.category_set_path
            else None
        )
    @property
    def category_set(self) -> pd.DataFrame:
        if self._category_set is None and self.category_set_path:
            self._category_set = pd.read_csv(self.category_set_path, sep='\t')
        return self._category_set

    @property
    def eff_num_test(self) -> bool:
        return self.args.eff_num_test

    @property
    def tag(self) -> str:
        return self.args.tag

    @property
    def replace_term(self) -> str:
        if self._replace_term is None:
            if '.zscores.txt.gz' in str(self.input_path):
                self._replace_term = '.zscores.txt.gz'
            elif 'correlation_matrix.pkl' in str(self.input_path):
                self._replace_term = '.correlation_matrix.pkl'
            elif 'intersection_matrix.pkl' in str(self.input_path):
                self._replace_term = '.intersection_matrix.pkl'
        return self._replace_term

    @property
    def corr_mat_path(self) -> Path:
        return Path(
            str(self.input_path).replace(self.replace_term, f'.correlation_matrix_{self.tag}.pickle')
        )

    @property
    def neg_lap_path(self) -> Path:
        return Path(
            str(self.input_path).replace(self.replace_term, f'.neg_lap_{self.tag}.pickle')
        )

    @property
    def eig_val_path(self) -> Path:
        return Path(
            str(self.input_path).replace(self.replace_term, f'.eig_vals_{self.tag}.pickle')
        )

    @property
    def eig_vec_path(self) -> Path:
        return Path(
            str(self.input_path).replace(self.replace_term, f'.eig_vecs_{self.tag}.txt.gz')
        )

    def run(self):
        print_arg("Number of simulations", self.num_sim)
        if self.eff_num_test:
            self.eign_decomposition()
            self.get_n_etests()
            self.update_env()
        else:
            self.eign_decomposition()
        print_progress("Done")
    
    def get_n_etests(self):
        """Get the number of effective tests """
        print_progress(self.get_n_etests.__doc__)
        if os.path.isfile(self.eig_val_path):
            print_log(
                "NOTICE",
                "You already have eign values for the selected categories.",
                False,
            )
        else:
            self.eign_decomposition(save_vecs=False)
            
        with self.eig_val_path.open('rb') as f:
            eig_vals = pickle.load(f)
        
        e = 1e-12
        eig_vals = sorted(eig_vals, key=np.linalg.norm, reverse=True)
        num_eig_val = self.num_sim
        clean_eig_vals = np.array(eig_vals[:num_eig_val])
        clean_eig_vals = clean_eig_vals[clean_eig_vals >= e]
        clean_eig_val_total_sum = np.sum(clean_eig_vals)
        clean_eig_val_sum = 0
        eff_num_test = 0
        
        for i in range(len(clean_eig_vals)):
            clean_eig_val_sum += clean_eig_vals[i]

            if clean_eig_val_sum / clean_eig_val_total_sum >= 0.99:
                eff_num_test = i + 1
                break
        
        print_log("RESULT", f"The number of effective tests is {eff_num_test}.", False)
        
        self.eff_num_test_value = eff_num_test
        
    def eign_decomposition(self, save_vecs: bool = True):
        print_progress(f"Calculate eign values")
        if self.category_set_path:
            filtered_combs = self.category_set["Category"]
        else:
            if self.input_format == 'corr':
                filtered_combs = self.correlation_matrix.columns.tolist()
            elif self.input_format == 'inter':
                filtered_combs = self.intersection_matrix.columns.tolist()
            elif self.input_format == 'zscores':
                if self.zscore_df.columns[0] == 'Simulation':
                    filtered_combs = self.zscore_df.columns[1:].tolist()
                else:
                    filtered_combs = self.zscore_df.columns.tolist()

        if self.input_format == 'corr':
            intermediate_mat = self.correlation_matrix.loc[filtered_combs,filtered_combs]
        elif self.input_format == 'inter':
            print_progress("Generating a covariance matrix")
            intermediate_mat = self.intersection_matrix.loc[filtered_combs,filtered_combs]
            intermediate_mat = intermediate_mat.mul((self.binom_p)*(1-self.binom_p))
        elif self.input_format == 'zscores':
            if not os.path.isfile(self.corr_mat_path):
                filtered_zscore_df = self.zscore_df[filtered_combs]
            
                # Convert infinite Z to finite Z
                if (np.isinf(filtered_zscore_df).any().sum()>0):
                    print_warn("The z-score matrix contains infinite Z. Infinite Z will be replaced with finite Z.")
                    for_inf_z = norm.ppf(1 - 0.9999999999999999)
                    filtered_zscore_df = filtered_zscore_df.replace([-np.inf], for_inf_z)            
                
                intermediate_mat = np.corrcoef(filtered_zscore_df.values.T)
                if np.isnan(intermediate_mat).any():
                    print_warn("The correlation matrix contains NaN. NaN will be replaced with 0,1.")
                    for i in range(intermediate_mat.shape[0]):
                        if np.isnan(intermediate_mat[i, i]):
                            intermediate_mat[i, i] = 1.0
                    np.nan_to_num(intermediate_mat, copy=False)
                    
                corr_mat2 = pd.DataFrame(intermediate_mat, columns=filtered_zscore_df.columns, index=filtered_zscore_df.columns)
                print_progress("Writing the correlation matrix to file")
                pickle.dump(corr_mat2, open(self.corr_mat_path, 'wb'), protocol=5)
            else:
                with self.corr_mat_path.open('rb') as f:
                    intermediate_mat = pickle.load(f)

        if not os.path.isfile(self.neg_lap_path):
            print_progress("Generating the negative laplacian matrix")
            if self.input_format == 'zscores':
                neg_lap = np.abs(intermediate_mat)
            else:
                neg_lap = np.abs(intermediate_mat.values)
            degrees = np.sum(neg_lap, axis=0)
            for i in tqdm(range(neg_lap.shape[0])):
                neg_lap[i, :] = neg_lap[i, :] / np.sqrt(degrees)
                neg_lap[:, i] = neg_lap[:, i] / np.sqrt(degrees)
            print_progress("Writing the negative laplacian matrix to file")
            pickle.dump(neg_lap, open(self.neg_lap_path, 'wb'), protocol=5)
        else:
            with self.neg_lap_path.open('rb') as f:
                neg_lap = pickle.load(f)

        if not os.path.isfile(self.eig_val_path) or not os.path.isfile(self.eig_vec_path):
            print_progress("Calculating the eigenvalues of the negative laplacian matrix")
            eig_vals, eig_vecs = np.linalg.eig(neg_lap)
            print_progress("Writing the eigenvalues to file")
            pickle.dump(eig_vals, open(self.eig_val_path, 'wb'), protocol=5)
            if save_vecs:
                print_progress("Writing the eigenvectors to file")
                pd.DataFrame(eig_vecs.real, index=filtered_combs).to_csv(self.eig_vec_path, sep='\t', index=True, header=False)
                
                
    def update_env(self):
        self.set_env("N_EFFECTIVE_TEST", self.eff_num_test_value)
        self.save_env()
    
    