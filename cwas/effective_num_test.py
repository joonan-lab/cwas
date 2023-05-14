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
        self._category_set_path = None
        self._category_set = None
        self._tag = None
        self._num_sim = None
        
    @staticmethod
    def _print_args(args: argparse.Namespace):
        print_arg("Concatenated z-scores", args.zscore_df_path)
        print_arg("Output directory", args.output_dir_path)
        if args.tag:
            print_arg("Output tag (prefix of output files)", args.tag)
        if args.category_set_path :
            print_arg("Category set file", args.category_set_path)

    @staticmethod
    def _check_args_validity(args: argparse.Namespace):
        check_is_file(args.zscore_df_path)
        check_is_dir(args.output_dir_path)
        if args.category_set_path :
            check_is_file(args.category_set_path)

    @property
    def zscore_df_path(self) -> Path:
        return self.args.zscore_df_path.resolve()
    
    @property
    def zscore_df(self) -> pd.DataFrame:
        if self._zscore_df is None:
            self._zscore_df = pd.read_table(self.zscore_df_path, index_col='Simulation')
        return self._zscore_df    

    @property
    def num_sim(self) -> int:
        if self._num_sim is None:
            self._num_sim = self.zscore_df.shape[0]
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
    def corr_mat_path(self) -> Path:
        return Path(
            str(self.zscore_df_path).replace('.zscores.txt.gz', f'.corr_mat_{self.tag}.pickle')
        )

    @property
    def neg_lap_path(self) -> Path:
        return Path(
            str(self.zscore_df_path).replace('.zscores.txt.gz', f'.neg_lap_{self.tag}.pickle')
        )

    @property
    def eig_val_path(self) -> Path:
        return Path(
            str(self.zscore_df_path).replace('.zscores.txt.gz', f'.eig_vals_{self.tag}.pickle')
        )

    @property
    def eig_vec_path(self) -> Path:
        return Path(
            str(self.zscore_df_path).replace('.zscores', f'.eig_vecs_{self.tag}')
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
        
        self.eff_num_test = eff_num_test
        
    def eign_decomposition(self, save_vecs: bool = True):
        print_progress(f"Calculate eign values")      
        if self.category_set_path:
            filtered_combs = self.category_set["Category"]
        else:
            if self.zscore_df.columns[0] == 'Simulation':
                filtered_combs = self.zscore_df.columns[1:]
            else:
                filtered_combs = self.zscore_df.columns
        
        if not os.path.isfile(self.corr_mat_path):
            filtered_zscore_df = self.zscore_df[filtered_combs]

            # Convert infinite Z to finite Z
            if (np.isinf(filtered_zscore_df).any().sum()>0):
                print_warn("The z-score matrix contains infinite Z. Infinite Z will be replaced with finite Z.")
                for_inf_z = norm.ppf(1 - 0.9999999999999999)
                filtered_zscore_df = filtered_zscore_df.replace([-np.inf], for_inf_z)            
            
            corr_mat = np.corrcoef(filtered_zscore_df.values.T)
            if np.isnan(corr_mat).any():
                print_warn("The correlation matrix contains NaN. NaN will be replaced with 0,1.")
                for i in range(corr_mat.shape[0]):
                    if np.isnan(corr_mat[i, i]):
                        corr_mat[i, i] = 1.0
                np.nan_to_num(corr_mat, copy=False)
                
            print_progress("Writing the correlation matrix to file")
            pickle.dump(corr_mat, open(self.corr_mat_path, 'wb'), protocol=5)
        else:
            with self.corr_mat_path.open('rb') as f:
                corr_mat = pickle.load(f)

        if not os.path.isfile(self.neg_lap_path):
            print_progress("Generating the negative laplacian matrix")
            neg_lap = np.abs(corr_mat)
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
        self.set_env("N_EFFECTIVE_TEST", self.eff_num_test)
        self.save_env()
    
    