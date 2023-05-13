import argparse, os, sys, pickle
from pathlib import Path
import numpy as np
import pandas as pd

from cwas.utils.log import print_progress, print_arg, print_warn, print_log
from cwas.runnable import Runnable

class EffNumTest(Runnable):
    def __init__(self, args: argparse.Namespace):
        super().__init__(args)
        self._in_vcf = None
        self._sample_info = None
        
    
    
    
    
    def get_n_etests(self):
        """Get the number of effective tests """
        print_progress(self.get_n_etests.__doc__)
        if self.eig_val_path("all").is_file():
            print_log(
                "NOTICE",
                "You already have eign values for all categories.",
                False,
            )
        else:
            self.eign_decomposition("all", save_vecs=False)
            
        with self.eig_val_path("all").open('rb') as f:
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
        
        log.print_log("RESULT", f"The number of effective tests is {eff_num_test}.", False)
        
        self.eff_num_test = eff_num_test
        
    def eign_decomposition(self, target_term: str, save_vecs: bool = True):
        print_progress(f"Calculate eign values for {target_term} ")      
        filtered_combs = self.target_combinations[target_term]
        
        if not self.corr_mat_path(target_term).is_file():
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
            pickle.dump(corr_mat, open(self.corr_mat_path(target_term), 'wb'), protocol=5)
        else:
            with self.corr_mat_path(target_term).open('rb') as f:
                corr_mat = pickle.load(f)

        if not self.neg_lap_path(target_term).is_file():
            print_progress("Generating the negative laplacian matrix")
            neg_lap = np.abs(corr_mat)
            degrees = np.sum(neg_lap, axis=0)
            for i in tqdm(range(neg_lap.shape[0])):
                neg_lap[i, :] = neg_lap[i, :] / np.sqrt(degrees)
                neg_lap[:, i] = neg_lap[:, i] / np.sqrt(degrees)
            print_progress("Writing the negative laplacian matrix to file")
            pickle.dump(neg_lap, open(self.neg_lap_path(target_term), 'wb'), protocol=5)
        else:
            with self.neg_lap_path(target_term).open('rb') as f:
                neg_lap = pickle.load(f)

        if not self.eig_val_path(target_term).is_file() or not self.eig_vec_path(target_term).is_file():
            print_progress("Calculating the eigenvalues of the negative laplacian matrix")
            eig_vals, eig_vecs = np.linalg.eig(neg_lap)
            print_progress("Writing the eigenvalues to file")
            pickle.dump(eig_vals, open(self.eig_val_path(target_term), 'wb'), protocol=5)
            if save_vecs:
                print_progress("Writing the eigenvectors to file")
                pd.DataFrame(eig_vecs.real, index=filtered_combs).to_csv(self.eig_vec_path(target_term), sep='\t', index=True, header=False)
                
                
                