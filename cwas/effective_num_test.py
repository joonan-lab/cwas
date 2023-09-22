import argparse, os, sys, pickle
from pathlib import Path
import numpy as np
import pandas as pd
from tqdm import tqdm
from typing import Optional
import re

from cwas.utils.log import print_progress, print_arg, print_log
from cwas.runnable import Runnable
from scipy.stats import norm
from cwas.utils.check import check_is_file, check_is_dir
from scipy.stats import binomtest
import zarr

class EffectiveNumTest(Runnable):
    def __init__(self, args: argparse.Namespace):
        super().__init__(args)
        self._intersection_matrix = None
        self._correlation_matrix = None
        self._category_set_path = None
        self._category_set = None
        self._num_eig = None
        self.eff_num_test_value = None
        self._binom_p = None
        self._sample_info = None
        
    @staticmethod
    def _print_args(args: argparse.Namespace):
        print_arg("Input file", args.input_path)
        print_arg("Input format", args.input_format)
        print_arg("Output directory", args.output_dir_path)
        if args.sample_info_path:
            print_arg("Sample information file", args.sample_info_path)
        if args.category_count_file:
            print_arg("Using variant (or sample) counts file: ", args.category_count_file)
        if args.domain_list:
            print_arg("Domain list", args.domain_list)  
        if args.tag:
            print_arg("Output tag (prefix of output files)", args.tag)
        if args.category_set_path :
            print_arg("Category set file", args.category_set_path)

    @staticmethod
    def _check_args_validity(args: argparse.Namespace):
        check_is_dir(args.input_path)
        check_is_dir(args.output_dir_path)
        if args.sample_info_path:
            check_is_file(args.sample_info_path)
        if args.category_count_file:
            check_is_file(args.category_count_file)
        if args.category_set_path :
            check_is_file(args.category_set_path)
        if (args.sample_info_path is None) & (args.count_thres is None):
            raise Exception("Enter -s (--sample_info) or -thr (--threshold) to calculate threshold of the number of variants (or samples) of given categories")

    @property
    def input_path(self) -> Path:
        return self.args.input_path.resolve()

    @property
    def output_dir_path(self) -> Path:
        return(self.args.output_dir_path.resolve())

    @property
    def input_format(self) -> str:
        return self.args.input_format

    @property
    def category_count(self):
        category_count_ = pd.read_table(self.args.category_count_file, sep="\t")
        return category_count_

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
    def intersection_matrix(self) -> pd.DataFrame:
        if self._intersection_matrix is None and self.input_format == 'inter':
            root = zarr.open(self.input_path, mode='r')
            self._intersection_matrix = pd.DataFrame(data=root['data'],
                              index=root['metadata'].attrs['category'],
                              columns=root['metadata'].attrs['category'])
        return self._intersection_matrix

    @property
    def correlation_matrix(self) -> pd.DataFrame:
        if self._correlation_matrix is None and self.input_format == 'corr':
            root = zarr.open(self.input_path, mode='r')
            self._correlation_matrix = pd.DataFrame(data=root['data'],
                              index=root['metadata'].attrs['category'],
                              columns=root['metadata'].attrs['category'])
        return self._correlation_matrix

    @property
    def sample_info_path(self) -> Optional[Path]:
        return self.args.sample_info_path.resolve()

    @property
    def sample_info(self) -> pd.DataFrame:
        if self._sample_info is None:
            self._sample_info = pd.read_table(
                self.sample_info_path, index_col="SAMPLE", dtype={"SAMPLE": str}, sep="\t"
            )
        return self._sample_info

    @property
    def binom_p(self) -> float:
        if self._binom_p is None:
            self._binom_p = (self.sample_info["PHENOTYPE"] == "case").sum() / np.isin(self.sample_info["PHENOTYPE"], ["case", "ctrl"]).sum()
        return self._binom_p

    @property
    def count_thres(self) -> int:
        if self.args.count_thres is None:
            m = 1
            while True:
                p_value = binomtest(m-1, m, self.binom_p, alternative='greater').pvalue
                if p_value < 0.05:
                    return m
                m += 1
        else:
            return self.args.count_thres

    @property
    def num_eig(self) -> int:
        if self._num_eig is None:
            self._num_eig = self.args.num_eig
            return self._num_eig

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
    def neg_lap_path(self) -> Path:
        replace_term = r'\.intersection_matrix|\.correlation_matrix'
        f_name = re.sub(replace_term, '.neg_lap', self.input_path.name)
        save_name = '.zarr' if self.tag is None else f'.{self.tag}.zarr'
        f_name = re.sub('.zarr', save_name, f_name)
        return Path(
            f"{self.output_dir_path}/" +
            f"{f_name}"
        )

    @property
    def eig_val_path(self) -> Path:
        replace_term = r'\.intersection_matrix|\.correlation_matrix'
        f_name = re.sub(replace_term, '.eig_vals', self.input_path.name)
        save_name = '.zarr' if self.tag is None else f'.{self.tag}.zarr'
        f_name = re.sub('.zarr', save_name, f_name)
        return Path(
            f"{self.output_dir_path}/" +
            f"{f_name}"
        )

    @property
    def eig_vec_path(self) -> Path:
        replace_term = r'\.intersection_matrix|\.correlation_matrix'
        f_name = re.sub(replace_term, '.eig_vecs', self.input_path.name)
        save_name = '.zarr' if self.tag is None else f'.{self.tag}.zarr'
        f_name = re.sub('.zarr', save_name, f_name)
        return Path(
            f"{self.output_dir_path}/" +
            f"{f_name}"
        )

    def run(self):
        print_arg("Number of simulations", self.num_eig)
        for i in self.domain_list:
            print_progress(f"Eigen decomposition for domain: {i}")
            self._domain = i
            if self.eff_num_test:
                self.get_n_etests()
                self.update_env()
            else:
                self.eigen_decomposition()
        print_progress("Done")
    
    def get_n_etests(self):
        """Get the number of effective tests """
        print_progress(self.get_n_etests.__doc__)
        file_extension = '' if self._domain == 'all' else f'.{self._domain}'
        if os.path.isdir(Path(str(self.eig_val_path).replace('.zarr', f'.{self._domain}.zarr'))):
            print_log(
                "NOTICE",
                "You already have eigen values for the selected categories.",
                False,
            )
        else:
            self.eigen_decomposition(save_vecs=False)
            
        root = zarr.open(Path(str(self.eig_val_path).replace('.zarr', f'.{file_extension}.zarr')), mode='r')
        eig_vals = root['data']
        
        e = 1e-12
        eig_vals = sorted(eig_vals, key=np.linalg.norm, reverse=True)
        num_eig_val = self.num_eig
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
        
    def eigen_decomposition(self, save_vecs: bool = True):
        print_progress(f"Calculate eigen values")
        file_extension = '' if self._domain == 'all' else f'.{self._domain}'
        if self.category_set_path:
            c1 = self.category_set["Category"].tolist()
        else:
            if self.input_format == 'corr':
                c1 = self.correlation_matrix.columns.tolist()
            elif self.input_format == 'inter':
                c1 = self.intersection_matrix.columns.tolist()

        print_progress(f"Categories with at least {self.count_thres} counts will be used")
        c2 = self.category_count[self.category_count['Raw_counts'] >= self.count_thres]['Category'].tolist()
        
        filtered_combs1 = [x for x in c1 if x in c2]
        if self._domain != 'all':
            filtered_combs2 = self.category_set.loc[self.category_set['is_'+self._domain]==1]["Category"]
        else:
            filtered_combs2 = pd.Series(filtered_combs1)
        filtered_combs = list(set(filtered_combs1) & set(filtered_combs2))
        
        print_progress(f"Use # of categories: {len(filtered_combs)}")

        if self.input_format == 'corr':
            intermediate_mat = self.correlation_matrix.loc[filtered_combs,filtered_combs]
        elif self.input_format == 'inter':
            print_progress("Generating a covariance matrix")
            intermediate_mat = self.intersection_matrix.loc[filtered_combs,filtered_combs]
            intermediate_mat = intermediate_mat.mul((self.binom_p)*(1-self.binom_p))

        domain_neg_lap_path = Path(str(self.neg_lap_path).replace('.neg_lap', f'{file_extension}.neg_lap'))
        if not os.path.isdir(domain_neg_lap_path):
            print_progress("Generating the negative laplacian matrix")
            neg_lap = np.abs(intermediate_mat.values)
            degrees = np.sum(neg_lap, axis=0)
            for i in tqdm(range(neg_lap.shape[0])):
                neg_lap[i, :] = neg_lap[i, :] / np.sqrt(degrees)
                neg_lap[:, i] = neg_lap[:, i] / np.sqrt(degrees)
            print_progress("Writing the negative laplacian matrix to file")
            root = zarr.open(domain_neg_lap_path, mode='w')
            root.create_dataset('data', data=neg_lap, dtype='float64')
        else:
            root = zarr.open(domain_neg_lap_path, mode='r')
            neg_lap = root['data']

        domain_eig_val_path = Path(str(self.eig_val_path).replace('.eig_vals', f'{file_extension}.eig_vals'))
        domain_eig_vec_path = Path(str(self.eig_vec_path).replace('.eig_vecs', f'{file_extension}.eig_vecs'))
        if not os.path.isdir(domain_eig_val_path) or not os.path.isdir(domain_eig_vec_path):
            print_progress("Calculating the eigenvalues of the negative laplacian matrix")
            eig_vals, eig_vecs = np.linalg.eig(neg_lap)
            print_progress("Writing the eigenvalues to file")
            
            root = zarr.open(domain_eig_val_path, mode='w')
            root.create_dataset('data', data=eig_vals, dtype='float64')

            if save_vecs:
                print_progress("Writing the eigenvectors to file")
                root = zarr.open(domain_eig_vec_path, mode='w')
                root.create_group('metadata')
                root['metadata'].attrs['category'] = filtered_combs
                root.create_dataset('data', data=eig_vecs.real, chunks=(1000, 1000), dtype='float64')
                
    def update_env(self):
        self.set_env("N_EFFECTIVE_TEST", self.eff_num_test_value)
        self.save_env()
    
    