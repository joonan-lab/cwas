import argparse
import multiprocessing as mp
from pathlib import Path
import parmap
from cwas.core.common import chunk_list
from tqdm import tqdm
import re
import zarr
import yaml

import pandas as pd
import numpy as np
from functools import partial

import cwas.utils.log as log
from cwas.core.categorization.categorizer import Categorizer
from cwas.runnable import Runnable
from cwas.utils.check import check_num_proc, check_is_file, check_is_dir

from cwas.core.categorization.parser import (
    parse_annotated_vcf,
    parse_gene_matrix,
)

class Correlation(Runnable):
    def __init__(self, args: argparse.Namespace):
        super().__init__(args)
        self._annotated_vcf = None
        self._gene_matrix = None
        self._category_domain = None
        self._categorization_result = None
        self._categorization_root = None
        self._correlation_matrix = None
        self._intersection_matrix = None
        self._sample_ids = None
        self._categories = None
        self._adj_factor = None
        self._category_set_path = None
        self._category_set = None

    @staticmethod
    def _print_args(args: argparse.Namespace):
        if args.generate_corr_matrix == 'variant':
            log.print_arg("Annotated VCF file", args.annot_path)
        log.print_arg("Categorized file", args.cat_path)
        log.print_arg(
            "No. worker processes for the categorization",
            f"{args.num_proc: ,d}",
        )
        log.print_arg("Genereate an intersection matrix", (False if args.generate_inter_matrix is None else True))
        log.print_arg("Genereate a correlation matrix", (False if args.generate_corr_matrix is None else True))

    @staticmethod
    def _check_args_validity(args: argparse.Namespace):
        check_num_proc(args.num_proc)
        if args.generate_corr_matrix == 'variant':
            check_is_file(args.annot_path)
        check_is_dir(args.cat_path)
        check_is_dir(args.output_dir_path)

    @property
    def annot_path(self):
        return self.args.annot_path.resolve()

    @property
    def annotated_vcf(self) -> pd.DataFrame:
        if self._annotated_vcf is None:
            log.print_progress("Parse the annotated VCF")
            self._annotated_vcf = parse_annotated_vcf(
                Path(self.annot_path)
            )
        return self._annotated_vcf

    @property
    def gene_matrix(self) -> dict:
        if self._gene_matrix is None:
            self._gene_matrix = parse_gene_matrix(
                Path(self.get_env("GENE_MATRIX"))
            )
        return self._gene_matrix

    @property
    def category_domain(self) -> dict:
        if self._category_domain is None:
            with Path(self.get_env("CATEGORY_DOMAIN")).open(
                "r"
            ) as category_domain_file:
                self._category_domain = yaml.safe_load(category_domain_file)
        return self._category_domain

    @property
    def categorizer(self) -> Categorizer:
        categorizer = Categorizer(self.category_domain, self.gene_matrix, self.mis_info_key, self.mis_thres)
        return categorizer

    @property
    def mis_info_key(self) -> str:
        return self.get_env("VEP_MIS_INFO_KEY")

    @property
    def mis_thres(self) -> float:
        return float(self.get_env("VEP_MIS_THRES"))

    @property
    def num_proc(self):
        return self.args.num_proc

    @property
    def cat_path(self) -> Path:
        return self.args.cat_path.resolve()

    @property
    def category_info_path(self) -> Path:
        return self.args.category_info_path.resolve()

    @property
    def categorization_root(self):
        if self._categorization_root is None:
            self._categorization_root = zarr.open(self.cat_path, mode='r')
        return self._categorization_root

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
    def category_set(self) -> pd.DataFrame:
        if self._category_set is None:
            self._category_set = pd.read_csv(self.category_info_path, sep='\t')
            self._category_set = self._category_set.loc[self._category_set["Category"].isin(self.categories)]
            self._category_set['Category'] = pd.Categorical(self._category_set['Category'],
                                                            categories=self.categories,
                                                            ordered=True)
            self._category_set.sort_values('Category', ignore_index=True, inplace=True)
        return self._category_set

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
    def output_dir_path(self):
        return self.args.output_dir_path.resolve()

    @property
    def generate_corr_matrix(self):
        return self.args.generate_corr_matrix

    @property
    def generate_inter_matrix(self):
        return self.args.generate_inter_matrix

    @property
    def matrix_path(self) -> Path:
        f_name = re.sub(r'categorization_result\.zarr\.gz|categorization_result\.zarr', 'correlation_matrix.zarr', self.cat_path.name)
        return Path(
            f"{self.output_dir_path}/" +
            f"{f_name}"
        )

    @property
    def intersection_matrix_path(self) -> Path:
        f_name = re.sub(r'categorization_result\.zarr\.gz|categorization_result\.zarr', 'intersection_matrix.zarr', self.cat_path.name)
        return Path(
            f"{self.output_dir_path}/" +
            f"{f_name}"
        )

    def run(self):
        for i in self.domain_list:
            self._domain = i
            log.print_progress(f"Generate correlation matrix for domain: {i}")
            self.generate_correlation_matrix()
            self.save_result()
        log.print_progress("Done")

    def generate_correlation_matrix(self):
        self.filtered_combs = self.category_set.loc[self.category_set['is_'+self._domain]==1]['Category'] if self._domain != 'all' else pd.Series(self.categories)
        if self._domain != 'all':
            column_indices = [self.categories.index(col) for col in self.filtered_combs]
            self.categorization_result = pd.DataFrame(self.categorization_root['data'][:, column_indices].astype(np.float64),
                                                      index=self.sample_ids,
                                                      columns=self.filtered_combs)
            self.categorization_result.index.name = 'SAMPLE'
        else:
            self.categorization_result = pd.DataFrame(self.categorization_root['data'].astype(np.float64),
                                                      index=self.sample_ids,
                                                      columns=self.categories)
            self.categorization_result.index.name = 'SAMPLE'
        if self.generate_corr_matrix == "sample":
            log.print_progress("Get an intersection matrix between categories using the number of samples")

            if self.num_proc == 1:
                intersection_matrix = self.process_columns_single(column_range = range(self.categorization_result.shape[1]), matrix=self.categorization_result)
            else:
                # Split the column range into evenly sized chunks based on the number of workers
                chunks = chunk_list(range(self.categorization_result.shape[1]), self.num_proc)
                result = parmap.map(self.process_columns, chunks, matrix=self.categorization_result, pm_pbar=True, pm_processes=self.num_proc)
                # Concatenate the count values
                intersection_matrix = pd.concat([pd.concat(chunk_results, axis=1) for chunk_results in result], axis=1)

        elif self.generate_corr_matrix == "variant":
            log.print_progress("Get an intersection matrix between categories using the number of variants")
            #pre_intersection_matrix = self.categorizer.get_intersection_variant_level(self.annotated_vcf, self.categorization_result.columns.tolist())
            intersection_matrix = (
                self.get_intersection_matrix(self.annotated_vcf, self.categorizer, self.categorization_result.columns)
                if self.num_proc == 1
                else self.get_intersection_matrix_with_mp()
            )
            #if self.num_proc == 1:
            #    intersection_matrix = self.process_columns_single(column_range = range(pre_intersection_matrix.shape[1]), matrix=pre_intersection_matrix)
            #else:
            #    # Split the column range into evenly sized chunks based on the number of workers
            #    log.print_progress(f"This step will use only {self.num_proc//3 + 1} worker processes to avoid memory error")
            #    chunks = chunk_list(range(pre_intersection_matrix.shape[1]), (self.num_proc//3 + 1))
            #    result = parmap.map(self.process_columns, chunks, matrix=pre_intersection_matrix, pm_pbar=True, pm_processes=(self.num_proc//3 + 1))
            #    # Concatenate the count values
            #    intersection_matrix = pd.concat([pd.concat(chunk_results, axis=1) for chunk_results in result], axis=1)
        
        diag_sqrt = np.sqrt(np.diag(intersection_matrix))
        log.print_progress("Calculate a correlation matrix")
        self._intersection_matrix = intersection_matrix
        self._correlation_matrix = intersection_matrix/np.outer(diag_sqrt, diag_sqrt)

    @staticmethod
    def process_columns(column_range, matrix: pd.DataFrame) -> list:
        results = []
    
        # Iterate over the column range
        for i in column_range:
            # Multiply the i-th column with values in the matrix
            df_multiplied = matrix.mul(matrix.iloc[:, i], axis=0)
            
            # Count the number of values greater than 0 in each column
            count_values_gt_zero = (df_multiplied > 0).sum(axis=0)
            
            # Assign the column name to count_values_gt_zero
            count_values_gt_zero.name = matrix.columns[i]
            
            results.append(count_values_gt_zero)
        
        return results

    @staticmethod
    def process_columns_single(column_range, matrix: pd.DataFrame) -> pd.DataFrame:
        # Initialize an empty DataFrame to store the concatenated results
        result = pd.DataFrame()

        # Define the progress bar
        pbar = tqdm(column_range, desc='Processing')

        # Perform the multiplication in a loop
        for i in pbar:
            # Multiply the i-th column with values in the matrix
            df_multiplied = matrix.mul(matrix.iloc[:, i], axis=0)
            
            # Count the number of values greater than 0 in each column
            count_values_gt_zero = (df_multiplied > 0).sum(axis=0)
            
            # Assign the column name to count_values_gt_zero
            count_values_gt_zero.name = matrix.columns[i]
            
            # Concatenate the count values to the 'result' DataFrame
            result = pd.concat([result, count_values_gt_zero], axis=1)

        # Close the progress bar
        pbar.close()
        
        return result

    def get_intersection_matrix_with_mp(self):
        ## use only one third of the cores to avoid memory error
        #log.print_progress(f"This step will use only {self.num_proc//3 + 1} worker processes to avoid memory error")
        split_vcfs = np.array_split(self.annotated_vcf, self.num_proc)
        _get_intersection_matrix = partial(self.get_intersection_matrix,
                                           categorizer=self.categorizer, 
                                           categories=self.categorization_result.columns)
        
        with mp.Pool(self.num_proc) as pool:
            return sum(pool.map(
                _get_intersection_matrix,
                split_vcfs
            ))
        
    @staticmethod
    def get_intersection_matrix(annotated_vcf: pd.DataFrame, categorizer: Categorizer, categories: pd.Index): 
        return pd.DataFrame(
            categorizer.get_intersection(annotated_vcf), 
            index=categories, 
            columns=categories
        ).fillna(0).astype(int)

    def save_result(self):
        if self.generate_inter_matrix == True:
            log.print_progress("Save the intersection matrix to file")

            domain_int_path = Path(str(self.intersection_matrix_path).replace('.zarr', f'.{self._domain}.zarr'))
            root = zarr.open(domain_int_path, mode='w')
            root.create_group('metadata')
            root['metadata'].attrs['category'] = self._intersection_matrix.columns.tolist()
            root.create_dataset('data', data=self._intersection_matrix, chunks=(1000, 1000), dtype='i4')

        log.print_progress("Save the correlation matrix to file")
        domain_corr_path = Path(str(self.matrix_path).replace('.zarr', f'.{self._domain}.zarr'))
        root = zarr.open(domain_corr_path, mode='w')
        root.create_group('metadata')
        root['metadata'].attrs['category'] = self._correlation_matrix.columns.tolist()
        root.create_dataset('data', data=self._correlation_matrix, chunks=(1000, 1000), dtype='float64')

