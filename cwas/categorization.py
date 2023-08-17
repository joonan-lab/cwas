import argparse
import multiprocessing as mp
import multiprocessing.shared_memory as shm
from functools import reduce
from itertools import product
from math import ceil
from pathlib import Path
import parmap
from cwas.core.common import chunk_list
from tqdm import tqdm
import re

import pandas as pd
import numpy as np
import yaml, pickle
from functools import partial

import cwas.utils.log as log
from cwas.core.categorization.categorizer import Categorizer
from cwas.core.categorization.parser import (
    parse_annotated_vcf,
    parse_gene_matrix,
)
from cwas.runnable import Runnable
from cwas.utils.check import check_num_proc
from cwas.utils.check import check_is_file
from cwas.utils.check import check_is_dir

lock = mp.Lock()

class Categorization(Runnable):
    def __init__(self, args: argparse.Namespace):
        super().__init__(args)
        self._annotated_vcf = None
        self._gene_matrix = None
        self._category_domain = None
        self._annotated_vcf_groupby_sample = None
        self._sample_ids = None
        self._result = None
        self._correlation_matrix = None
        self._intersection_matrix = None
        self._shm = None

    @staticmethod
    def _print_args(args: argparse.Namespace):
        log.print_arg(
            "No. worker processes for the categorization",
            f"{args.num_proc: ,d}",
        )
        log.print_arg("Annotated VCF file", args.input_path)
        log.print_arg("Genereate a correlation matrix and an intersection matrix", (False if args.generate_matrix is None else True))

    @staticmethod
    def _check_args_validity(args: argparse.Namespace):
        check_num_proc(args.num_proc)
        check_is_file(args.input_path)
        check_is_dir(args.output_dir_path)

    @property
    def num_proc(self):
        return self.args.num_proc

    @property
    def input_path(self):
        return self.args.input_path.resolve()

    @property
    def output_dir_path(self):
        return self.args.output_dir_path.resolve()

    @property
    def generate_matrix(self):
        return self.args.generate_matrix

    @property
    def mis_info_key(self) -> str:
        return self.get_env("VEP_MIS_INFO_KEY")

    @property
    def mis_thres(self) -> float:
        return float(self.get_env("VEP_MIS_THRES"))

    @property
    def result_path(self) -> Path:
        f_name = re.sub(r'annotated\.vcf\.gz|annotated\.vcf', 'categorization_result.txt.gz', self.input_path.name)
        return Path(
            f"{self.output_dir_path}/" + 
            f"{f_name}"
        )

    @property
    def matrix_path(self) -> Path:
        f_name = re.sub(r'annotated\.vcf\.gz|annotated\.vcf', 'correlation_matrix.pkl', self.input_path.name)
        return Path(
            f"{self.output_dir_path}/" +
            f"{f_name}"
        )

    @property
    def intersection_matrix_path(self) -> Path:
        f_name = re.sub(r'annotated\.vcf\.gz|annotated\.vcf', 'intersection_matrix.pkl', self.input_path.name)
        return Path(
            f"{self.output_dir_path}/" +
            f"{f_name}"
        )

    @property
    def annotated_vcf(self) -> pd.DataFrame:
        if self._annotated_vcf is None:
            log.print_progress("Parse the annotated VCF")
            self._annotated_vcf = parse_annotated_vcf(
                Path(self.input_path)
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
    def annotated_vcf_groupby_sample(self):
        if self._annotated_vcf_groupby_sample is None:
            self._annotated_vcf_groupby_sample = self.annotated_vcf.groupby(
                "SAMPLE"
            )
        return self._annotated_vcf_groupby_sample

    @property 
    def sample_ids(self) -> list:
        if self._sample_ids is None:
            self._sample_ids = list(self.annotated_vcf_groupby_sample.groups)
        return self._sample_ids

    @property
    def annotated_vcf_split_by_sample(self) -> list:
        return [
            self.annotated_vcf_groupby_sample.get_group(sample_id)
            for sample_id in self.sample_ids
        ]

    @property
    def redundant_categories(self) -> set:
        redundant_category_table = pd.read_table(
            self.get_env("REDUNDANT_CATEGORY")
        )
        return reduce(
            lambda x, y: x.union(y),
            [
                self._get_redundant_categories_from_row(row)
                for _, row in redundant_category_table.iterrows()
            ],
            set(),
        )

    def _get_redundant_categories_from_row(self, row: pd.Series) -> list:
        annotation_term_lists = {}
        for k, v in row.items():
            if v == "*":
                annotation_term_lists[k] = self._category_domain[k]
            else:
                annotation_term_lists[k] = [v]

        return {
            "_".join(combination)
            for combination in product(
                annotation_term_lists["variant_type"],
                annotation_term_lists["gene_set"],
                annotation_term_lists["functional_score"],
                annotation_term_lists["gencode"],
                annotation_term_lists["functional_annotation"],
            )
        }

    def run(self):
        self.categorize_vcf()
        self.remove_redundant_category()
        self.generate_correlation_matrix()
        self.save_result()
        self.update_env()
        log.print_progress("Done")

    def categorize_vcf(self):
        results_each_sample = (
            self.categorize_vcf_for_each_sample()
            if self.num_proc == 1
            else self.categorize_vcf_for_each_sample_with_mp()
        )
        log.print_progress("Organize the results")
        self._result = pd.DataFrame(results_each_sample).fillna(0)
        self._result = self._result.astype(int)
        self._result["SAMPLE"] = self.sample_ids
        self._result.set_index("SAMPLE", inplace=True)

    def remove_redundant_category(self):
        log.print_progress("Remove redundant categories from the result")
        self._result.drop(
            self.redundant_categories,
            axis="columns",
            inplace=True,
            errors="ignore",
        )
        log.print_progress(
            f"{len(self._result.columns):,d} categories have remained."
        )

    def categorize_vcf_for_each_sample(self):
        result = []
        for i, sample_vcf in enumerate(self.annotated_vcf_split_by_sample, 1):
            if i % 100 == 0:
                log.print_progress(f"Categorize variants of {i} samples")
            result.append(self.categorizer.categorize_variant(sample_vcf))
        return result

    def categorize_vcf_for_each_sample_with_mp(self):
        sample_vcfs = self.annotated_vcf_split_by_sample
        with mp.Pool(self.num_proc) as pool:
            log.print_progress("Categorize your input variants")
            return pool.map(
                self.categorizer.categorize_variant,
                sample_vcfs,
                chunksize=ceil(len(sample_vcfs) / self.num_proc),
            )

    def generate_correlation_matrix(self):
        if self.generate_matrix is None:
            return

        if self.generate_matrix == "sample":
            log.print_progress("Get an intersection matrix between categories using the number of samples")

            if self.num_proc == 1:
                intersection_matrix = self.process_columns_single(column_range = range(self._result.shape[1]), matrix=self._result)
            else:
                # Split the column range into evenly sized chunks based on the number of workers
                chunks = chunk_list(range(self._result.shape[1]), self.num_proc)
                result = parmap.map(self.process_columns, chunks, matrix=self._result, pm_pbar=True, pm_processes=self.num_proc)
                # Concatenate the count values
                intersection_matrix = pd.concat([pd.concat(chunk_results, axis=1) for chunk_results in result], axis=1)

        elif self.generate_matrix == "variant":
            log.print_progress("Get an intersection matrix between categories using the number of variants")
            intersection_matrix = (
                self.get_intersection_matrix(self.annotated_vcf, self.categorizer, self._result.columns)
                if self.num_proc == 1
                else self.get_intersection_matrix_with_mp()
            )
        
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
        log.print_progress(f"This step will use only {self.num_proc//3 + 1} worker processes to avoid memory error")
        split_vcfs = np.array_split(self.annotated_vcf, self.num_proc//3 + 1)
        _get_intersection_matrix = partial(self.get_intersection_matrix,
                                           categorizer=self.categorizer, 
                                           categories=self._result.columns)
        
        with mp.Pool(self.num_proc//3 + 1) as pool:
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
        log.print_progress(f"Save the result to the file {self.result_path}")
        self._result.to_csv(self.result_path, sep="\t")
        if self._correlation_matrix is not None:
            log.print_progress("Save the intersection matrix to file")
            pickle.dump(self._intersection_matrix, open(self.intersection_matrix_path, 'wb'), protocol=5)
            log.print_progress("Save the correlation matrix to file")
            pickle.dump(self._correlation_matrix, open(self.matrix_path, 'wb'), protocol=5)

    def update_env(self):
        self.set_env("CATEGORIZATION_RESULT", self.result_path)
        self.save_env()
