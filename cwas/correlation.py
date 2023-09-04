import argparse
import multiprocessing as mp
from pathlib import Path
import parmap
from cwas.core.common import chunk_list
from tqdm import tqdm
import re
import polars as pl

import pandas as pd
import numpy as np
import pickle
from functools import partial

import cwas.utils.log as log
from cwas.core.categorization.categorizer import Categorizer
from cwas.runnable import Runnable
from cwas.utils.check import check_num_proc
from cwas.utils.check import check_is_file
from cwas.utils.check import check_is_dir

class Correlation(Runnable):
    def __init__(self, args: argparse.Namespace):
        super().__init__(args)
        self._categorization_result = None
        self._correlation_matrix = None
        self._intersection_matrix = None

    @staticmethod
    def _print_args(args: argparse.Namespace):
        log.print_arg(
            "No. worker processes for the categorization",
            f"{args.num_proc: ,d}",
        )
        log.print_arg("Categorized file", args.cat_path)
        log.print_arg("Genereate a correlation matrix and an intersection matrix", (False if args.generate_matrix is None else True))

    @staticmethod
    def _check_args_validity(args: argparse.Namespace):
        check_num_proc(args.num_proc)
        check_is_file(args.cat_path)
        check_is_dir(args.output_dir_path)

    @property
    def num_proc(self):
        return self.args.num_proc

    @property
    def cat_path(self) -> Path:
        return self.args.cat_path.resolve()

    @property
    def categorization_result(self) -> pd.DataFrame:
        if self._categorization_result is None:
            log.print_progress("Load the categorization result")
            self._categorization_result = pl.read_csv(
                self.cat_path, separator="\t", dtypes={"SAMPLE": str}
            )
            self._categorization_result = self._categorization_result.to_pandas().set_index("SAMPLE")
        return self._categorization_result

    @property
    def output_dir_path(self):
        return self.args.output_dir_path.resolve()

    @property
    def generate_matrix(self):
        return self.args.generate_matrix

    @property
    def matrix_path(self) -> Path:
        f_name = re.sub(r'annotated\.vcf\.gz|annotated\.vcf', 'correlation_matrix.pkl', self.cat_path.name)
        return Path(
            f"{self.output_dir_path}/" +
            f"{f_name}"
        )

    @property
    def intersection_matrix_path(self) -> Path:
        f_name = re.sub(r'annotated\.vcf\.gz|annotated\.vcf', 'intersection_matrix.pkl', self.cat_path.name)
        return Path(
            f"{self.output_dir_path}/" +
            f"{f_name}"
        )

    def run(self):
        self.generate_correlation_matrix()
        self.save_result()
        self.update_env()
        log.print_progress("Done")

    def generate_correlation_matrix(self):
        if self.generate_matrix is None:
            return

        if self.generate_matrix == "sample":
            log.print_progress("Get an intersection matrix between categories using the number of samples")

            if self.num_proc == 1:
                intersection_matrix = self.process_columns_single(column_range = range(self.categorization_result.shape[1]), matrix=self.categorization_result)
            else:
                start = 0
                end = self.categorization_result.shape[1]
                step = 300
                
                chunks = [range(i, min(i + step, end)) for i in range(start, end, step)]
                
                # Split the column range into evenly sized chunks based on the number of workers
                #chunks = chunk_list(range(self.categorization_result.shape[1]), self.num_proc)
                result = parmap.map(self.process_columns, chunks, matrix=self.categorization_result, pm_pbar=True, pm_processes=self.num_proc)
                # Concatenate the count values
                intersection_matrix = pd.concat([pd.concat(chunk_results, axis=1) for chunk_results in result], axis=1)

        elif self.generate_matrix == "variant":
            log.print_progress("Not prepared")
            return
        
        diag_sqrt = np.sqrt(np.diag(intersection_matrix))
        log.print_progress("Calculate a correlation matrix")
        self._intersection_matrix = intersection_matrix
        self._correlation_matrix = intersection_matrix/np.outer(diag_sqrt, diag_sqrt)

    @staticmethod
    def process_columns(column_range, matrix: pd.DataFrame) -> list:
        results = []

        pbar = tqdm(column_range, desc='Processing')
    
        # Iterate over the column range
        for i in pbar:
            # Multiply the i-th column with values in the matrix
            df_multiplied = matrix.mul(matrix.iloc[:, i], axis=0)
            
            # Count the number of values greater than 0 in each column
            count_values_gt_zero = (df_multiplied > 0).sum(axis=0)
            
            # Assign the column name to count_values_gt_zero
            count_values_gt_zero.name = matrix.columns[i]
            
            results.append(count_values_gt_zero)
        
        pbar.close()

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
        
        _get_intersection_matrix = partial(get_intersection_matrix_,
                                           categorizer=self.categorizer)

        #pool = mp.Pool(processes=self.num_proc//3 + 1)

        with mp.Pool(processes=self.num_proc//3 + 1) as pool:
            split_results = pool.map(_get_intersection_matrix, split_vcfs)
        
        #split_results = pool.map(_get_intersection_matrix, split_vcfs)

        # Initialize a variable to store the final summed DataFrame
        summed_df = None

        log.print_progress(f"Gather multiprocessed outputs")

        # Loop through the list of dictionaries
        for result_dict in split_results:
            # Convert the dictionary into a DataFrame
            df = pd.DataFrame(result_dict,
                              index=self.categorization_result.columns,
                              columns=self.categorization_result.columns).fillna(0).astype(int)
            
            # Sum the DataFrame with the existing summed DataFrame
            if summed_df is None:
                summed_df = df
            else:
                summed_df = summed_df.add(df, fill_value=0)
        
        return summed_df
        #with mp.Pool(self.num_proc//3 + 1) as pool:
        #    return sum(pool.map(
        #        _get_intersection_matrix,
        #        split_vcfs
        #    ))
        
    @staticmethod
    def get_intersection_matrix(annotated_vcf: pd.DataFrame, categorizer: Categorizer, categories: pd.Index): 
        return pd.DataFrame(
            categorizer.get_intersection(annotated_vcf), 
            index=categories, 
            columns=categories
        ).fillna(0).astype(int)

    def save_result(self):
        if self._correlation_matrix is not None:
            log.print_progress("Save the intersection matrix to file")
            pickle.dump(self._intersection_matrix, open(self.intersection_matrix_path, 'wb'), protocol=5)
            log.print_progress("Save the correlation matrix to file")
            pickle.dump(self._correlation_matrix, open(self.matrix_path, 'wb'), protocol=5)

    def update_env(self):
        self.save_env()

def get_intersection_matrix_(annotated_vcf: pd.DataFrame, categorizer: Categorizer): 
    return categorizer.get_intersection(annotated_vcf)