import argparse, os, sys, gzip, glob
import multiprocessing as mp
from pathlib import Path
import numpy as np
import pandas as pd
from typing import Optional

from cwas.utils.log import print_progress, print_arg, print_warn, print_log
from cwas.runnable import Runnable
from functools import partial
from cwas.utils.check import check_is_file, check_is_dir, check_num_proc
from scipy.stats import binom_test, norm

class ConcatZscore(Runnable):
    def __init__(self, args: argparse.Namespace):
        super().__init__(args)
        self._in_dir = None
        self._out_dir = None
        self._burden_test_paths = None
        self._sample_info = None
        self._category_set_path = None
        self._category_set = None
        self._tag = None
        self._binom_p = None
        self._zscore_df = None

    @staticmethod
    def _print_args(args: argparse.Namespace):
        print_arg("Input directory", args.input_dir_path)
        print_arg("Sample information file", args.sample_info_path)
        print_arg("Output directory", args.output_dir_path)
        if args.tag:
            print_arg("Output tag (prefix of output files)", args.tag)
        print_arg("Category set file", args.category_set_path)
        print_arg(
            "No. worker processes for concatenation",
            f"{args.num_proc: ,d}",
        )
    
    @staticmethod
    def _check_args_validity(args: argparse.Namespace):
        check_is_dir(args.input_dir_path)
        check_is_dir(args.output_dir_path)
        check_is_file(args.sample_info_path)
        if args.category_set_path :
            check_is_file(args.category_set_path)
        check_num_proc(args.num_proc)

    @property
    def in_dir(self) -> Path:
        return (
            self.args.in_dir_path.resolve()
        )

    @property
    def out_dir(self) -> Path:
        return(self.args.output_dir_path.resolve())

    @property
    def num_proc(self) -> int:
        return self.args.num_proc            

    @property
    def sample_info_path(self) -> Path:
        return self.args.sample_info_path.resolve()

    @property
    def tag(self) -> str:
        return self.args.tag
    
    @property
    def sample_info(self) -> pd.DataFrame:
        if self._sample_info is None:
            self._sample_info = pd.read_table(
                self.sample_info_path, index_col="SAMPLE"
            )
        return self._sample_info

    @property
    def burden_test_paths(self) -> list:
        if self._burden_test_paths is None:
            self._burden_test_paths = sorted(glob(f'{self.in_dir}/*.burden_test.txt.gz'))
        return self._burden_test_paths

    @property
    def zscore_df_path(self) -> Path:
        if self.tag is None:
            save_name = 'output.zscores.txt.gz'
        else:
            save_name = '.'.join(['output', self.tag, 'zscores.txt.gz'])
        return Path(
            f"{self.out_dir}/" +
            save_name
        )

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
    def binom_p(self) -> float:
        if self._binom_p is None:
            self._binom_p = (self.sample_info["PHENOTYPE"] == "case").sum() / np.isin(self.sample_info["PHENOTYPE"], ["case", "ctrl"]).sum()
        return self._binom_p
    
    @staticmethod
    def _get_zscore_dict(burden_test_path: Path, z_dict: dict) -> dict:
        with gzip.open(burden_test_path, mode='rt') as f:
            _ = f.readline()
            for line in f:
                fields = line.strip().split('\t')
            ## Get z score for each category
                if z_dict.get(fields[0]) is not None:
                    z_dict[fields[0]] = float(fields[11])
                    
        return z_dict
    
    def concat_zscores(self):
        """Concatenate zscores for each result """
        print_progress(self.concat_zscores.__doc__)
        
        if self.zscore_df_path.is_file():
            print_log(
                "NOTICE",
                "You already have a Z score table. Skip this step.",
                False,
            )
            return

        default_p = binom_test(x=1, n=2, p=self.binom_p, alternative='greater')
        default_z = norm.ppf(1 - default_p)
        default_z_dict = {comb: default_z for comb in self.category_set['Category']}
        
        _get_zscore_dict_partial = partial(self._get_zscore_dict, z_dict=default_z_dict)
        
        if self.num_proc == 1:
            z_dicts = []
            for burden_test_path in self.burden_test_paths:
                z_dict = _get_zscore_dict_partial(burden_test_path)
                z_dicts.append(z_dict)
        else:
            with mp.Pool(self.num_proc) as pool:
                z_dicts = pool.map(
                    _get_zscore_dict_partial,
                    self.burden_test_paths,
                )
                
        zscore_df = pd.DataFrame.from_records(z_dicts)
        zscore_df.index.name = 'Simulation'
        zscore_df.index += 1
        
        print_progress("Save a Z score table")
        zscore_df.to_csv(self.zscore_df_path, sep='\t')
        self._zscore_df = zscore_df
    
    def run(self):
        self.concat_zscores()
        print_progress("Done")