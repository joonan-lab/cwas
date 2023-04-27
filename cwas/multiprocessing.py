import yaml, os, gzip, sys, argparse
import multiprocessing as mp
from functools import partial
import re
import numpy as np
import pandas as pd
from pathlib import Path
from cwas.annotation import Annotation
from cwas.categorization import Categorization
from cwas.binomial_test import BinomialTest
from typing import Optional
from cwas.runnable import Runnable
from cwas.core.categorization.parser import parse_annotated_vcf
import cwas.utils.log as log
from cwas.utils.check import check_is_file, check_num_proc, check_same_n_lines, check_is_dir
from tqdm import tqdm

class Multiprocessing(Runnable):
    def __init__(self, args: argparse.Namespace):
        super().__init__(args)
        self._in_dir = None
        self._sample_info = None
        self._out_dir = None
        self._rand_mut_paths = None
        self._resume = None
        self._num_sim = None
    
    @staticmethod
    def _print_args(args: argparse.Namespace):
        log.print_arg(f"CWAS step:", str(args.step))
        log.print_arg("Input directory", args.in_vcf_path if args.in_vcf_path else "Not specified: $ANNOTATED_VCF will be used")
        if args.step == 'binomial_test':
            log.print_arg("Sample information file", args.sample_info_path)
            if args.adj_file_path :
                log.print_arg("Adjust factor file", args.adj_file_path)
        log.print_arg("Output directory", args.out_dir)
        log.print_arg("Number of simulations", args.num_sim)
        log.print_arg(
            "No. worker processes for this step",
            f"{args.num_mp: ,d}",
        )
        log.print_arg("No. processes for each process in this step", f"{args.num_proc: ,d}")
    
    @staticmethod
    def _check_args_validity(args: argparse.Namespace):
        check_is_dir(args.in_dir_path)
        check_is_file(args.sample_info_path)
        check_num_proc(args.num_mp) # multitply?
        check_num_proc(args.num_proc)
        if args.step == 'binomial_test':
            check_is_file(args.sample_info_path)
            if args.adj_file_path :
                check_is_file(args.adj_file_path)

    @property
    def in_dir(self) -> Path:
        return (
            self.args.in_dir_path.resolve()
        )
    
    @property
    def out_dir(self) -> Path:
        if self.args.output_dir_path:
            return(self.args.output_dir_path.resolve())
        else:
            if self.args.step == 'annotation':
                dir_name = 'vep'
            if self.args.step == 'categorization':
                dir_name = 'cat'
            if self.args.step == 'binomial_test':
                dir_name = 'burden'
            return(
                Path('/'.join([self.get_env("CWAS_WORKSPACE"), '_'.join(["random-mutations", dir_name])]))
            )

    @property
    def num_mp(self) -> int:
        return self.args.num_mp

    @property
    def num_proc(self) -> int:
        return self.args.num_proc

    @property
    def step(self) -> str:
        return self.args.step

    @property
    def resume(self) -> bool:
        if self._resume is None:
            self._resume = self.args.resume
        return self._resume

    @property
    def rand_mut_paths(self) -> list:
        if self._rand_mut_paths is None:
            self._rand_mut_paths = sorted(self.in_dir.glob(f'{self.in_dir}/*.vcf.gz'))
        return self._rand_mut_paths

    @property
    def annot_vcf_paths(self) -> list:
        if self._annot_vcf_paths is None:
            self._annot_vcf_paths = sorted(self.in_dir.glob(f'{self.in_dir}/*.annotated.vcf.gz'))
        return self._annot_vcf_paths

    @property
    def cat_result_paths(self) -> list:
        if self._cat_result_paths is None:
            self._cat_result_paths = sorted(self.in_dir.glob(f'{self.in_dir}/*.categorization_result.txt.gz'))
        return self._cat_result_paths
    
    @property
    def burden_test_paths(self) -> list:
        if self._burden_test_paths is None:
            self._burden_test_paths = sorted(self.in_dir.glob(f'{self.in_dir}/*.burden_test.txt.gz'))
        return self._burden_test_paths

    @property
    def num_sim(self) -> int:
        if self._num_sim is None:
            if self.step == 'annotation':
                self._num_sim = len(sorted(self.in_dir.glob(f'{self.in_dir}/*.vcf.gz')))
            if self.step == 'categorization':
                self._num_sim = len(self.annot_vcf_paths)
            if self.step == 'binomial_test':
                self._num_sim = len(self.cat_result_paths)
        return self._num_sim

    def multi_annotate(self):
        """Annotation for random mutations """
        log.print_progress(self.multi_annotate.__doc__)
        pre_files = sorted(self.out_dir.glob(f'{self.out_dir}/*.annotated.vcf.gz'))
        if len(pre_files) == 0:
            target_inputs = self.rand_mut_paths
        elif len(pre_files) == self.num_sim:
            log.print_progress("Checking the number of lines for the last 100 VCFs...")
            check_same_n_lines(pre_files[-100:], gzip_file=True)
            log.print_log(
                "NOTICE",
                "You already have annotated vcfs. Skip this step.",
                False,
            )
            return
        elif self.resume & (len(pre_files) < self.num_sim):
            log.print_progress("Checking the number of lines for the last 100 VCFs...")
            check_same_n_lines(pre_files[-100:], gzip_file=True)
            log.print_log(
                "NOTICE",
                f"You have some annotated vcfs ({len(pre_files)}). Resume this step.",
                False,
            )
            target_inputs = sorted(list(set(self.rand_mut_paths) - set([Path(str(path).replace('.annotated.vcf.gz', '.vcf.gz')) for path in pre_files])))
            self._resume = False
        else:
            raise RuntimeError(
                "The number of annotated vcfs is not the same as the number of simulations."
                " Check and remove the files to rerun this step."
            )

        _annotate_one_partial = partial(self._annotate_one, 
                                       num_proc=self.num_proc,
                                       out_dir = self.out_dir)
        
        if self.num_proc == 1:
            for rand_mut_path in tqdm(target_inputs):
                _annotate_one_partial(
                    rand_mut_path,
                )
        else:
            def mute():
                sys.stderr = open(os.devnull, 'w')
            with mp.Pool(self.num_proc, initializer=mute) as pool:
                for _ in tqdm(pool.imap_unordered(_annotate_one_partial, target_inputs), total=len(target_inputs)):
                    pass

    @staticmethod
    def _annotate_one(rand_mut_path: Path, num_proc: int, out_dir: Path):
        annotator = Annotation.get_instance(['-v', str(rand_mut_path), '-p', str(num_proc), '-o_dir', str(out_dir)])
        annotator.vep_output_vcf_path = str(rand_mut_path).replace('.vcf.gz', '.vep.vcf')
        annotator.annotate_using_bigwig()
        annotator.process_vep_vcf()
        annotator.annotate_using_bed()

    def multi_categorize(self):
        """Categorize random mutations """
        log.print_progress(self.multi_categorize.__doc__)
        pre_files = sorted(self.out_dir.glob(f'{self.out_dir}/*.categorization_result.txt.gz'))
        if len(pre_files) == 0:
            target_inputs = self.annot_vcf_paths
        elif len(pre_files) == self.num_sim:
            log.print_progress("Checking the number of lines for the last 100 result files...")
            check_same_n_lines(pre_files[-100:], gzip_file=True)
            log.print_log(
                "NOTICE",
                "You already have categorization results. Skip this step.",
                False,
            )
            return
        elif self.resume & (len(pre_files) < self.num_sim):
            log.print_progress("Checking the number of lines for the last 100 result files...")
            check_same_n_lines(pre_files[-100:], gzip_file=True)
            log.print_log(
                "NOTICE",
                f"You have some categorization results ({len(pre_files)}). Resume this step.",
                False,
            )
            target_inputs = sorted(list(set(self.annot_vcf_paths) - set([Path(str(path).replace('.categorization_result.txt.gz', '.annotated.vcf.gz')) for path in pre_files])))
            self._resume = False
        else:
            raise RuntimeError(
                "The number of categorization results is not the same as the number of simulations."
                " Check and remove the files to rerun this step."
            )

        _categorize_one_partial = partial(self._categorize_one, 
                                       num_proc=self.num_proc,
                                       out_dir = self.out_dir)
        
        if self.num_proc == 1:
            for annot_vcf_path in tqdm(target_inputs):
                _categorize_one_partial(
                    annot_vcf_path,
                )
        else:
            def mute():
                sys.stderr = open(os.devnull, 'w')
            with mp.Pool(self.num_proc, initializer=mute) as pool:
                for _ in tqdm(pool.imap_unordered(_categorize_one_partial, target_inputs), total=len(target_inputs)):
                    pass

    @staticmethod
    def _categorize_one(annot_vcf_path: Path, num_proc: int, out_dir: Path):
        categorizer = Categorization.get_instance(['-i', str(annot_vcf_path), '-p', str(num_proc), '-o_dir', str(out_dir)])
        categorizer.annotated_vcf = parse_annotated_vcf(annot_vcf_path)
        categorizer.result_path = Path(str(annot_vcf_path).replace('.annotated.vcf.gz', '.categorization_result.txt.gz'))
        categorizer.categorize_vcf()
        categorizer.remove_redundant_category()
        categorizer.save_result()

    def multi_binomial_test(self):
        """Burden tests for random mutations """
        log.print_progress(self.multi_binomial_test.__doc__)
        pre_files = self.burden_test_paths
        if len(pre_files) == 0:
            target_inputs = self.cat_result_paths
        elif len(pre_files) == self.num_sim:
            log.print_log(
                "NOTICE",
                "You already have burden test results.",
                False,
            )
            return
        elif self.resume & (len(pre_files) < self.num_sim):
            log.print_log(
                "WARNING",
                f"You have some burden test results ({len(pre_files)}). "
                f"Resume this step. Check whether the files are truncated.",
                False,
            )
            target_inputs = sorted(list(set(self.cat_result_paths) - set([Path(str(path).replace('.burden_test.txt.gz', '.categorization_result.txt.gz')) for path in pre_files])))
            self._resume = False
        else:
            raise RuntimeError(
                "The number of burden test results is not the same as the number of simulations."
                " Check and remove the files to rerun this step."
            )

        _burden_test_partial = partial(self._burden_test, 
                                       sample_info_path=self.sample_info_path,
                                       use_n_carrier=self.use_n_carrier,
                                       adj_factor_path=self.adj_factor_path,
                                       out_dir = self.out_dir)
        
        if self.num_proc == 1:
            for cat_result_path in target_inputs:
                _burden_test_partial(
                    cat_result_path,
                )
        else:
            def mute():
                sys.stderr = open(os.devnull, 'w')
            with mp.Pool(self.num_proc, initializer=mute) as pool:
                for _ in tqdm(pool.imap_unordered(_burden_test_partial, target_inputs), total=len(target_inputs)):
                    pass


    @staticmethod
    def _burden_test(cat_result_path: Path, sample_info_path: Path, out_dir: Path, adj_factor_path: Optional[Path], use_n_carrier: bool):
        argv = ['-i', str(cat_result_path), '-s', str(sample_info_path), '-o_dir', str(out_dir)]
        if adj_factor_path is not None:
            argv.extend(['-a', str(adj_factor_path)])
        if use_n_carrier:
            argv.append('-u')
        tester = BinomialTest.get_instance(argv=argv)
        tester.result_path = Path(str(cat_result_path).replace('.categorization_result.txt.gz', '.burden_test.txt.gz'))
        
        if tester.use_n_carrier:
            tester.count_carrier_for_each_category()
            tester.calculate_relative_risk_with_n_carrier()
        else:
            tester.count_variant_for_each_category()
            tester.calculate_relative_risk()

        tester.run_burden_test()
        tester.concat_category_info()
        tester.save_result()
    
    def run(self):
        if self.step == 'annotation':
            self.multi_annotate()
        if self.cws_process == 'categorization':
            self.multi_categorize()
        if self.step == 'binomial_test':
            self.multi_binomial_test()
        if self.step == 'all':
            self.step = 'annotation'
            self.multi_annotate()
            self.in_dir = self.out_dir
            self.out_dir = Path('/'.join([self.out_dir.parent, "random-mutations_cat"]))
            self.step = 'categorization'
            self.multi_categorize()
            self.in_dir = self.out_dir
            self.out_dir = Path('/'.join([self.out_dir.parent, "random-mutations_burden"]))
            self.step = 'binomial_test'
            self.multi_binomial_test()
        log.print_progress("Done")
    