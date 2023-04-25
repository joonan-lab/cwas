import yaml, os, gzip, sys, argparse, pickle
import multiprocessing as mp
from functools import partial
import re
from tqdm import tqdm
import numpy as np
import pandas as pd
from pathlib import Path
from cwas.runnable import Runnable
from cwas.core.common import cmp_two_arr
from cwas.core.categorization.parser import parse_annotated_vcf
import cwas.utils.log as log
from cwas.core.simulation.fastafile import FastaFile
from cwas.core.simulation.randomize import label_variant, pick_mutation
from cwas.utils.check import check_is_file, check_num_proc, check_same_n_lines
from typing import Optional

class Simulation(Runnable):
    def __init__(self, args: argparse.Namespace):
        super().__init__(args)
        self._in_vcf = None
        self._sample_info = None
        self._fam_to_label_cnt = None
        self._fam_to_sample_set = None
        self._filepath_dict = None
        self._chrom_size_df = None
        self._fa_file_dict = None
        self._rand_mut_paths = None
        self._resume = None
    
    @staticmethod
    def _create_arg_parser() -> argparse.ArgumentParser:
        parser = argparse.ArgumentParser(
            description="Arguments of CWAS simulation step",
            formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        )
        parser.add_argument('-i', '--in_vcf', dest='in_vcf_path', required=True, type=Path,
                            help='Input VCF file which is referred to generate random mutations')
        parser.add_argument('-s', '--sample_info', dest='sample_info_path', required=True, type=Path,
                            help='File listing sample IDs with their families and sample_types (case or ctrl)')
        parser.add_argument('-o', '--out_dir', dest='out_dir', required=False, type=Path,
                            help='Directory of outputs that lists random mutations. '
                                'The number of outputs will be the same with the number of simulations. '
                                '(default: $CWAS_WORKSPACE/random-mutation)')
        parser.add_argument('-t', '--out_tag', dest='out_tag', required=False, type=str,
                            help='Prefix of output files. Each output file name will start with this tag.', default='rand_mut')
        parser.add_argument('-n', '--num_sim', dest='num_sim', required=False, type=int,
                            help='Number of simulations to generate random mutations', default=1)
        parser.add_argument('-p', '--num_proc', dest='num_proc', required=False, type=int,
                            help='Number of processes for this script (only necessary for split VCF files)', default=1)
        parser.add_argument(
            "-r", "--resume",
            dest="resume", required=False, default=False, action="store_true",
            help="Resume the simulation from the last step. Assume some generated output files are not truncated.",
        )
        
    @staticmethod
    def _print_args(args: argparse.Namespace):
        log.print_arg("Input VCF file", args.in_vcf_path if args.in_vcf_path else "Not specified: $ANNOTATED_VCF will be used")
        log.print_arg("Sample information file", args.sample_info_path)
        log.print_arg("Output directory", args.out_dir)
        log.print_arg("Output tag (prefix of output files)", args.out_tag)
        log.print_arg("Number of simulations", args.num_sim)
        log.print_arg("File listing adjustment factors of each sample", args.adj_factor_path)
        log.print_arg(
            "No. worker processes for simulations",
            f"{args.num_proc: ,d}",
        )

    @staticmethod
    def _check_args_validity(args: argparse.Namespace):
        check_is_file(args.sample_info_path)
        check_num_proc(args.num_proc)
        if args.in_vcf_path:
            check_is_file(args.in_vcf_path)
        if args.adj_factor_path is not None:
            check_is_file(args.adj_factor_path)
            
    @property
    def in_vcf_path(self) -> Path:
        return (
            self.args.in_vcf_path.resolve()
            if self.args.in_vcf_path
            else Path(self.get_env("ANNOTATED_VCF"))
        )
    
    @property
    def sample_info_path(self) -> Path:
        return self.args.sample_info_path.resolve()
    
    @property
    def out_dir(self) -> Path:
        return (
            self.args.out_dir.resolve()
            if self.args.out_dir
            else self.workspace / "random-mutations"
        )

    @property
    def out_tag(self) -> str:
        return self.args.out_tag

    @property
    def num_sim(self) -> int:
        return self.args.num_sim
    
    @property
    def num_proc(self) -> int:
        return self.args.num_proc

    @property
    def resume(self) -> bool:
        if self._resume is None:
            self._resume = self.args.resume
        return self._resume

    @property
    def adj_factor_path(self) -> Optional[Path]:
        return (
            self.args.adj_factor_path.resolve()
            if self.args.adj_factor_path
            else None
        )
        
    @property
    def in_vcf(self) -> pd.DataFrame:
        if self._in_vcf is None:
            self._in_vcf = parse_annotated_vcf(
                self.in_vcf_path
            )[["REF", "ALT", "SAMPLE"]] 

        return self._in_vcf

    @property
    def sample_info(self) -> pd.DataFrame:
        if self._sample_info is None:
            self._sample_info = pd.read_table(
                self.sample_info_path, index_col="SAMPLE"
            )
        return self._sample_info

    @property
    def fam_to_label_cnt(self) -> dict:
        if self._fam_to_label_cnt is None:
            sample_to_fam = self.sample_info.to_dict()['FAMILY']
            variant_labels = np.vectorize(label_variant)(self.in_vcf.REF.values, self.in_vcf.ALT.values)
            
            fam_to_label_cnt = {}
            for sample, variant_label in zip(self.in_vcf["SAMPLE"].values, variant_labels):
                family = sample_to_fam[sample]
                label_cnt_arr = fam_to_label_cnt.get(family, np.zeros(4, dtype=int))
                label_cnt_arr[variant_label] += 1
                fam_to_label_cnt[family] = label_cnt_arr

            self._fam_to_label_cnt = fam_to_label_cnt
            
        return self._fam_to_label_cnt

    @property
    def fam_to_sample_set(self) -> dict:
        if self._fam_to_sample_set is None:
            sample_to_fam = self.sample_info.to_dict()['FAMILY']
            
            fam_to_sample_set = {}
            for sample in self.in_vcf["SAMPLE"].values:
                family = sample_to_fam[sample]
                sample_set = fam_to_sample_set.get(family, set())
                sample_set.add(sample)
                fam_to_sample_set[family] = sample_set

            self._fam_to_sample_set = fam_to_sample_set
            
        return self._fam_to_sample_set
    
    @property
    def filepath_dict(self) -> dict:
        if self._filepath_dict is None:
            try:
                self.sim_data_dir = Path(self.get_env("SIMULATION_DATA_DIR"))
                self.annot_data_dir = Path(self.get_env("SIMULATION_PATHS"))
            except TypeError:
                raise RuntimeError(
                    "Failed to get one of CWAS environment variable."
                    " Maybe you omitted to run Configuration step."
                )
            with open(self.annot_data_dir) as filepath_conf_file:
                filepath_conf = yaml.safe_load(filepath_conf_file)
                filepath_dict = {path_key: (self.sim_data_dir / filepath_conf[path_key]) for path_key in filepath_conf}
            for filepath in filepath_dict.values(): check_is_file(filepath)

            self._filepath_dict = filepath_dict
            
        return self._filepath_dict

    @property
    def chrom_size_df(self) -> pd.DataFrame:
        if self._chrom_size_df is None:
            self._chrom_size_df = pd.read_table(
                self.filepath_dict["chrom_size"], index_col="Chrom"
            )
        return self._chrom_size_df
    
    @property
    def fa_file_dict(self) -> dict:
        if self._fa_file_dict is None:
            # Make a dictionary for FASTA files masked in the previous preparation step.
            unq_chroms = np.unique(self.chrom_size_df.index.values)
            fa_file_dict = {}

            for chrom in unq_chroms:
                fa_file_path = self.filepath_dict[f'{chrom}']
                fa_file_dict[chrom] = FastaFile(fa_file_path)
        
            self._fa_file_dict = fa_file_dict
            
        return self._fa_file_dict

    @property
    def rand_mut_paths(self) -> list:
        if self._rand_mut_paths is None:
            rand_mut_paths = []
            for i in range(self.num_sim):
                str_i = str(i+1).zfill(len(str(self.num_sim)))
                output_filename = f'{self.out_tag}.{str_i}.vcf.gz'
                output_path = self.out_dir / output_filename
                rand_mut_paths.append(output_path)
            
            self._rand_mut_paths = rand_mut_paths
            
        return self._rand_mut_paths

    def run(self):
        self.prepare()
        self.make_rand_mut_files()
        log.print_progress("Done")
    
    def prepare(self):
        if not cmp_two_arr(self.sample_info.index, self.in_vcf['SAMPLE'].unique()):
            log.print_warn("The sample IDs in the sample info file and the VCF file are not the same.")

        log.print_progress("Make an output directory")
        self.out_dir.mkdir(parents=True, exist_ok=True)

    def make_rand_mut_files(self):
        """Make VCF files listing random mutations"""
        log.print_progress(self.make_rand_mut_files.__doc__)
        pre_files = sorted(self.out_dir.glob(f'{self.out_tag}.{"?"*len(str(self.num_sim))}.vcf.gz'))
        if len(pre_files) == 0:
            target_files = self.rand_mut_paths
        elif len(pre_files) == self.num_sim:
            log.print_progress("Checking the number of lines for the last 100 VCFs...")
            check_same_n_lines(pre_files[-100:], gzip_file=True)
            log.print_log(
                "NOTICE",
                "You already have random mutation vcfs. Skip this step.",
                False,
            )
            return
        elif self.resume & (len(pre_files) < self.num_sim):
            log.print_progress("Checking the number of lines for the last 100 VCFs...")
            check_same_n_lines(pre_files[-100:], gzip_file=True)
            log.print_log(
                "NOTICE",
                f"You have some random mutation vcfs ({len(pre_files)}). Resume this step.",
                False,
            )
            target_files = sorted(list(set(self.rand_mut_paths) - set(pre_files)))
            self._resume = False
        else:
            raise RuntimeError(
                "The number of random mutation vcfs is not the same as the number of simulations."
                " Check and remove the files to rerun this step."
            )
            
        if self.num_proc == 1:
            for output_path in tqdm(target_files):
                self.make_rand_mut_file(output_path)
        else:
            with mp.Pool(self.num_proc) as pool:
                for _ in tqdm(pool.imap_unordered(self.make_rand_mut_file, target_files), total=len(target_files)):
                    pass
        
        for fasta_file in self.fa_file_dict.values():
            fasta_file.close()


    def make_rand_mut_file(self, output_path):
        """Make a VCF file listing random mutations """
            
        rand_variants = []

        for fam in self.fam_to_label_cnt:
            label_cnt_arr = self.fam_to_label_cnt[fam]
            sample_ids = list(self.fam_to_sample_set[fam])

            for label, label_cnt in enumerate(label_cnt_arr):
                for _ in range(label_cnt):
                    rand_variant = self.make_random_mutation(label, sample_ids)
                    rand_variants.append(rand_variant)

        rand_variants.sort(key=lambda x: (int(x.get('chrom').replace('chr', '')), x.get('pos')))
        
        with gzip.open(output_path, 'wt') as outfile:
            print('#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', sep='\t', file=outfile)

            for variant in rand_variants:
                print(variant['chrom'], variant['pos'], variant['id'], variant['ref'], variant['alt'],
                    variant['qual'], variant['filter'], variant['info'], sep='\t', file=outfile)
            
       
    def make_random_mutation(self, label: int, sample_ids: list) -> dict:
        """Generate and return a random mutation in the VCF format """
        
        chrom_eff_sizes = self.chrom_size_df['Effective'].values
        chrom_probs = chrom_eff_sizes / np.sum(chrom_eff_sizes)  # Normalization
        chrom_sizes = self.chrom_size_df['Size'].values

        np.random.seed()
        sample_id = np.random.choice(sample_ids)
        ref = None
        alt = None

        while True:
            chrom_idx = np.random.choice(range(len(chrom_probs)), p=chrom_probs)
            chrom_size = chrom_sizes[chrom_idx]
            chrom = f'chr{chrom_idx + 1}'
            fa_file = self.fa_file_dict[chrom]

            pos = np.random.randint(chrom_size)
            base = fa_file.get_base(chrom, pos).upper()

            if base == 'N':
                continue

            ref, alt = pick_mutation()

            if base != ref:
                continue

            break

        alt += 'A' * label

        variant = {
            'chrom': chrom,
            'pos': pos + 1,  # 0-based -> 1-based
            'id': f'{chrom}:{pos + 1}:{ref}:{alt}',
            'ref': ref,
            'alt': alt,
            'qual': '.',
            'filter': '.',
            'info': f'SAMPLE={sample_id}'
        }
        return variant
    