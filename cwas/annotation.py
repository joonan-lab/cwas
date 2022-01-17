"""
CWAS Annotation Step

This step annotate user's VCF file using annotation data specified 
in the CWAS configuration step. This step mainly uses 
Variant Effect Predictor (VEP) to annotate user's VCF file.
"""
import argparse
from pathlib import Path

import yaml

from cwas.runnable import Runnable
from cwas.utils.check import check_is_file, check_num_proc
from cwas.utils.log import print_arg, print_log


class Annotation(Runnable):
    def __init__(self, args: argparse.Namespace):
        super().__init__(args)

    @staticmethod
    def _create_arg_parser() -> argparse.ArgumentParser:
        parser = argparse.ArgumentParser(
            description="Arguments of CWAS annotation step",
            formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        )
        parser.add_argument(
            "-v",
            "--vcf_file",
            dest="vcf_path",
            required=True,
            type=Path,
            help="Target VCF file",
        )
        parser.add_argument(
            "-p",
            "--num_process",
            dest="num_proc",
            required=False,
            type=int,
            help="Number of worker processes to use",
            default=1,
        )
        return parser

    @staticmethod
    def _print_args(args: argparse.Namespace):
        print_arg("Target VCF file", args.vcf_path)
        print_arg("Number of worker processes", args.num_proc)

    @staticmethod
    def _check_args_validity(args: argparse.Namespace):
        check_is_file(args.vcf_path)
        check_num_proc(args.num_proc)

    def run(self):
        print_log("Notice", "Not implemented yet.")

    @property
    def bw_custom_annotations(self):
        with open(self.get_env("ANNOTATION_BW_KEY")) as infile:
            bw_custom_path_dict = yaml.safe_load(infile)
        annotation_data_dir = self.get_env("ANNOTATION_DATA")

        for bw_filename, bw_annotation_key in bw_custom_path_dict.items():
            yield (f"{annotation_data_dir}/{bw_filename}", bw_annotation_key)
