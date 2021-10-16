"""
CWAS Annotation Step

This step annotate user's VCF file using annotation data specified 
in the CWAS configuration step. This step mainly uses 
Variant Effect Predictor (VEP) to annotate user's VCF file.
"""
import argparse
from pathlib import Path

from cwas.runnable import Runnable
from cwas.utils.error import check_is_file, check_num_proc
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
