import argparse
from pathlib import Path

import pandas as pd
import yaml

import cwas.utils.log as log
from cwas.core.categorization.parser import (
    parse_annotated_vcf,
    parse_gene_matrix,
)
from cwas.runnable import Runnable
from cwas.utils.check import check_num_proc


class Categorization(Runnable):
    def __init__(self, args: argparse.Namespace):
        super().__init__(args)

    @staticmethod
    def _create_arg_parser() -> argparse.ArgumentParser:
        parser = argparse.ArgumentParser(description=__doc__)
        parser.add_argument(
            "-p",
            "--num_proc",
            dest="num_proc",
            required=False,
            type=int,
            help="Number of worker processes for the categorization",
            default=1,
        )
        return parser

    @staticmethod
    def _print_args(args: argparse.Namespace):
        log.print_arg(
            "No. worker processes for the categorization",
            f"{args.num_proc: ,d}",
        )

    @staticmethod
    def _check_args_validity(args: argparse.Namespace):
        check_num_proc(args.num_proc)

    @property
    def annotated_vcf_path(self) -> Path:
        return Path(self.get_env("ANNOTATED_VCF"))

    @property
    def gene_matrix(self) -> Path:
        return Path(self.get_env("GENE_MATRIX"))

    @property
    def category_domain(self) -> Path:
        return Path(self.get_env("CATEGORY_DOMAIN"))

    @property
    def result_path(self) -> Path:
        return (
            Path(self.get_env("CWAS_WORKSPACE")) / "categorization_result.txt"
        )

    @property
    def variant_df(self) -> pd.DataFrame:
        log.print_progress("Parse the annotated VCF")
        return parse_annotated_vcf(self.annotated_vcf_path)

    @property
    def gene_matrix_dict(self) -> dict:
        log.print_progress("Parse the gene matrix")
        return parse_gene_matrix(self.gene_matrix)

    @property
    def category_domain_dict(self) -> dict:
        log.print_progress("Load the category domain")
        with self.category_domain as category_domain_file:
            return yaml.safe_load(category_domain_file)

    def run(self):
        log.print_log("Notice", "Not implemented yet.")

    def update_env(self):
        self.set_env("CATEGORIZATION_RESULT", self.result_path)
        self.save_env()
