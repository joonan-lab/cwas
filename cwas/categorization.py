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
        self._annotated_vcf = None
        self._gene_matrix = None
        self._category_domain = None

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
    def result_path(self) -> Path:
        return (
            Path(self.get_env("CWAS_WORKSPACE")) / "categorization_result.txt"
        )

    @property
    def annotated_vcf(self) -> pd.DataFrame:
        if self._annotated_vcf is None:
            log.print_progress("Parse the annotated VCF")
            self._annotated_vcf = parse_annotated_vcf(
                Path(self.get_env("ANNOTATED_VCF"))
            )
        return self._annotated_vcf

    @property
    def gene_matrix(self) -> dict:
        if self._gene_matrix is None:
            log.print_progress("Parse the gene matrix")
            self._gene_matrix = parse_gene_matrix(
                Path(self.get_env("GENE_MATRIX"))
            )
        return self._gene_matrix

    @property
    def category_domain(self) -> dict:
        if self._category_domain is None:
            log.print_progress("Load the category domain")
            with Path(self.get_env("CATEGORY_DOMAIN")) as category_domain_file:
                self._category_domain = yaml.safe_load(category_domain_file)
        return self._category_domain

    def run(self):
        log.print_log("Notice", "Not implemented yet.")

    def categorize_vcf_for_each_sample(self):
        return [
            _categorize_variant(
                vcf_for_sample, self.category_domain, self.gene_matrix
            )
            for vcf_for_sample in self.split_annotated_vcf_by_sample()
        ]

    def split_annotated_vcf_by_sample(self) -> list:
        groupby_sample = self.annotated_vcf.groupby("SAMPLE")
        sample_ids = list(groupby_sample.groups)
        return [groupby_sample.get_group(sample_id) for sample_id in sample_ids]

    def update_env(self):
        self.set_env("CATEGORIZATION_RESULT", self.result_path)
        self.save_env()
