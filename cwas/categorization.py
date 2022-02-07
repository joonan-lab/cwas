import argparse
from pathlib import Path

import pandas as pd
import yaml

import cwas.utils.log as log
from cwas.core.categorization.categorizer import Categorizer
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
        self._annotation_bw_cutoff = None
        self._annotated_vcf_groupby_sample = None
        self._sample_ids = None

    @staticmethod
    def _create_arg_parser() -> argparse.ArgumentParser:
        parser = argparse.ArgumentParser(
            description="Arguments of CWAS categorization step",
            formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        )
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
            with Path(self.get_env("CATEGORY_DOMAIN")).open(
                "r"
            ) as category_domain_file:
                self._category_domain = yaml.safe_load(category_domain_file)
        return self._category_domain

    @property
    def annotation_bw_cutoff(self) -> dict:
        if self._annotation_bw_cutoff is None:
            log.print_progress("Load the annotation BigWig cutoff")
            with Path(self.get_env("ANNOTATION_BW_CUTOFF")).open(
                "r"
            ) as annotation_bw_cutoff_file:
                self._annotation_bw_cutoff = yaml.safe_load(
                    annotation_bw_cutoff_file
                )
        return self._annotation_bw_cutoff

    @property
    def categorizer(self) -> Categorizer:
        categorizer = Categorizer(self.category_domain, self.gene_matrix)
        categorizer.phastcons_cutoff = self.annotation_bw_cutoff[
            "phastCons46way"
        ]
        categorizer.phylop_cutoff = self.annotation_bw_cutoff["phyloP46way"]
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

    def run(self):
        self.categorize_vcf()
        self.save_result()
        self.update_env()

    def categorize_vcf(self):
        log.print_progress("Categorize variants of the annotated VCF")
        self._result = pd.DataFrame(
            self.categorize_vcf_for_each_sample()
        ).fillna(0)
        self._result = self._result.astype(int)
        self._result["SAMPLE"] = self.sample_ids
        self._result.set_index("SAMPLE", inplace=True)

    def categorize_vcf_for_each_sample(self):
        return [
            self.categorizer.categorize_variant(sample_vcf)
            for sample_vcf in self.annotated_vcf_split_by_sample
        ]

    def save_result(self):
        log.print_progress(f"Save the result as the file {self.result_path}")
        self._result.to_csv(self.result_path, sep="\t")

    def update_env(self):
        self.set_env("CATEGORIZATION_RESULT", self.result_path)
        self.save_env()
