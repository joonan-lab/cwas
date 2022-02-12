import argparse
from abc import abstractmethod
from pathlib import Path

import numpy as np
import pandas as pd

from cwas.core.categorization.category import Category
from cwas.core.common import cmp_two_arr
from cwas.runnable import Runnable
from cwas.utils.check import check_is_file
from cwas.utils.log import print_arg, print_progress


class BurdenTest(Runnable):
    def __init__(self, args: argparse.Namespace):
        super().__init__(args)
        self._sample_info = None
        self._adj_factor = None
        self._categorization_result = None
        self._result = None
        self._phenotypes = None
        self._case_variant_cnt = None
        self._ctrl_variant_cnt = None

    @staticmethod
    def _create_arg_parser() -> argparse.ArgumentParser:
        parser = argparse.ArgumentParser(
            description="Arguments of Burden Tests",
            formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        )
        parser.add_argument(
            "-s",
            "--sample_info",
            dest="sample_info_path",
            required=True,
            type=Path,
            help="File listing information of your samples",
        )
        parser.add_argument(
            "-a",
            "--adjustment_factor",
            dest="adj_factor_path",
            required=False,
            default=None,
            type=Path,
            help="File listing adjustment factors of each sample",
        )
        return parser

    @staticmethod
    def _print_args(args: argparse.Namespace):
        print_arg("Sample information file", args.sample_info_path)
        print_arg("Adjustment factor list", args.adj_factor_path)

    @staticmethod
    def _check_args_validity(args: argparse.Namespace):
        check_is_file(args.sample_info_path)
        if args.adj_factor_path is not None:
            check_is_file(args.adj_factor_path)

    @property
    def result_path(self) -> Path:
        return Path(
            self.get_env("ANNOTATED_VCF").replace(
                "annotated.vcf", "burden_test.txt"
            )
        )

    @property
    def sample_info(self) -> pd.DataFrame:
        if self._sample_info is None:
            self._sample_info = pd.read_table(
                self.sample_info_path, index_col="SAMPLE"
            )
        return self._sample_info

    @property
    def adj_factor(self) -> pd.DataFrame:
        if self.adj_factor_path is None:
            return None
        if self._adj_factor is None:
            self._adj_factor = pd.read_table(
                self.adj_factor_path, index_col="SAMPLE"
            )
        return self._adj_factor

    @property
    def categorization_result(self) -> pd.DataFrame:
        if self._categorization_result is None:
            print_progress("Load the categorization result")
            self._categorization_result = pd.read_table(
                self.get_env("CATEGORIZATION_RESULT"), index_col="SAMPLE"
            )
            if self.adj_factor is not None:
                self._adjust_categorization_result()
        return self._categorization_result

    def _adjust_categorization_result(self):
        if not _contain_same_index(
            self._categorization_result, self.adj_factor
        ):
            raise ValueError(
                "The sample IDs from the adjustment factor list are "
                "not the same with the sample IDs "
                "from the categorization result."
            )
        adj_factors = [
            self.adj_factor.to_dict()["AdjustFactor"][sample_id]
            for sample_id in self._categorization_result.index.values
        ]
        self._categorization_result = self._categorization_result.multiply(
            adj_factors, axis="index"
        )

    @property
    def phenotypes(self) -> np.ndarray:
        if self._phenotypes is None:
            self._phenotypes = np.vectorize(
                lambda sample_id: self.sample_info.to_dict()["PHENOTYPE"][
                    sample_id
                ]
            )(self.categorization_result.index.values)
        return self._phenotypes

    @property
    def case_variant_cnt(self) -> np.ndarray:
        if self._case_variant_cnt is None:
            self._case_variant_cnt = self.categorization_result.values[
                self.phenotypes == "case", :
            ].sum(axis=0)
        return self._case_variant_cnt

    @property
    def ctrl_variant_cnt(self) -> np.ndarray:
        if self._ctrl_variant_cnt is None:
            self._ctrl_variant_cnt = self.categorization_result.values[
                self.phenotypes == "ctrl", :
            ].sum(axis=0)
        return self._ctrl_variant_cnt

    @property
    def category_table(self) -> pd.DataFrame:
        categories = [
            Category.from_str(category_str).to_dict()
            for category_str in self.categorization_result.columns.values
        ]
        return pd.DataFrame(
            categories, index=self.categorization_result.columns.values
        )

    def run(self):
        if not _contain_same_index(
            self.categorization_result, self.sample_info
        ):
            raise ValueError(
                "The sample IDs from the sample information are "
                "not the same with the sample IDs "
                "from the categorization result."
            )
        self.count_variant_for_each_category()
        self.calculate_relative_risk()
        self.run_burden_test()
        self.concat_category_info()
        self.save_result()
        self.update_env()

    def count_variant_for_each_category(self):
        print_progress("Count the number of variants for each category")
        variant_cnt_arr = np.concatenate(
            [
                self.case_variant_cnt[:, np.newaxis],
                self.ctrl_variant_cnt[:, np.newaxis],
            ],
            axis=1,
        )
        self._result = pd.DataFrame(
            variant_cnt_arr,
            index=self.categorization_result.columns.values,
            columns=["Case_DNV_Count", "Ctrl_DNV_Count"],
        )

    def calculate_relative_risk(self):
        print_progress("Calculate relative risks for each category")
        self._result["Relative_Risk"] = (
            self.case_variant_cnt / self.ctrl_variant_cnt
        )

    @abstractmethod
    def run_burden_test(self):
        raise RuntimeError(
            "This method cannot be called via the instance of BurdenTest."
        )

    def concat_category_info(self):
        self._result = pd.concat([self.category_table, self._result], axis=1)
        self._result.index.name = "Category"

    def save_result(self):
        print_progress(f"Save the result to the file {self.result_path}")
        self._result.to_csv(self.result_path, sep="\t")

    def update_env(self):
        self.set_env("BURDEN_TEST_RESULT", self.result_path)
        self.save_env()


def _contain_same_index(table1: pd.DataFrame, table2: pd.DataFrame) -> bool:
    return cmp_two_arr(table1.index.values, table2.index.values)
