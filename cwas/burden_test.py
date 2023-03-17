import argparse
from abc import abstractmethod
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd

from cwas.core.categorization.category import Category
from cwas.core.common import cmp_two_arr
from cwas.runnable import Runnable
from cwas.utils.check import check_is_file
from cwas.utils.check import check_is_dir
from cwas.utils.log import print_arg, print_progress


class BurdenTest(Runnable):
    def __init__(self, args: Optional[argparse.Namespace] = None):
        super().__init__(args)
        self._sample_info = None
        self._adj_factor = None
        self._categorization_result = None
        self._result = None
        self._phenotypes = None
        self._case_variant_cnt = None
        self._ctrl_variant_cnt = None
        self._case_carrier_cnt = None
        self._ctrl_carrier_cnt = None

    @staticmethod
    def _create_arg_parser() -> argparse.ArgumentParser:
        parser = argparse.ArgumentParser(
            description="Arguments of Burden Tests",
            formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        )
        default_workspace = Path.home() / ".cwas"
        parser.add_argument(
            "-i",
            "--input_file",
            dest="cat_path",
            required=True,
            type=Path,
            help="Categorized file (gzipped)",
        )
        parser.add_argument(
            "-o_dir",
            "--output_directory",
            dest="output_dir_path",
            required=False,
            default=default_workspace,
            type=Path,
            help="Directory where output file will be saved",
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
        parser.add_argument(
            "-u",
            "--use_n_carrier",
            dest="use_n_carrier",
            required=False,
            default=False,
            action="store_true",
            help="Use the number of samples with variants in each category for burden test instead of the number of variants",
        )
        return parser

    @staticmethod
    def _print_args(args: argparse.Namespace):
        print_arg("Categorized file", args.cat_path)
        print_arg("Output directory", args.output_dir_path)
        print_arg("Sample information file", args.sample_info_path)
        print_arg("Adjustment factor list", args.adj_factor_path)
        print_arg("If the number of carrier is used for burden test or not", args.use_n_carrier)

    @staticmethod
    def _check_args_validity(args: argparse.Namespace):
        check_is_file(args.cat_path)
        check_is_dir(args.output_dir_path)
        check_is_file(args.sample_info_path)
        if args.adj_factor_path is not None:
            check_is_file(args.adj_factor_path)

    @property
    def cat_path(self) -> Path:
        return self.args.cat_path.resolve()

    @property
    def output_dir_path(self) -> Path:
        return self.args.output_dir_path.resolve()

    @property
    def sample_info_path(self) -> Path:
        return self.args.sample_info_path.resolve()

    @property
    def adj_factor_path(self) -> Optional[Path]:
        return (
            self.args.adj_factor_path.resolve()
            if self.args.adj_factor_path
            else None
        )

    @property
    def result_path(self) -> Path:
        return Path(
            f"{self.output_dir_path}/"
            f"{self.cat_path.name.replace('categorization_result.txt', 'burden_test.txt')}"
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
        if self._adj_factor is None and self.adj_factor_path:
            self._adj_factor = pd.read_table(
                self.adj_factor_path, index_col="SAMPLE"
            )
        return self._adj_factor

    @property
    def categorization_result(self) -> pd.DataFrame:
        if self._categorization_result is None:
            print_progress("Load the categorization result")
            self._categorization_result = pd.read_table(
                self.cat_path, index_col="SAMPLE", compression='gzip'
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
    def case_cnt(self) -> int:
        return (self.phenotypes == "case").sum()

    @property
    def ctrl_cnt(self) -> int:
        return (self.phenotypes == "ctrl").sum()

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
    def case_carrier_cnt(self) -> np.ndarray:
        if self._case_carrier_cnt is None:
            var_counts = self.categorization_result.values
            is_carrier = np.where(var_counts > 0, 1, 0)
            self._case_carrier_cnt = is_carrier[
                self.phenotypes == "case", :
            ].sum(axis=0)
        return self._case_carrier_cnt

    @property
    def ctrl_carrier_cnt(self) -> np.ndarray:
        if self._ctrl_carrier_cnt is None:
            var_counts = self.categorization_result.values
            is_carrier = np.where(var_counts > 0, 1, 0)
            self._ctrl_carrier_cnt = is_carrier[
                self.phenotypes == "ctrl", :
            ].sum(axis=0)
        return self._ctrl_carrier_cnt

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
        if self.use_n_carrier:
            self.count_carrier_for_each_category()
            self.calculate_relative_risk_with_n_carrier()
        else:
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

    def count_carrier_for_each_category(self):
        print_progress("Count the number of carriers in each category")
        carrier_cnt_arr = np.concatenate(
            [
                self.case_carrier_cnt[:, np.newaxis],
                self.ctrl_carrier_cnt[:, np.newaxis],
            ],
            axis=1,
        )
        self._result = pd.DataFrame(
            carrier_cnt_arr,
            index=self.categorization_result.columns.values,
            columns=["Case_Carrier_Count", "Ctrl_Carrier_Count"],
        )

    def calculate_relative_risk(self):
        print_progress("Calculate relative risks for each category")
        normalized_case_variant_cnt = self.case_variant_cnt / self.case_cnt
        normalized_ctrl_variant_cnt = self.ctrl_variant_cnt / self.ctrl_cnt
        self._result["Relative_Risk"] = (
            normalized_case_variant_cnt / normalized_ctrl_variant_cnt
        )

    def calculate_relative_risk_with_n_carrier(self):
        print_progress("Calculate relative risks for each category")
        case_carrier_rate = self.case_carrier_cnt / self.case_cnt
        ctrl_carrier_rate = self.ctrl_carrier_cnt / self.ctrl_cnt
        self._result["Relative_Risk"] = (
            case_carrier_rate / ctrl_carrier_rate
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
