import argparse
from abc import abstractmethod
from pathlib import Path

import pandas as pd

from cwas.runnable import Runnable
from cwas.utils.check import check_is_file
from cwas.utils.log import print_arg, print_log, print_progress


class BurdenTest(Runnable):
    def __init__(self, args: argparse.Namespace):
        super().__init__(args)
        self._categorization_result = None
        self._sample_info = None

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
    def categorization_result(self) -> pd.DataFrame:
        if self._categorization_result is None:
            print_progress("Load the categorization result")
            self._categorization_result = pd.read_table(
                self.get_env("CATEGORIZATION_RESULT"), index_col="SAMPLE"
            )
        return self._categorization_result

    @property
    def sample_info(self) -> pd.DataFrame:
        if self._sample_info is None:
            self._sample_info = pd.read_table(
                self.sample_info_path, index_col="SAMPLE"
            )
        return self._sample_info

    def run(self):
        print_log("Notice", "Not implemented yet.")

    @abstractmethod
    def test(self):
        raise RuntimeError(
            "This method cannot be called via the instance of BurdenTest."
        )
