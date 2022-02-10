import argparse
from abc import abstractmethod
from pathlib import Path

from cwas.runnable import Runnable
from cwas.utils.check import check_is_file
from cwas.utils.log import print_arg, print_log


class BurdenTest(Runnable):
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

    def run(self):
        print_log("Notice", "Not implemented yet.")

    @abstractmethod
    def test(self):
        pass
