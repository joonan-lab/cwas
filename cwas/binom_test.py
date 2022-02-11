import argparse

from cwas.burden_test import BurdenTest


class BinomTest(BurdenTest):
    @staticmethod
    def _create_arg_parser() -> argparse.ArgumentParser:
        pass

    @staticmethod
    def _print_args(args: argparse.Namespace):
        pass

    @staticmethod
    def _check_args_validity(args: argparse.Namespace):
        pass

    def run_burden_test(self):
        pass
