import argparse

from runnable import Runnable


class BurdenTest(Runnable):
    @staticmethod
    def _create_arg_parser() -> argparse.ArgumentParser:
        pass

    @staticmethod
    def _print_args(args: argparse.Namespace):
        pass

    @staticmethod
    def _check_args_validity(args: argparse.Namespace):
        pass

    def run(self):
        pass
