import argparse

from cwas.runnable import Runnable
from cwas.utils.log import print_log


class Configuration(Runnable):
    def __init__(self, args: argparse.Namespace):
        super().__init__(args)

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
        print_log('LOG', 'CWAS Configuration')
        print_log('Notice', 'Not implemented yet.')
