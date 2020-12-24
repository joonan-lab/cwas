from __future__ import annotations

import argparse
from abc import ABC, abstractmethod


class Runnable(ABC):
    def __init__(self, args: argparse.Namespace):
        self.args = args

    @classmethod
    def get_instance(cls) -> Runnable:
        arg_parser = cls._create_arg_parser()
        args = arg_parser.parse_args()
        cls._print_args(args)
        cls._check_args_validity(args)
        return cls(args)

    @staticmethod
    @abstractmethod
    def _create_arg_parser() -> argparse.ArgumentParser:
        pass

    @staticmethod
    @abstractmethod
    def _print_args(args: argparse.Namespace):
        pass

    @staticmethod
    @abstractmethod
    def _check_args_validity(args: argparse.Namespace):
        pass

    @abstractmethod
    def run(self):
        pass
