from __future__ import annotations

import argparse
from abc import ABC, abstractmethod


class Runnable(ABC):
    def __init__(self, args: argparse.Namespace):
        self.args = args

    @classmethod
    def get_instance(cls) -> Runnable:
        arg_parser = cls.create_arg_parser()
        args = arg_parser.parse_args()
        cls.print_args(args)
        cls.check_args_validity(args)
        return cls(args)

    @staticmethod
    @abstractmethod
    def create_arg_parser() -> argparse.ArgumentParser:
        pass

    @staticmethod
    @abstractmethod
    def print_args(args: argparse.Namespace):
        pass

    @staticmethod
    @abstractmethod
    def check_args_validity(args: argparse.Namespace):
        pass

    @abstractmethod
    def run(self):
        pass
