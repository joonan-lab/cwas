from __future__ import annotations

import argparse
from abc import ABC, abstractmethod
from importlib.resources import path
from pathlib import Path

from cwas.utils.log import print_err


class Runnable(ABC):
    def __init__(self, args: argparse.Namespace):
        arg_dict = vars(args)
        for arg in arg_dict:
            arg_val = arg_dict.get(arg)
            if isinstance(arg_val, Path):
                arg_val = arg_val.resolve()
            setattr(self, arg, arg_val)

    @classmethod
    def get_instance(cls, argv: list = None) -> Runnable:
        arg_parser = cls._create_arg_parser()
        args = arg_parser.parse_args(argv)
        cls._print_args(args)
        cls._check_args_validity(args)
        return cls(args)

    @staticmethod
    @abstractmethod
    def _create_arg_parser() -> argparse.ArgumentParser:
        return argparse.ArgumentParser()

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
