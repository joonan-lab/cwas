from __future__ import annotations

import argparse
from abc import ABC, abstractmethod
from importlib.resources import path

from cwas.utils.log import print_err


class Runnable(ABC):
    def __init__(self, args: argparse.Namespace):
        self.args = args

    def _assign_config_to_attr(self, attr_name: str, config_filename: str):
        """
        Assign the path of a configuration file in cwas.config as an attribute
        of the instance of this class. The type of the attribute is
        pathlib.Path.
        """
        try:
            with path('cwas.config', config_filename) as p:
                setattr(self, attr_name, p)
        except FileNotFoundError:
            print_err(f'The configuration file \'{config_filename}\''
                      f'does not exist in cwas.config.')
            raise

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
