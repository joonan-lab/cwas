from __future__ import annotations

import argparse
from abc import ABC, abstractmethod

from cwas.utils.log import print_err


class Runnable(ABC):
    def __init__(self, args: argparse.Namespace):
        self.args = args
        try:
            self._set_attr()
        except FileNotFoundError:
            print_err('One of configuration files or resources does not exist.')
            raise

    @abstractmethod
    def _set_attr(self):
        """
        Set attributes of the instance of this class.
        """
        pass



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
