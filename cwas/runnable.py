from __future__ import annotations

import argparse
from abc import ABC, abstractmethod


class Runnable(ABC):
    def __init__(self, args: argparse.Namespace):
        self.args = args
        self._set_env()
        self._check_env()

    @abstractmethod
    def _set_env(self):
        """
        Set paths of configuration files and resources to run.
        These paths are assigned to attributes of this instance.
        """
        pass

    @abstractmethod
    def _check_env(self):
        """
        This method must be defined to verify your environment.
        This method may include exception or error handling.
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
