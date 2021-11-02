from __future__ import annotations

import argparse
from abc import ABC, abstractmethod
from pathlib import Path
from typing import Any, Optional

from cwas.env import Env


class Runnable(ABC):
    def __init__(self, args: argparse.Namespace):
        self.env = Env()
        arg_dict = vars(args)
        for arg in arg_dict:
            arg_val = arg_dict.get(arg)
            if isinstance(arg_val, Path):
                arg_val = arg_val.resolve()
            setattr(self, arg, arg_val)

    @classmethod
    def get_instance(cls, argv: list = []) -> Runnable:
        arg_parser = cls._create_arg_parser()
        args = arg_parser.parse_args(argv)
        cls._print_args(args)
        cls._check_args_validity(args)
        return cls(args)

    def set_env_path(self, path: Path):
        self.env.set_path(path)

    def set_env(self, env_key: str, env_value: Any):
        self.env.set_env(env_key, env_value)

    def get_env(self, env_key: str) -> Optional[str]:
        return self.env.get_env(env_key)

    def save_env(self):
        self.env.save()

    @staticmethod
    def _create_arg_parser() -> argparse.ArgumentParser:
        return argparse.ArgumentParser()

    @staticmethod
    def _print_args(args: argparse.Namespace):
        pass

    @staticmethod
    def _check_args_validity(args: argparse.Namespace):
        pass

    @abstractmethod
    def run(self):
        pass
