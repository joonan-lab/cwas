import argparse
from abc import ABC, abstractmethod
from pathlib import Path

import cwas.utils.check as check
import cwas.utils.log as log


class CWASArgumentManager(ABC):
    """Abstract class of an argument manager for each CWAS module"""

    def __init__(self) -> None:
        super().__init__()
        self._args = None
        self._parser = self._create_argument_parser()

    @staticmethod
    @abstractmethod
    def _create_argument_parser() -> argparse.ArgumentParser:
        pass

    def parse_args(self, argv: list = []) -> argparse.Namespace:
        result = self._parser.parse_args(argv)
        self._check_args(result)
        self._args = result
        return self._args

    @staticmethod
    def _check_args(args: argparse.Namespace):
        pass

    def print_args(self):
        pass


class CategorizationArgumentManager(CWASArgumentManager):
    def _create_argument_parser() -> argparse.ArgumentParser:
        result = argparse.ArgumentParser(
            description="Arguments of CWAS categorization step",
            formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        )
        result.add_argument(
            "-i",
            "--input_file",
            dest="input_path",
            required=True,
            type=Path,
            help="Annotated VCF file",
        )
        result.add_argument(
            "-p",
            "--num_proc",
            dest="num_proc",
            required=False,
            type=int,
            help="Number of worker processes for the categorization",
            default=1,
        )
        return result

    @staticmethod
    def _check_args(args: argparse.Namespace):
        check.check_num_proc(args.num_proc)
        check.check_is_file(args.input_path)

    def print_args(self):
        if self._args:
            log.print_arg(
                "No. worker processes for the categorization",
                f"{self._args.num_proc: ,d}",
            )
            log.print_arg("Annotated VCF file", self._args.input_path)
        else:
            log.print_warn("No parsed arguments")
