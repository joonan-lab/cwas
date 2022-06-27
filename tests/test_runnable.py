import argparse
import pathlib
import random
import unittest
from unittest.mock import patch

from cwas.runnable import Runnable


class RunnableMock(Runnable):
    @staticmethod
    def _create_arg_parser() -> argparse.ArgumentParser:
        parser = argparse.ArgumentParser()
        parser.add_argument("-d", dest="test_dir", type=pathlib.Path)
        parser.add_argument("-f", dest="test_file", type=pathlib.Path)
        parser.add_argument("-n", dest="test_int", type=int)
        return parser

    @staticmethod
    def _check_args_validity(args: argparse.Namespace):
        pass

    @staticmethod
    def _print_args(args: argparse.Namespace):
        pass

    def run(self):
        pass


class TestRunnable(unittest.TestCase):
    # Enable the Runnable class to be instantiated
    # by assign an empty set to Runnable.__abstractmethods__
    @patch.multiple(Runnable, __abstractmethods__=set())
    def test_get_instance(self):
        inst = Runnable.get_instance(argv=list())
        assert isinstance(inst, Runnable)

    @patch.multiple(Runnable, __abstractmethods__=set())
    def test_run(self):
        inst = Runnable.get_instance(argv=list())
        inst.run()
        assert 1

    def test_argument_assignment(self):
        randint = random.randint(1, 1000000)
        test_dir = pathlib.Path.home() / f".cwas-runnable-test-{randint}"
        test_dir.mkdir()
        test_file = test_dir / "test.txt"
        test_file.touch()

        argv = ["-d", str(test_dir), "-f", str(test_file), "-n", str(randint)]
        inst = RunnableMock.get_instance(argv)
        assert hasattr(inst, "args")

        # Teardown
        test_file.unlink()
        test_dir.rmdir()
