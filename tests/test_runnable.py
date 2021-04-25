import argparse
import pathlib
import unittest
import pytest
import random
from unittest.mock import patch

from cwas.runnable import Runnable


class RunnableMock(Runnable):
    @staticmethod
    def _create_arg_parser() -> argparse.ArgumentParser:
        parser = argparse.ArgumentParser()
        parser.add_argument('-d', dest='test_dir', type=pathlib.Path)
        parser.add_argument('-f', dest='test_file', type=pathlib.Path)
        parser.add_argument('-n', dest='test_int', type=int)
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

    @patch.multiple(Runnable, __abstractmethods__=set())
    def test__assign_config_to_attr(self):
        inst = Runnable.get_instance(argv=list())
        inst._assign_config_to_attr('gene_mat_path', 'gene_matrix.txt')
        assert hasattr(inst, 'gene_mat_path')
        assert isinstance(getattr(inst, 'gene_mat_path'), pathlib.Path)

        # Set a fake file
        with pytest.raises(FileNotFoundError):
            inst._assign_config_to_attr('non_exist', 'non_exist.txt')

    def test_argument_assignment(self):
        randint = random.randint(1, 1000000)
        test_dir = pathlib.Path.home() / f'.cwas-test-{randint}'
        test_dir.mkdir()
        test_file = test_dir / 'test.txt'
        test_file.touch()

        argv = ['-d', str(test_dir), '-f', str(test_file), '-n', str(randint)]
        inst = RunnableMock.get_instance(argv)

        assert hasattr(inst, 'test_dir')
        assert hasattr(inst, 'test_file')
        assert hasattr(inst, 'test_int')
        assert isinstance(getattr(inst, 'test_dir'), pathlib.Path)
        assert isinstance(getattr(inst, 'test_file'), pathlib.Path)
        assert getattr(inst, 'test_int') == randint
