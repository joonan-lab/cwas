import argparse
import pathlib
import unittest
import pytest
from unittest.mock import patch

from cwas.runnable import Runnable


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
