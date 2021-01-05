import argparse
import unittest
from unittest.mock import patch

from cwas.runnable import Runnable


class TestRunnable(unittest.TestCase):
    # Enable the Runnable class to be instantiated
    # by assign an empty set to Runnable.__abstractmethods__
    @patch.multiple(Runnable, __abstractmethods__=set())
    def test_get_instance(self):
        inst = Runnable.get_instance(argv=list())
        assert hasattr(inst, 'args')
        assert isinstance(inst.args, argparse.Namespace)

    @patch.multiple(Runnable, __abstractmethods__=set())
    def test_run(self):
        inst = Runnable.get_instance(argv=list())
        inst.run()
        assert 1
