"""
Factory method for Runnable
Ref: https://realpython.com/python-import/
"""
import argparse
import importlib

import cwas.argparser
from cwas.runnable import Runnable


class CWASFactory:
    def __init__(
        self, argparser: argparse.ArgumentParser, runnable: Runnable
    ) -> None:
        self._argparser = argparser
        self._runnable = runnable

    @property
    def argparser(self) -> argparse.ArgumentParser:
        return self._argparser

    @property
    def runnable(self) -> Runnable:
        return self._runnable


def create(step_name: str) -> CWASFactory:
    try:
        runnable_module = importlib.import_module(f"cwas.{step_name}")
        runnable = getattr(runnable_module, make_class_name(step_name))
        argparser = getattr(cwas.argparser, step_name)
    except (ImportError, AttributeError):
        raise ValueError(f'CWAS does not support "{step_name}".')

    return CWASFactory(argparser, runnable)


def make_class_name(module_name: str) -> str:
    class_name = module_name.title()
    class_name = "".join(class_name.split("_"))
    return class_name
