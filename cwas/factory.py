"""
Factory method for Runnable
Ref: https://realpython.com/python-import/
"""
import importlib

from cwas.runnable import Runnable


def make_class_name(module_name: str) -> str:
    class_name = module_name.title()
    class_name = "".join(class_name.split("_"))
    return class_name


def get_runnable(module_name: str) -> Runnable:
    try:
        module = importlib.import_module(f"cwas.{module_name}")
        runnable = getattr(module, make_class_name(module_name))
    except (ImportError, AttributeError):
        raise ValueError(f'CWAS does not support "{module_name}".')

    return runnable
