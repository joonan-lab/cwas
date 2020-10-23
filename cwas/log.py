"""
Functions for printing logs for this project.
"""
import sys
from datetime import datetime


def print_arg(arg_name: str, arg_val: str):
    _eprint(f'[ARG] {arg_name}: {arg_val}')


def print_progress(msg: str):
    _eprint(f'[{_get_curr_time()}, Progress] {msg}')


def print_warn(msg: str):
    _eprint(f'[{_get_curr_time()}, WARNING] {msg}')


def print_err(msg: str):
    _eprint(f'[{_get_curr_time()}, ERROR] {msg}')


# Ref:
# https://stackoverflow.com/questions/5574702/how-to-print-to-stderr-in-python
def _eprint(*args, **kwargs):
    """
    print to the stderr
    """
    print(*args, file=sys.stderr, **kwargs)


def _get_curr_time() -> str:
    now = datetime.now()
    curr_time = now.strftime('%H:%M:%S %m/%d/%y')
    return curr_time
