"""
Functions for printing logs for this project.
"""
import sys
from datetime import datetime


def print_log(log_type: str, log_msg: str, print_time: bool = False):
    log_header = (
        f"[{_get_curr_time()}, {log_type}]" if print_time else f"[{log_type}]"
    )
    _eprint(f"{log_header} {log_msg}")


def print_arg(arg_name: str, arg_val: str):
    msg = f"{arg_name}: {arg_val}"
    print_log("ARG", msg)


def print_progress(msg: str):
    print_log("PROGRESS", msg, True)


def print_warn(msg: str):
    print_log("WARNING", msg, True)


def print_err(msg: str):
    print_log("ERROR", msg, True)


# Ref:
# https://stackoverflow.com/questions/5574702/how-to-print-to-stderr-in-python
def _eprint(*args, **kwargs):
    """
    print to the stderr
    """
    print(*args, file=sys.stderr, **kwargs)


def _get_curr_time() -> str:
    now = datetime.now()
    curr_time = now.strftime("%H:%M:%S %m/%d/%y")
    return curr_time
