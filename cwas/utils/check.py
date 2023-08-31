"""
Functions to check common conditions that cause cumbersome errors in this
project. Those conditions should be checked to prevent this project from
being shut down unexpectedly.
"""
import multiprocessing as mp
from pathlib import Path
from typing import Union
import gzip


def check_is_file(file_path: Union[Path, str]):
    if not file_path:
        raise ValueError("The file path cannot be None.")

    if not isinstance(file_path, Path):
        file_path = Path(file_path)

    if not file_path.is_file():
        raise FileNotFoundError(
            f'"{file_path}" cannot be found. Check this file path first.'
        )


def check_is_dir(dir_path: Union[Path, str]):
    if not dir_path:
        raise ValueError(f"dir_path: {dir_path}")

    if not isinstance(dir_path, Path):
        dir_path = Path(dir_path)

    if not dir_path.is_dir():
        raise NotADirectoryError(
            f'"{dir_path}" is not a valid directory. Check this path first.'
        )


def check_num_proc(num_proc: int):
    max_num_proc = mp.cpu_count()

    if num_proc < 1:
        raise ValueError(
            f"Your number of worker processes '{num_proc}' is invalid. "
            f"This number must be a positive integer."
        )
    elif num_proc > max_num_proc:
        raise ValueError(
            f"Your number of worker processes '{num_proc}' exceeds the "
            f"number of available CPUs in your computer ({max_num_proc})."
        )

def check_same_n_lines(file_paths: list, gzip_file: bool = False):
    n_lines_set = set()
    for file_path in file_paths:
        f = gzip.open(file_path, "r") if gzip_file else open(file_path, "r")
        for count, _ in enumerate(f):
                pass
        f.close()
        n_lines_set.add(count)
    if len(n_lines_set) > 1:
        raise ValueError(
            f"Some files in the previous step have different number of lines. "
            f"Check the files first."
        )
