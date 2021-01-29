"""
Functions to check common conditions that cause cumbersome errors in this
project. Those conditions should be checked to prevent this project from
being shut down unexpectedly.
"""
import multiprocessing as mp
import os


def check_is_file(file_path: str):
    if not file_path or not os.path.isfile(file_path):
        raise FileNotFoundError(
            f'"{file_path}" is not a valid file path. Check this path first.'
        )


def check_is_dir(dir_path: str):
    if dir_path and not os.path.isdir(dir_path):
        raise NotADirectoryError(
            f'"{dir_path}" is not a valid directory. Check this path first.'
        )


def check_num_proc(num_proc: int):
    max_num_proc = mp.cpu_count()

    if num_proc < 1:
        raise ValueError(
            f'Your number of worker processes \'{num_proc}\' is invalid. '
            f'This number must be a positive integer.'
        )
    elif num_proc > max_num_proc:
        raise ValueError(
            f'Your number of worker processes \'{num_proc}\' exceeds the '
            f'number of available CPUs in your computer ({max_num_proc}).'
        )
