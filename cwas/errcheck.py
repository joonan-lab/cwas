"""
Functions to check conditions that cause cumbersome errors in this project.
Those conditions should be checked to prevent this project from being shut down
unexpectedly.
"""
import os


def check_is_file(file_path: str):
    if not file_path or not os.path.isfile(file_path):
        raise FileNotFoundError(f'"{file_path}" is not a valid file path.')


def check_is_dir(dir_path: str):
    if dir_path and not os.path.isdir(dir_path):
        raise NotADirectoryError(f'"{dir_path}" is not a valid directory.')
