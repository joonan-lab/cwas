"""
Test the methods in cwas.utils.log
"""
from datetime import datetime

import cwas.utils.log as log


def get_curr_time():
    now = datetime.now()
    curr_time = now.strftime("%H:%M:%S %m/%d/%y")
    return curr_time


def test_print_log(capsys):
    test_log_type = "TEST"
    test_log_msg = "test"
    curr_time = get_curr_time()

    log.print_log(test_log_type, test_log_msg, True)
    captured = capsys.readouterr().err
    expected = f"[{curr_time}, {test_log_type}] {test_log_msg}\n"
    assert captured == expected

    log.print_log(test_log_type, test_log_msg, False)
    captured = capsys.readouterr().err
    expected = f"[{test_log_type}] {test_log_msg}\n"
    assert captured == expected

    log.print_log(test_log_type, test_log_msg)
    captured = capsys.readouterr().err
    expected = f"[{test_log_type}] {test_log_msg}\n"
    assert captured == expected


def test_print_arg(capsys):
    test_arg_name = "My argument"
    test_arg_val = "My value"
    log.print_arg(test_arg_name, test_arg_val)
    captured = capsys.readouterr().err
    expected = f"[ARG] {test_arg_name}: {test_arg_val}\n"
    assert captured == expected


def test_print_progress(capsys):
    test_msg = "My progress"
    curr_time = get_curr_time()
    log.print_progress(test_msg)
    captured = capsys.readouterr().err
    expected = f"[{curr_time}, PROGRESS] {test_msg}\n"
    assert captured == expected


def test_print_warn(capsys):
    test_msg = "My warning"
    curr_time = get_curr_time()
    log.print_warn(test_msg)
    captured = capsys.readouterr().err
    expected = f"[{curr_time}, WARNING] {test_msg}\n"
    assert captured == expected


def test_print_err(capsys):
    test_msg = "My error"
    curr_time = get_curr_time()
    log.print_err(test_msg)
    captured = capsys.readouterr().err
    expected = f"[{curr_time}, ERROR] {test_msg}\n"
    assert captured == expected
