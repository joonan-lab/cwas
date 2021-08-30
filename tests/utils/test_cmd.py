"""
Test cwas.utils.cmd
"""
import random
import subprocess

import cwas.utils.cmd as cmd
import pytest

_RAND_N = random.randint(1, 1000000)


def test_execute_basic():
    args = ["ls", "-l"]
    output = cmd.execute(args)
    assert isinstance(output, subprocess.CompletedProcess)
    assert output.returncode == 0
    assert output.args == args


def test_execute_fail():
    fake_filename = f"cwas_fake_{_RAND_N}"
    args = ["ls", "-l", fake_filename]
    output = cmd.execute(args, raise_err=False)
    assert output.returncode != 0

    with pytest.raises(subprocess.CalledProcessError):
        cmd.execute(args)


def test_execute_bin():
    binary = "ls"
    args = ["-l"]
    output = cmd.execute_bin(binary, args)
    assert isinstance(output, subprocess.CompletedProcess)
    assert output.returncode == 0


def test_execute_bin_fail():
    binary = "ls"
    fake_filename = f"cwas_fake_{_RAND_N}"
    args = ["-l", fake_filename]
    output = cmd.execute_bin(binary, args, raise_err=False)
    assert output.returncode != 0

    with pytest.raises(subprocess.CalledProcessError):
        cmd.execute_bin(binary, args)


def test_execute_bin_fake_bin():
    fake_binary = f"fake_bin_{_RAND_N}"
    args = ["-h"]

    with pytest.raises(FileNotFoundError):
        cmd.execute_bin(fake_binary, args)
