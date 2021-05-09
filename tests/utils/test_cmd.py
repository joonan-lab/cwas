"""
Test cwas.utils.cmd
"""
import subprocess
import random
import pytest

import cwas.utils.cmd as cmd


def test_execute_basic():
    args = ['ls', '-l']
    output = cmd.execute(args)
    assert isinstance(output, subprocess.CompletedProcess)
    assert output.returncode == 0
    assert output.args == args


def test_execute_fail():
    fake_filename = f'cwas_fake_{random.randint(1, 1000000)}'
    args = ['ls', '-l', fake_filename]
    output = cmd.execute(args, raise_err=False)
    assert output.returncode != 0

    with pytest.raises(subprocess.CalledProcessError):
        cmd.execute(args)
