"""
Test cwas.utils.cmd
"""
import random
import subprocess

import cwas.utils.cmd as cmd
import pytest

_RAND_N = random.randint(1, 1000000)


def test_cmd_executor():
    bin_name = "ls"
    args = ["-l"]
    inst = cmd.CmdExecutor(bin_name, args)
    return_code = inst.execute()
    assert return_code == 0


def test_cmd_executor_with_fail():
    fake_filename = f"cwas_fake_{_RAND_N}"
    bin_name = "ls"
    args = ["-l", fake_filename]
    inst = cmd.CmdExecutor(bin_name, args)
    return_code = inst.execute()
    assert return_code != 0


def test_cmd_executor_raising_err():
    fake_filename = f"cwas_fake_{_RAND_N}"
    bin_name = "ls"
    args = ["-l", fake_filename]
    inst = cmd.CmdExecutor(bin_name, args)

    with pytest.raises(subprocess.CalledProcessError):
        inst.execute_raising_err()


def test_cmd_executor_with_invalid_bin():
    with pytest.raises(FileNotFoundError):
        cmd.CmdExecutor("not_exists")


@pytest.fixture(scope="function")
def bin_path(cwas_workspace):
    cwas_workspace.mkdir()
    _bin_path = cwas_workspace / "my_bin"
    _bin_path.touch()
    yield _bin_path
    _bin_path.unlink()
    cwas_workspace.rmdir()


def test_cmd_executor_with_file(bin_path):
    assert cmd.CmdExecutor(bin_path).bin_path == str(bin_path)
