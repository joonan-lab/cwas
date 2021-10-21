"""
Tests of the 'Start' class
"""
import pytest
from cwas.start import Start


@pytest.fixture(scope="module")
def args(cwas_workspace):
    return ["-w", cwas_workspace]


def test_start_step(args, cwas_workspace):
    inst = Start.get_instance(args)
    inst.run()
    check_initial_file_exist(cwas_workspace)


def check_initial_file_exist(cwas_workspace):
    cwas_env_file_path = cwas_workspace / ".env"
    cwas_config_file_path = cwas_workspace / "configuration.txt"
    assert cwas_workspace.isdir()
    assert cwas_env_file_path.isfile()
    assert cwas_config_file_path.isfile()
