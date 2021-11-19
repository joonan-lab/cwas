"""
Tests the situation when we call CWAS Start again
"""
from pathlib import Path

import pytest
from cwas.env import Env
from cwas.start import Start


@pytest.fixture(scope="module")
def args(cwas_workspace: Path) -> list:
    return ["-w", str(cwas_workspace)]


@pytest.fixture(scope="module")
def config() -> dict:
    test_config = {"TEST": "HELLO"}
    return test_config


@pytest.fixture(scope="module")
def env() -> dict:
    test_env = {"FAKE_ENV": "HI"}
    return test_env


@pytest.fixture(scope="module")
def create_config_early(cwas_workspace: Path, config: dict):
    cwas_workspace.mkdir()
    cwas_config_path = cwas_workspace / "configuration.txt"

    with cwas_config_path.open("w") as config_file:
        for config_key, config_value in config.items():
            print(f"{config_key}={config_value}", file=config_file)


@pytest.fixture(scope="module")
def create_env_early(cwas_env_path: Path, env: dict):
    with cwas_env_path.open("w") as env_file:
        for env_key, env_value in env.items():
            print(f"{env_key}={env_value}", file=env_file)


@pytest.fixture(scope="module", autouse=True)
def setup(args: list, create_config_early, create_env_early):
    inst = Start.get_instance(args)
    inst.run()


@pytest.fixture(scope="module", autouse=True)
def teardown(cwas_workspace: Path, cwas_env_path: Path):
    yield
    cwas_env_path.unlink()
    for f in cwas_workspace.glob("*"):
        f.unlink()
    cwas_workspace.rmdir()


def test_cwas_config(cwas_workspace: Path, config: dict):
    cwas_config_path = cwas_workspace / "configuration.txt"
    actual_config = get_config_from_file(cwas_config_path)
    assert actual_config.items() == config.items()


def test_cwas_env_path(cwas_env_path: Path, env: dict):
    actual_env = get_config_from_file(cwas_env_path)
    expected_key_set = {
        "CWAS_WORKSPACE",
        "FAKE_ENV",
    }
    assert actual_env.items() >= env.items()
    assert set(actual_env.keys()) == expected_key_set


def get_config_from_file(config_path: Path) -> dict:
    config_from_file = {}
    with config_path.open() as config_file:
        for line in config_file:
            k, v = line.strip().split("=")
            config_from_file[k] = v
    return config_from_file
