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
def run_start(
    args: list, cwas_env_path: Path, create_config_early, create_env_early
):
    inst = Start.get_instance(args)
    inst.set_env_path(cwas_env_path)
    inst.env = Env()  # Reload
    inst.run()
    yield


def test_cwas_config(cwas_workspace: Path, config: dict):
    cwas_config_path = cwas_workspace / "configuration.txt"
    actual_config = get_config_from_file(cwas_config_path)
    assert actual_config.items() >= config.items()


def test_cwas_env_path(cwas_env_path: Path, env: dict):
    actual_env = get_config_from_file(cwas_env_path)
    assert actual_env.items() >= env.items()


def get_config_from_file(config_path: Path) -> dict:
    config_from_file = {}
    with config_path.open() as config_file:
        for line in config_file:
            k, v = line.strip().split("=")
            config_from_file[k] = v
    return config_from_file
