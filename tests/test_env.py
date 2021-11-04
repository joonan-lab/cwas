"""
Test cwas.env
"""
from __future__ import annotations

import os
from pathlib import Path
from typing import OrderedDict

import pytest
from cwas.env import Env


@pytest.fixture(scope="module")
def env_inst(cwas_env_path: Path) -> Env:
    env = Env(cwas_env_path)
    yield env
    if cwas_env_path.exists():
        cwas_env_path.unlink()


def test_env_init(env_inst: Env):
    assert isinstance(env_inst.env, OrderedDict)


def test_env_singleton(cwas_env_path: Path):
    # Test if the two instances are the same.
    env1 = Env(cwas_env_path)
    env2 = Env(cwas_env_path)

    assert env1 is env2

    env3 = Env(Path(str(cwas_env_path) + "_new"))
    env4 = Env()

    assert env1 is env3
    assert env1 is env4
    assert str(env1.path).endswith("_new")


def test_env_setting(env_inst: Env):
    env_key = "TEST"
    expected = "Hello!"
    env_inst.set_env(env_key, expected)
    env_value = env_inst.get_env(env_key)
    assert env_value == expected


def test_env_get_nonexist(env_inst: Env):
    env_value = env_inst.get_env("FAKE")
    assert env_value is None


def test_env_prev_set_exists(env_inst: Env):
    # Test if the environment variable exists.
    env_key = "TEST"
    expected = "Hello!"
    env_value = env_inst.get_env(env_key)
    assert env_value == expected


def test_env_reset(env_inst: Env):
    env_inst.reset()
    env_key = "TEST"
    env_value = env_inst.get_env(env_key)
    assert not env_inst.env
    assert env_value is None


def test_env_save(env_inst: Env):
    env_key = "TEST"
    expected = "Hello!"
    env_inst.set_env(env_key, expected)
    env_inst.save()

    cwas_env_path = env_inst.get_path()
    with cwas_env_path.open() as env_file:
        env_line = env_file.read()
        env_line = env_line.strip()
    assert env_line == "TEST=Hello!"

    env_inst.get_path().unlink()


def test_load_env_to_os(env_inst: Env):
    env_key = "CWAS_TEST"
    env_value = "HELLO"
    env_inst.set_env(env_key, env_value)
    env_inst.save()

    env_inst.load_env_to_os()
    assert os.getenv(env_key) == env_value

    env_inst.get_path().unlink()
    os.environ.pop(env_key)


def test_load_env_to_os_without_save(env_inst: Env):
    env_key = "CWAS_TEST"
    env_value = "HELLO"
    env_inst.set_env(env_key, env_value)
    env_inst.load_env_to_os()
    assert os.getenv(env_key) is None


def test_load_env_from_file(env_inst: Env):
    env_inst.set_env("TEST", "TEST")
    env_inst.save()
    env_inst.reset()
    assert not env_inst.env
    env_inst.load_env_from_file()
    assert env_inst.env
