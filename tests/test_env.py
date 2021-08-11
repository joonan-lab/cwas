"""
Test cwas.env
"""
from __future__ import annotations
from typing import OrderedDict
import pytest
import dotenv
from cwas.env import Env


class EnvMock(Env):
    def __init__(self, env_path: pathlib.Path):
        self.path = env_path
        if not self.path.exists():
            self.path.touch()
        self.env = dotenv.dotenv_values(dotenv_path=self.path)


@pytest.fixture(scope='module')
def env_path(cwas_workspace):
    return cwas_workspace / '.env'


@pytest.fixture(scope='module')
def env_mock(env_path):
    env = EnvMock(env_path)
    return env


def test_env_init(env_mock):
    assert isinstance(env_mock.env, OrderedDict)


def test_env_singleton(env_path):
    # Test if the two instances are the same.
    env1 = EnvMock(env_path)
    env2 = EnvMock(env_path)
    assert env1 is env2


def test_env_setting(env_mock):
    env_key = 'TEST'
    expected = 'Hello!'
    env_mock.set_env(env_key, expected)
    env_value = env_mock.get_env(env_key)
    assert env_value == expected


def test_env_prev_set_exists(env_mock):
    # Test if the environment variable exists.
    env_key = 'TEST'
    expected = 'Hello!'
    env_value = env_mock.get_env(env_key)
    assert env_value == expected


def test_env_reset(env_mock):
    env_mock.reset()
    env_key = 'TEST'
    env_value = env_mock.get_env(env_key)
    assert env_value is None


def test_env_save(env_mock, env_path):
    env_key = 'TEST'
    expected = 'Hello!'
    env_mock.set_env(env_key, expected)
    env_mock.save()
    with env_path.open() as env_file:
        env_line = env_file.read()
        env_line = env_line.strip()
    assert env_line == 'TEST=Hello!'
