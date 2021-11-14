"""
Tests of the 'Start' step
"""
from pathlib import Path

import pytest
from cwas.start import Start


@pytest.fixture(scope="module")
def args(cwas_workspace: Path) -> list:
    return ["-w", str(cwas_workspace)]


@pytest.fixture(scope="module", autouse=True)
def setup(args: list):
    inst = Start.get_instance(args)
    inst.run()


@pytest.fixture(scope="module", autouse=True)
def teardown(cwas_workspace: Path, cwas_env_path: Path):
    yield
    cwas_env_path.unlink()
    for f in cwas_workspace.glob("*"):
        f.unlink()
    cwas_workspace.rmdir()


def test_initial_file_exist(cwas_env_path: Path, cwas_workspace: Path):
    cwas_config_path = cwas_workspace / "configuration.txt"
    assert cwas_env_path.is_file()
    assert cwas_workspace.is_dir()
    assert cwas_config_path.is_file()


def test_config_keys(cwas_workspace: Path):
    config_key_set = set()
    cwas_config_path = cwas_workspace / "configuration.txt"
    with cwas_config_path.open() as config_file:
        for line in config_file:
            config_key, _ = line.strip().split("=")
            config_key_set.add(config_key)

    expected_key_set = {
        "ANNOTATION_DATA_DIR",
        "GENE_MATRIX",
        "ANNOTATION_KEY_CONFIG",
        "BIGWIG_CUTOFF_CONFIG",
        "VEP",
    }
    assert config_key_set == expected_key_set


def test_init_without_args():
    inst = Start.get_instance()
    expect_default_workspace = Path.home() / ".cwas"
    actual_workspace = getattr(inst, "workspace")
    assert expect_default_workspace == actual_workspace
