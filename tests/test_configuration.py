"""
Tests of the 'Configuration' step
"""
import os
import random
from pathlib import Path

import pytest
from cwas.configuration import Configuration


@pytest.fixture(scope="module", autouse=True)
def create_cwas_env_file(cwas_workspace: Path):
    cwas_env_path = Path.home() / ".cwas_env"
    with cwas_env_path.open("w") as cwas_env_file:
        print(f"CWAS_WORKSPACE={str(cwas_workspace)}", file=cwas_env_file)
    yield
    cwas_env_path.unlink()


@pytest.fixture(scope="module")
def load_env_to_os(cwas_workspace: Path):
    os.environ["CWAS_WORKSPACE"] = str(cwas_workspace)
    yield
    os.unsetenv("CWAS_WORKSPACE")


@pytest.fixture
def cwas_config(
    annotation_dir: Path,
    annotation_key_conf: Path,
    bw_cutoff_conf: Path,
    gene_matrix: Path,
):
    config = {
        "ANNOTATION_DATA_DIR": annotation_dir,
        "GENE_MATRIX": gene_matrix,
        "ANNOTATION_KEY_CONFIG": annotation_key_conf,
        "BIGWIG_CUTOFF_CONFIG": bw_cutoff_conf,
        "VEP": "VEP",
    }
    return config


@pytest.fixture
def create_cwas_config_file(cwas_workspace, cwas_config):
    _create_cwas_config_file(cwas_workspace, cwas_config)


@pytest.fixture
def create_cwas_config_file_incomplete(cwas_workspace, cwas_config):
    random_config_key = random.choice(list(cwas_config.keys()))
    cwas_config[random_config_key] = ""
    _create_cwas_config_file(cwas_workspace, cwas_config)


def _create_cwas_config_file(cwas_workspace, cwas_config):
    config_path = cwas_workspace / "configuration.txt"

    with config_path.open("w") as config_file:
        for k, v in cwas_config.items():
            print(f"{k}={str(v)}", file=config_file)


def test_get_inst_without_load_to_env():
    with pytest.raises(RuntimeError):
        Configuration.get_instance()


@pytest.fixture
def configuration_inst(load_env_to_os):
    return Configuration.get_instance()


def test_run_configuration_with_incomplete_config(
    configuration_inst, create_cwas_config_file_incomplete,
):
    with pytest.raises(ValueError):
        configuration_inst.run()


def test_run_configuration_make_files(
    cwas_workspace, configuration_inst, create_cwas_config_file
):
    configuration_inst.run()

    data_dir_symlink = cwas_workspace / "annotation-data"
    gene_matrix_symlink = cwas_workspace / "gene_matrix.txt"
    bed_key_list = cwas_workspace / "annotation_key_bed.yaml"
    bw_key_list = cwas_workspace / "annotation_key_bw.yaml"
    bw_cutoff_list = cwas_workspace / "annotation_cutoff_bw.yaml"
    category_domain_list = cwas_workspace / "category_domain.yaml"
    redundant_category_table = cwas_workspace / "redundant_category.txt"

    assert data_dir_symlink.is_dir() and data_dir_symlink.is_symlink()
    assert gene_matrix_symlink.is_file() and gene_matrix_symlink.is_symlink()
    assert bed_key_list.is_file()
    assert bw_key_list.is_file()
    assert bw_cutoff_list.is_file()
    assert category_domain_list.is_file()
    assert redundant_category_table.is_file()

    # Teardown
    data_dir_symlink.unlink()
    gene_matrix_symlink.unlink()
    bed_key_list.unlink()
    bw_key_list.unlink()
    bw_cutoff_list.unlink()
    category_domain_list.unlink()
    redundant_category_table.unlink()


def test_env_after_run_configuration(configuration_inst):
    configuration_inst.run()

    env_keys = [
        "VEP",
        "CWAS_WORKSPACE",
        "ANNOTATION_DATA",
        "GENE_MATRIX",
        "ANNOTATION_BED_KEY",
        "ANNOTATION_BW_KEY",
        "ANNOTATION_BW_CUTOFF",
        "CATEGORY_DOMAIN",
        "REDUNDANT_CATEGORY",
    ]
    for env_key in env_keys:
        assert configuration_inst.get_env(env_key)
