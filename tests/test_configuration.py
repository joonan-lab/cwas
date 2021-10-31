"""
Test cwas.configuration
"""
import os
from pathlib import Path

import pytest
from cwas.configuration import Configuration


@pytest.fixture(scope="module")
def set_cwas_env(cwas_workspace: Path):
    cwas_env_path = Path.home() / ".cwas_env"
    with cwas_env_path.open("w") as cwas_env_file:
        print(f"CWAS_WORKSPACE={str(cwas_workspace)}", file=cwas_env_file)
    os.environ["CWAS_WORKSPACE"] = str(cwas_workspace)


@pytest.fixture(scope="module")
def create_cwas_conf(
    cwas_workspace,
    annotation_dir,
    annotation_key_conf,
    bw_cutoff_conf,
    gene_matrix,
    vep,
):
    config = {
        "ANNOTATION_DATA_DIR": annotation_dir,
        "GENE_MATRIX": gene_matrix,
        "ANNOTATION_KEY_CONFIG": annotation_key_conf,
        "BIGWIG_CUTOFF_CONFIG": bw_cutoff_conf,
        "VEP": vep,
    }
    config_path = cwas_workspace / "configuration.txt"

    with config_path.open("w") as config_file:
        for k, v in config.items():
            print(f"{k}={str(v)}", file=config_file)


@pytest.fixture(scope="module")
def configuration_inst(set_cwas_env, create_cwas_conf):
    inst = Configuration.get_instance()
    return inst


def test_run_configuration_make_files(cwas_workspace, configuration_inst):
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
