"""
Test cwas.configuration
"""
import pytest
from cwas.configuration import Configuration


@pytest.fixture(scope="module")
def configuration_inst(
    cwas_workspace,
    annotation_dir,
    annotation_key_conf,
    bw_cutoff_conf,
    gene_matrix,
):
    argv = [
        "-d",
        str(annotation_dir),
        "-m",
        str(gene_matrix),
        "-k",
        str(annotation_key_conf),
        "-c",
        str(bw_cutoff_conf),
        "-w",
        str(cwas_workspace),
    ]
    inst = Configuration.get_instance(argv)
    return inst


def test_run_configuration_make_files(cwas_workspace, configuration_inst):
    env_file_path = cwas_workspace / ".config"
    configuration_inst.set_env_path(env_file_path)
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
    assert env_file_path.is_file()

    # Teardown
    data_dir_symlink.unlink()
    gene_matrix_symlink.unlink()
    bed_key_list.unlink()
    bw_key_list.unlink()
    bw_cutoff_list.unlink()
    category_domain_list.unlink()
    redundant_category_table.unlink()
    env_file_path.unlink()


def test_env_after_run_configuration(cwas_workspace, configuration_inst):
    env_file_path = cwas_workspace / ".config"
    configuration_inst.set_env_path(env_file_path)
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
