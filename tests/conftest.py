"""
CWAS test fixtures
"""
from pathlib import Path
from random import randint

import pytest
import yaml
from cwas.env import Env

_rand_n = randint(1, 1000000)


@pytest.fixture(scope="package")
def cwas_env_path():
    _tmp_env = Path.home() / f".cwas_env_{_rand_n}"
    return _tmp_env


@pytest.fixture(scope="package", autouse=True)
def cwas_env_inst(cwas_env_path):
    env_inst = Env(cwas_env_path)
    yield
    env_inst.remove_file()


@pytest.fixture(scope="package")
def cwas_workspace():
    _tmp_workspace = Path.home() / f".cwas-test-{_rand_n}"
    return _tmp_workspace


@pytest.fixture(scope="package")
def cwas_input_dir():
    _tmp_input_dir = Path.home() / f".cwas-test-input-{_rand_n}"
    return _tmp_input_dir


@pytest.fixture(scope="package")
def annotation_dir(cwas_input_dir):
    _annotation_dir = cwas_input_dir / "annotation-data"
    return _annotation_dir


@pytest.fixture(scope="package")
def annotation_key_conf(cwas_input_dir):
    _annotation_key_conf = cwas_input_dir / "annotation_key.yaml"
    return _annotation_key_conf


@pytest.fixture(scope="package")
def bw_cutoff_conf(cwas_input_dir):
    _bw_cutoff_conf = cwas_input_dir / "user_def_bigwig_cutoff.yaml"
    return _bw_cutoff_conf


@pytest.fixture(scope="package")
def gene_matrix(cwas_input_dir):
    _gene_matrix = cwas_input_dir / "gene_matrix.txt"
    return _gene_matrix


@pytest.fixture(scope="package", autouse=True)
def create_cwas_input_dir(
    cwas_input_dir, gene_matrix, annotation_key_conf, bw_cutoff_conf
):
    cwas_input_dir.mkdir()
    create_gene_matrix(gene_matrix)
    create_annotation_key_conf(annotation_key_conf)
    create_bw_cutoff_conf(bw_cutoff_conf)
    print("[TEST] Temporary CWAS input directory has created.")
    yield
    for f in cwas_input_dir.glob("*"):
        f.unlink()
    cwas_input_dir.rmdir()
    print("[TEST] Temporary CWAS input directory has deleted.")


@pytest.fixture(scope="package", autouse=True)
def create_annotation_dir(create_cwas_input_dir, annotation_dir):
    annotation_dir.mkdir()
    annot_filepaths = [
        annotation_dir / "bed_annot1.bed.gz",
        annotation_dir / "bed.annot2.bed.gz",
        annotation_dir / "bed.annot3.bed",
        annotation_dir / "bw_annot1.bw",
        annotation_dir / "bw.annot2.bw",
        annotation_dir / "bw.annot3.bw.gz",
    ]
    for annot_filepath in annot_filepaths:
        annot_filepath.touch()
    yield
    for f in annotation_dir.glob("*"):
        f.unlink()
    annotation_dir.rmdir()


def create_gene_matrix(gene_matrix):
    gene_matrix_header = [
        "gene_id",
        "gene_name",
        "gene1",
        "gene2",
        "gene3",
        "gene4",
        "gene5",
    ]
    with gene_matrix.open("w") as out_f:
        print(*gene_matrix_header, sep="\t", file=out_f)


def create_annotation_key_conf(annotation_key_conf):
    annot_key_dict = {
        "bed_annot1.bed.gz": "bed1",
        "bed.annot2.bed.gz": "bed2",
        "bed.annot3.bed": "bed3",
        "bw_annot1.bw": "bw1",
        "bw.annot2.bw": "bw2",
        "bw.annot3.bw.gz": "bw3",
    }
    with annotation_key_conf.open("w") as out_f:
        yaml.safe_dump(annot_key_dict, out_f)


def create_bw_cutoff_conf(bw_cutoff_conf):
    bw_cutoff_dict = {
        "bw_annot1.bw": 1.0,
        "bw.annot2.bw": -2.2,
        "bw.annot3.bw.gz": 3.3,
    }
    with bw_cutoff_conf.open("w") as out_f:
        yaml.safe_dump(bw_cutoff_dict, out_f)
