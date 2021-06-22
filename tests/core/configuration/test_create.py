"""
Test cwas.core.configuration.create
"""
from pathlib import Path
from random import randint

import pytest
import yaml

import cwas.core.configuration.create as create

_tmp_dir = Path.home() / f'.cwas_test_{randint(1, 1000000)}'
_annotation_dir = _tmp_dir / 'annotation-data'
_annotation_key_conf = _tmp_dir / 'annotation_key.yaml'
_gene_matrix_path = _tmp_dir / 'gene_matrix.txt'


@pytest.fixture(scope='module', autouse=True)
def create_workspace():
    _tmp_dir.mkdir()
    create_annotation_key_conf()
    create_gene_matrix()
    yield
    for f in _tmp_dir.glob('*'):
        f.unlink()
    _tmp_dir.rmdir()


@pytest.fixture(scope='module', autouse=True)
def create_annotation_dir(create_workspace):
    _annotation_dir.mkdir()
    annot_filepaths = [
        _annotation_dir / 'bed_annot1.bed.gz',
        _annotation_dir / 'bed.annot2.bed.gz',
        _annotation_dir / 'bed.annot3.bed',
        _annotation_dir / 'bw_annot1.bw',
        _annotation_dir / 'bw.annot2.bw',
        _annotation_dir / 'bw.annot3.bw.gz',
    ]
    for annot_filepath in annot_filepaths:
        annot_filepath.touch()
    yield
    for f in _annotation_dir.glob('*'):
        f.unlink()
    _annotation_dir.rmdir()


def create_annotation_key_conf():
    annot_key_dict = {
        'bed_annot1.bed.gz': 'bed1',
        'bed.annot2.bed.gz': 'bed2',
        'bed.annot3.bed': 'bed3',
        'bw_annot1.bw': 'bw1',
        'bw.annot2.bw': 'bw2',
        'bw.annot3.bw.gz': 'bw3',
    }
    with _annotation_key_conf.open('w') as out_f:
        yaml.dump(annot_key_dict, out_f)


def create_gene_matrix():
    gene_matrix_header = [
        'gene_id',
        'gene_name',
        'gene1',
        'gene2',
        'gene3',
        'gene4',
        'gene5'
    ]
    with _gene_matrix_path.open('w') as out_f:
        print(*gene_matrix_header, sep='\t', file=out_f)


def test_create_annotation_key_bed():
    bed_key_conf = _tmp_dir / 'annotation_key_bed.yaml'
    create.create_annotation_key(bed_key_conf, _annotation_dir, 'bed')
    assert bed_key_conf.exists()
    with bed_key_conf.open() as f:
        bed_key = yaml.safe_load(f)
    assert 'bed_annot1.bed.gz' in bed_key
    assert 'bed.annot2.bed.gz' in bed_key
    assert 'bed.annot.bed' not in bed_key
    assert bed_key['bed_annot1.bed.gz'] == 'bed_annot1'
    assert bed_key['bed.annot2.bed.gz'] == 'bed_annot2'


def test_create_annotation_key_bw():
    bw_key_conf = _tmp_dir / 'annotation_key_bw.yaml'
    create.create_annotation_key(bw_key_conf, _annotation_dir, 'bw')
    assert bw_key_conf.exists()
    with bw_key_conf.open() as f:
        bw_key = yaml.safe_load(f)
    assert 'bw_annot1.bw' in bw_key
    assert 'bw.annot2.bw' in bw_key
    assert 'bw.annot3.bw.gz' not in bw_key
    assert bw_key['bw_annot1.bw'] == 'bw_annot1'
    assert bw_key['bw.annot2.bw'] == 'bw_annot2'


def test_split_annotation_key():
    bed_key_conf = _tmp_dir / 'split_key_bed.yaml'
    bw_key_conf = _tmp_dir / 'split_key_bw.yaml'
    create.split_annotation_key(bed_key_conf, bw_key_conf, _annotation_key_conf)

    assert bed_key_conf.exists()
    with bed_key_conf.open() as f:
        bed_key = yaml.safe_load(f)
    assert 'bed_annot1.bed.gz' in bed_key and bed_key['bed_annot1.bed.gz'] ==\
           'bed1'
    assert 'bed.annot2.bed.gz' in bed_key and bed_key['bed.annot2.bed.gz'] == \
           'bed2'
    assert 'bed.annot3.bed' not in bed_key

    assert bw_key_conf.exists()
    with bw_key_conf.open() as f:
        bw_key = yaml.safe_load(f)
    assert 'bw_annot1.bw' in bw_key and bw_key['bw_annot1.bw'] == 'bw1'
    assert 'bw.annot2.bw' in bw_key and bw_key['bw.annot2.bw'] == 'bw2'
    assert 'bw.annot3.bw.gz' not in bw_key
