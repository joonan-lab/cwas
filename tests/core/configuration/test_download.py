"""
Test cwas.core.configuration.download
"""
from pathlib import Path
from random import randint

import pytest

import cwas.core.configuration.download as download

_tmp_dir = Path.home() / f'.cwas_test_{randint(1, 1000000)}'


@pytest.fixture(scope='session', autouse=True)
def run_around_tests():
    _tmp_dir.mkdir()
    yield
    for f in _tmp_dir.glob('*'):
        f.unlink()
    _tmp_dir.rmdir()


def test_download_bed_tar_gz():
    download.download_bed_tar_gz(str(_tmp_dir))
    expect_filepath = _tmp_dir / 'bed.tar.gz'
    assert expect_filepath.exists()


def test_download_bw_tar_gz():
    download.download_bw_tar_gz(str(_tmp_dir))
    expect_filepath = _tmp_dir / 'bw.tar.gz'
    assert expect_filepath.exists()
