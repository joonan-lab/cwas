"""
Tests for cwas.dawn module
"""
import pytest
import argparse
import numpy as np
import pandas as pd
from pathlib import Path
from unittest.mock import patch

try:
    from cwas.dawn import Dawn
    HAS_RPY2 = True
except Exception:
    HAS_RPY2 = False

pytestmark = pytest.mark.skipif(not HAS_RPY2, reason="rpy2/R not available")


@pytest.fixture
def make_dawn():
    def _make(**overrides):
        defaults = dict(
            num_proc=2,
            eig_vector_file="/tmp/eig_vector.zarr",
            corr_mat_file="/tmp/corr_mat.zarr",
            permut_test_file="/tmp/permut_test.txt.gz",
            category_count_file="/tmp/category_count.txt",
            output_dir_path=Path("/tmp/output"),
            leiden_clustering=None,
            lambda_val=1.0,
            count_threshold=5,
            corr_threshold=0.5,
            size_threshold=10,
            k_val=None,
            k_range=(2, 10),
            seed=42,
            resolution=1.0,
            tsne_method='barnes_hut',
            tag='test',
            input_dir_path=Path("/tmp/input"),
        )
        defaults.update(overrides)
        args = argparse.Namespace(**defaults)
        with patch.object(Dawn, '_print_args'), \
             patch.object(Dawn, '_check_args_validity'), \
             patch('cwas.dawn.importr'):
            return Dawn(args)
    return _make


class TestDawnInit:
    def test_init(self, make_dawn):
        dawn = make_dawn()
        assert dawn._eig_vector is None
        assert dawn._corr_mat is None
        assert dawn._permut_test is None
        assert dawn._category_set is None
        assert dawn._k_val is None


class TestDawnProperties:
    def test_num_proc(self, make_dawn):
        dawn = make_dawn(num_proc=4)
        assert dawn.num_proc == 4

    def test_leiden_clustering_none(self, make_dawn):
        dawn = make_dawn(leiden_clustering=None)
        assert dawn.leiden_clustering is None

    def test_leiden_clustering_value(self, make_dawn):
        dawn = make_dawn(leiden_clustering='eigen_vector')
        assert dawn.leiden_clustering == 'eigen_vector'

    def test_lambda_val(self, make_dawn):
        dawn = make_dawn(lambda_val=2.5)
        assert dawn.lambda_val == 2.5

    def test_seed(self, make_dawn):
        dawn = make_dawn(seed=123)
        assert dawn.seed == 123

    def test_resolution(self, make_dawn):
        dawn = make_dawn(resolution=1.5)
        assert dawn.resolution == 1.5

    def test_tsne_method(self, make_dawn):
        dawn = make_dawn(tsne_method='barnes_hut')
        assert dawn.tsne_method == 'barnes_hut'

    def test_tag(self, make_dawn):
        dawn = make_dawn(tag='mytest')
        assert dawn.tag == 'mytest'

    def test_count_threshold(self, make_dawn):
        dawn = make_dawn(count_threshold=10)
        assert dawn.count_threshold == 10

    def test_corr_threshold(self, make_dawn):
        dawn = make_dawn(corr_threshold=0.7)
        assert dawn.corr_threshold == 0.7

    def test_size_threshold(self, make_dawn):
        dawn = make_dawn(size_threshold=20)
        assert dawn.size_threshold == 20

    def test_k_range(self, make_dawn):
        dawn = make_dawn(k_range=(3, 15))
        assert dawn.k_range == (3, 15)

    def test_output_dir_path(self, make_dawn):
        dawn = make_dawn(output_dir_path=Path("/tmp/output"))
        assert dawn.output_dir_path == Path("/tmp/output")

    def test_eig_vector_file(self, make_dawn):
        dawn = make_dawn(eig_vector_file="/tmp/eig.zarr")
        assert dawn.eig_vector_file == "/tmp/eig.zarr"

    def test_corr_mat_file(self, make_dawn):
        dawn = make_dawn(corr_mat_file="/tmp/corr.zarr")
        assert dawn.corr_mat_file == "/tmp/corr.zarr"


class TestDawnKValue:
    def test_k_val_when_provided(self, make_dawn):
        dawn = make_dawn(k_val=7)
        dawn._tsne_out = pd.DataFrame(
            np.random.rand(20, 2), columns=['t-SNE1', 't-SNE2']
        )
        assert dawn.k_val == 7
