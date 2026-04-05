"""
Tests for cwas.permutation_test module
"""
import pytest
import argparse
import pandas as pd
import numpy as np
from pathlib import Path
from unittest.mock import patch
from cwas.permutation_test import PermutationTest


def _make_pt(**kwargs):
    args = argparse.Namespace(**kwargs)
    with patch.object(PermutationTest, '_print_args'), \
         patch.object(PermutationTest, '_check_args_validity'):
        return PermutationTest(args)


class TestPermutationTestProperties:
    def test_cat_path_property(self):
        pt = _make_pt(cat_path=Path('/test/categorization.zarr'))
        assert isinstance(pt.cat_path, Path)

    def test_output_dir_path_property(self):
        pt = _make_pt(output_dir_path=Path('/test/output'))
        assert isinstance(pt.output_dir_path, Path)

    def test_burden_shift_property(self):
        pt = _make_pt(burden_shift=True)
        assert pt.burden_shift is True

    def test_burden_shift_property_false(self):
        pt = _make_pt(burden_shift=False)
        assert pt.burden_shift is False

    def test_use_n_carrier_property(self):
        pt = _make_pt(use_n_carrier=True)
        assert pt.use_n_carrier is True


class TestResultPaths:
    def test_result_path_basic(self):
        pt = _make_pt(
            cat_path=Path('/output/categorization_result.zarr'),
            output_dir_path=Path('/output'),
        )
        assert 'permutation_test.txt.gz' in str(pt.result_path)

    def test_perm_rrs_path(self):
        pt = _make_pt(
            cat_path=Path('/output/categorization_result.zarr'),
            output_dir_path=Path('/output'),
        )
        assert 'permutation_RRs.txt.gz' in str(pt.perm_rrs_path)

    def test_binom_pvals_path(self):
        pt = _make_pt(
            cat_path=Path('/output/categorization_result.zarr'),
            output_dir_path=Path('/output'),
        )
        assert 'binom_pvals.txt.gz' in str(pt.binom_pvals_path)

    def test_result_path_with_gzip(self):
        pt = _make_pt(
            cat_path=Path('/output/categorization_result.zarr.gz'),
            output_dir_path=Path('/output'),
        )
        assert 'permutation_test.txt.gz' in str(pt.result_path)


class TestGetPermPval:
    def test_get_perm_pval_high_rr(self):
        pt = _make_pt()
        perm_rrs = np.array([
            [0.5, 0.6, 0.7],
            [0.8, 0.9, 1.0],
            [1.1, 1.2, 1.3],
        ])
        rr = np.array([1.2, 1.0, 0.8])
        result = pt.get_perm_pval(perm_rrs, rr)
        assert all(0.0 <= v <= 1.0 for v in result)

    def test_get_perm_pval_exact_match(self):
        pt = _make_pt()
        perm_rrs = np.array([
            [1.0, 1.0],
            [1.0, 1.0],
            [1.0, 1.0],
        ])
        rr = np.array([1.0, 1.0])
        result = pt.get_perm_pval(perm_rrs, rr)
        assert np.allclose(result, 1.0)


class TestBurdenTestReproducibility:
    def test_burden_test_same_seed_same_result(self):
        case_cnt = 3
        ctrl_cnt = 3
        var_counts = np.array([
            [1, 0],
            [2, 1],
            [0, 1],
            [1, 1],
            [0, 2],
            [1, 0],
        ])
        result1 = PermutationTest._burden_test(
            seed_range=(10001, 10002),
            case_cnt=case_cnt, ctrl_cnt=ctrl_cnt,
            var_counts=var_counts,
            use_n_carrier=False, burden_shift=False,
        )
        result2 = PermutationTest._burden_test(
            seed_range=(10001, 10002),
            case_cnt=case_cnt, ctrl_cnt=ctrl_cnt,
            var_counts=var_counts,
            use_n_carrier=False, burden_shift=False,
        )
        np.testing.assert_array_almost_equal(result1[0], result2[0])
