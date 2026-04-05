"""
Tests for cwas.burden_shift module
"""
import pytest
import argparse
from pathlib import Path
from unittest.mock import Mock, patch
from cwas.burden_shift import BurdenShift


class TestBurdenShiftProperties:
    @pytest.fixture
    def burden_shift(self):
        args = Mock(spec=argparse.Namespace)
        args.input_path = Path('/test/burden_test.txt')
        args.burden_res = Path('/test/burden_res.txt')
        args.output_dir_path = Path('/test/output')
        args.cat_set_file = Path('/test/cat_set.txt')
        args.cat_count_file = Path('/test/cat_count.txt')
        args.count_cutoff = 7
        args.pval = 0.05
        args.tag = 'test_tag'
        args.cat_set_list = None
        args.n_cat_sets = 10
        args.fontsize = 10
        args.plot_title = 'Test Plot'
        with patch.object(BurdenShift, '_print_args'), \
             patch.object(BurdenShift, '_check_args_validity'):
            return BurdenShift(args)

    def test_input_file_property(self, burden_shift):
        assert isinstance(burden_shift.input_file, Path)

    def test_plot_title_property(self, burden_shift):
        assert burden_shift.plot_title == 'Test Plot'

    def test_output_dir_path_property(self, burden_shift):
        assert isinstance(burden_shift.output_dir_path, Path)

    def test_tag_property(self, burden_shift):
        assert burden_shift.tag == 'test_tag'

    def test_pval_property(self, burden_shift):
        assert burden_shift.pval == 0.05

    def test_n_cat_sets_property(self, burden_shift):
        assert burden_shift.n_cat_sets == 10

    def test_cat_set_list_none(self, burden_shift):
        assert burden_shift.cat_set_list is None


class TestCountCutoffValidation:
    @pytest.fixture
    def make_burden_shift(self):
        def _make(cutoff):
            args = Mock(spec=argparse.Namespace)
            args.input_path = Path('/test/a.txt')
            args.burden_res = Path('/test/b.txt')
            args.output_dir_path = Path('/test/output')
            args.cat_set_file = Path('/test/c.txt')
            args.cat_count_file = Path('/test/d.txt')
            args.count_cutoff = cutoff
            args.pval = 0.05
            args.tag = None
            args.cat_set_list = None
            args.n_cat_sets = 10
            args.fontsize = 10
            args.plot_title = 'Test'
            with patch.object(BurdenShift, '_print_args'), \
                 patch.object(BurdenShift, '_check_args_validity'):
                return BurdenShift(args)
        return _make

    def test_positive_count_cutoff(self, make_burden_shift):
        bs = make_burden_shift(7)
        assert bs.c_cutoff == 7

    def test_zero_count_cutoff_accepted(self, make_burden_shift):
        bs = make_burden_shift(0)
        assert bs.c_cutoff == 0

    def test_negative_count_cutoff_raises(self, make_burden_shift):
        bs = make_burden_shift(-1)
        with pytest.raises(ValueError):
            _ = bs.c_cutoff


class TestCountCats:
    @pytest.fixture
    def burden_shift(self):
        args = Mock(spec=argparse.Namespace)
        args.input_path = Path('/test/a.txt')
        args.burden_res = Path('/test/b.txt')
        args.output_dir_path = Path('/test/output')
        args.cat_set_file = Path('/test/c.txt')
        args.cat_count_file = Path('/test/d.txt')
        args.count_cutoff = 7
        args.pval = 0.05
        args.tag = None
        args.cat_set_list = None
        args.n_cat_sets = 10
        args.fontsize = 10
        args.plot_title = 'Test'
        with patch.object(BurdenShift, '_print_args'), \
             patch.object(BurdenShift, '_check_args_validity'):
            return BurdenShift(args)

    def test_count_cats_positive_values(self, burden_shift):
        pvals = [0.01, 0.02, 0.03, 0.04, 0.05, 0.06]
        nCase, nCtrl = burden_shift._count_cats(pvals, 0.05)
        assert nCase == 5
        assert nCtrl == 0

    def test_count_cats_negative_values(self, burden_shift):
        pvals = [-0.01, -0.02, -0.03, -0.04, -0.05, -0.06]
        nCase, nCtrl = burden_shift._count_cats(pvals, 0.05)
        assert nCase == 0
        assert nCtrl == 5

    def test_count_cats_mixed_values(self, burden_shift):
        pvals = [0.01, 0.02, -0.01, -0.02, 0.06, -0.06]
        nCase, nCtrl = burden_shift._count_cats(pvals, 0.05)
        assert nCase == 2
        assert nCtrl == 2


class TestBurdenShiftSize:
    @pytest.fixture
    def burden_shift(self):
        args = Mock(spec=argparse.Namespace)
        args.input_path = Path('/test/a.txt')
        args.burden_res = Path('/test/b.txt')
        args.output_dir_path = Path('/test/output')
        args.cat_set_file = Path('/test/c.txt')
        args.cat_count_file = Path('/test/d.txt')
        args.count_cutoff = 7
        args.pval = 0.05
        args.tag = None
        args.cat_set_list = None
        args.n_cat_sets = 10
        args.fontsize = 10
        args.plot_title = 'Test'
        with patch.object(BurdenShift, '_print_args'), \
             patch.object(BurdenShift, '_check_args_validity'):
            return BurdenShift(args)

    def test_size_below_first_bin(self, burden_shift):
        bins = [10, 20, 30, 40, 50, 60]
        assert burden_shift._burden_shift_size(5, bins) == 1

    def test_size_in_first_range(self, burden_shift):
        bins = [10, 20, 30, 40, 50, 60]
        assert burden_shift._burden_shift_size(15, bins) == 3

    def test_size_above_last_bin(self, burden_shift):
        bins = [10, 20, 30, 40, 50, 60]
        assert burden_shift._burden_shift_size(70, bins) == 13
