"""
Tests for cwas.effective_num_test module
"""
import pytest
import argparse
import pandas as pd
import numpy as np
from pathlib import Path
from unittest.mock import patch
from cwas.effective_num_test import EffectiveNumTest


def _make_ent(**kwargs):
    args = argparse.Namespace(**kwargs)
    with patch.object(EffectiveNumTest, '_print_args'), \
         patch.object(EffectiveNumTest, '_check_args_validity'):
        return EffectiveNumTest(args)


class TestEffectiveNumTestStaticMethods:
    def test_print_args_minimal(self, capsys):
        args = argparse.Namespace(
            input_path=Path('/test/input'),
            input_format='corr',
            output_dir_path=Path('/test/output'),
            sample_info_path=None,
            category_count_file=None,
            domain_list=None,
            tag=None,
            category_set_path=None,
        )
        EffectiveNumTest._print_args(args)
        captured = capsys.readouterr()
        assert captured.out or captured.err

    def test_check_args_validity_valid(self):
        with patch('cwas.effective_num_test.check_is_dir'), \
             patch('cwas.effective_num_test.check_is_file'):
            args = argparse.Namespace(
                input_path=Path('/test/input'),
                output_dir_path=Path('/test/output'),
                sample_info_path=Path('/test/sample.txt'),
                category_count_file=None,
                category_set_path=None,
                count_thres=10,
            )
            EffectiveNumTest._check_args_validity(args)

    def test_check_args_validity_missing_sample_and_threshold(self):
        with patch('cwas.effective_num_test.check_is_dir'):
            args = argparse.Namespace(
                input_path=Path('/test/input'),
                output_dir_path=Path('/test/output'),
                sample_info_path=None,
                category_count_file=None,
                category_set_path=None,
                count_thres=None,
            )
            with pytest.raises(Exception):
                EffectiveNumTest._check_args_validity(args)


class TestEffectiveNumTestProperties:
    def test_input_path_property(self):
        ent = _make_ent(input_path=Path('/test/input'))
        assert isinstance(ent.input_path, Path)

    def test_output_dir_path_property(self):
        ent = _make_ent(output_dir_path=Path('/test/output'))
        assert isinstance(ent.output_dir_path, Path)

    def test_input_format_property(self):
        ent = _make_ent(input_format='corr')
        assert ent.input_format == 'corr'

    def test_num_eig_property(self):
        ent = _make_ent(num_eig=100)
        assert ent.num_eig == 100

    def test_tag_property(self):
        ent = _make_ent(tag='mytag')
        assert ent.tag == 'mytag'

    def test_sample_info_path_property(self):
        ent = _make_ent(sample_info_path=Path('/test/sample.txt'))
        assert isinstance(ent.sample_info_path, Path)

    def test_count_thres_explicit_value(self):
        ent = _make_ent(count_thres=15, sample_info_path=None)
        assert ent.count_thres == 15


class TestPathProperties:
    def test_neg_lap_path_default_tag(self):
        ent = _make_ent(
            input_path=Path('/output/data.intersection_matrix.zarr'),
            output_dir_path=Path('/output'),
            tag=None,
        )
        result = ent.neg_lap_path
        assert 'neg_lap' in str(result)
        assert str(result).endswith('.zarr')

    def test_neg_lap_path_with_tag(self):
        ent = _make_ent(
            input_path=Path('/output/data.intersection_matrix.zarr'),
            output_dir_path=Path('/output'),
            tag='mytag',
        )
        assert 'mytag' in str(ent.neg_lap_path)

    def test_eig_val_path(self):
        ent = _make_ent(
            input_path=Path('/output/data.correlation_matrix.zarr'),
            output_dir_path=Path('/output'),
            tag=None,
        )
        assert 'eig_vals' in str(ent.eig_val_path)

    def test_eig_vec_path(self):
        ent = _make_ent(
            input_path=Path('/output/data.correlation_matrix.zarr'),
            output_dir_path=Path('/output'),
            tag=None,
        )
        assert 'eig_vecs' in str(ent.eig_vec_path)


class TestCheckDomainList:
    def test_check_domain_list_valid(self):
        ent = _make_ent()
        assert ent._check_domain_list('domain1', ['domain1', 'domain2']) == 'domain1'

    def test_check_domain_list_exact_match(self):
        ent = _make_ent()
        assert ent._check_domain_list('domain1', ['domain1', 'domain2']) == 'domain1'

    def test_check_domain_list_invalid(self):
        ent = _make_ent()
        with pytest.raises(ValueError, match="Invalid domain name"):
            ent._check_domain_list('invalid_domain', ['domain1', 'domain2'])
