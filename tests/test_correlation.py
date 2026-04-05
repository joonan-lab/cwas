"""
Tests for cwas.correlation module
"""
import pytest
import argparse
import pandas as pd
import numpy as np
from pathlib import Path
from unittest.mock import patch
from cwas.correlation import Correlation


def _make_corr(**kwargs):
    args = argparse.Namespace(**kwargs)
    with patch.object(Correlation, '_print_args'), \
         patch.object(Correlation, '_check_args_validity'):
        return Correlation(args)


class TestCorrelationStaticMethods:
    def test_print_args_with_variant_generation(self, capsys):
        args = argparse.Namespace(
            generate_corr_matrix='variant',
            annot_path=Path('/test/path.vcf'),
            cat_path=Path('/test/cat'),
            num_proc=4,
            generate_inter_matrix=True,
        )
        Correlation._print_args(args)
        captured = capsys.readouterr()
        assert captured.out or captured.err

    def test_print_args_with_sample_generation(self, capsys):
        args = argparse.Namespace(
            generate_corr_matrix='sample',
            cat_path=Path('/test/cat'),
            num_proc=2,
            generate_inter_matrix=None,
        )
        Correlation._print_args(args)
        captured = capsys.readouterr()
        assert captured.out or captured.err

    def test_check_args_validity_valid_args(self):
        with patch('cwas.correlation.check_num_proc'), \
             patch('cwas.correlation.check_is_file'), \
             patch('cwas.correlation.check_is_dir'):
            args = argparse.Namespace(
                generate_corr_matrix='variant',
                annot_path=Path('/test/annot.vcf'),
                cat_path=Path('/test/cat'),
                num_proc=4,
                output_dir_path=Path('/test/output'),
            )
            Correlation._check_args_validity(args)

    def test_check_args_validity_invalid_num_proc(self):
        with patch('cwas.correlation.check_num_proc', side_effect=ValueError("Invalid")):
            args = argparse.Namespace(
                generate_corr_matrix='variant',
                annot_path=Path('/test/annot.vcf'),
                cat_path=Path('/test/cat'),
                num_proc=-1,
                output_dir_path=Path('/test/output'),
            )
            with pytest.raises(ValueError):
                Correlation._check_args_validity(args)


class TestCorrelationProperties:
    def test_annot_path_property(self):
        corr = _make_corr(annot_path=Path('/test/path.vcf'))
        assert isinstance(corr.annot_path, Path)
        assert corr.annot_path.name == 'path.vcf'

    def test_num_proc_property(self):
        corr = _make_corr(num_proc=8)
        assert corr.num_proc == 8

    def test_cat_path_property(self):
        corr = _make_corr(cat_path=Path('/test/categorization'))
        assert isinstance(corr.cat_path, Path)

    def test_generate_corr_matrix_property(self):
        corr = _make_corr(generate_corr_matrix='variant')
        assert corr.generate_corr_matrix == 'variant'

    def test_generate_inter_matrix_property(self):
        corr = _make_corr(generate_inter_matrix=True)
        assert corr.generate_inter_matrix is True

    def test_output_dir_path_property(self):
        corr = _make_corr(output_dir_path=Path('/test/output'))
        assert isinstance(corr.output_dir_path, Path)


class TestCorrelationDomainList:
    @pytest.fixture
    def correlation_with_mock_category_set(self):
        corr = _make_corr(domain_list='all')
        corr._category_set = pd.DataFrame({
            'Category': ['cat1', 'cat2'],
            'is_domain1': [1, 0],
            'is_domain2': [0, 1],
        })
        return corr

    def test_domain_list_all(self, correlation_with_mock_category_set):
        assert correlation_with_mock_category_set.domain_list == ['all']

    def test_domain_list_run_all(self, correlation_with_mock_category_set):
        corr = correlation_with_mock_category_set
        corr._args.domain_list = 'run_all'
        result = corr.domain_list
        assert 'all' in result
        assert 'domain1' in result
        assert 'domain2' in result

    def test_check_domain_list_valid(self, correlation_with_mock_category_set):
        result = correlation_with_mock_category_set._check_domain_list('domain1', ['domain1', 'domain2'])
        assert result == 'domain1'

    def test_check_domain_list_exact_match(self, correlation_with_mock_category_set):
        result = correlation_with_mock_category_set._check_domain_list('domain1', ['domain1', 'domain2'])
        assert result == 'domain1'

    def test_check_domain_list_invalid(self, correlation_with_mock_category_set):
        with pytest.raises(ValueError, match="Invalid domain name"):
            correlation_with_mock_category_set._check_domain_list('invalid', ['domain1', 'domain2'])


class TestMatrixPaths:
    def test_matrix_path_with_zarr_extension(self):
        corr = _make_corr(
            cat_path=Path('/output/categorization_result.zarr'),
            output_dir_path=Path('/output'),
        )
        assert 'correlation_matrix.zarr' in str(corr.matrix_path)

    def test_matrix_path_with_gzip_extension(self):
        corr = _make_corr(
            cat_path=Path('/output/categorization_result.zarr.gz'),
            output_dir_path=Path('/output'),
        )
        assert 'correlation_matrix.zarr' in str(corr.matrix_path)

    def test_intersection_matrix_path_with_zarr_extension(self):
        corr = _make_corr(
            cat_path=Path('/output/categorization_result.zarr'),
            output_dir_path=Path('/output'),
        )
        assert 'intersection_matrix.zarr' in str(corr.intersection_matrix_path)


class TestProcessColumns:
    def test_process_columns_single_basic(self):
        matrix = pd.DataFrame({
            'col1': [1, 2, 3],
            'col2': [0, 2, 4],
            'col3': [1, 0, 3],
        })
        result = Correlation.process_columns_single(
            column_range=range(matrix.shape[1]),
            matrix=matrix
        )
        assert isinstance(result, pd.DataFrame)
        assert result.shape[1] == 3

    def test_process_columns(self):
        matrix = pd.DataFrame({
            'col1': [1, 2],
            'col2': [2, 3],
        })
        result = Correlation.process_columns(
            column_range=range(2),
            matrix=matrix
        )
        assert isinstance(result, list)
        assert len(result) == 2
        assert all(isinstance(item, pd.Series) for item in result)
