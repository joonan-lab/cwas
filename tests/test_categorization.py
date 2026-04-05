"""
Tests for cwas.categorization module (top-level wrapper)
"""
import pytest
import argparse
from pathlib import Path
from unittest.mock import Mock, patch
from cwas.categorization import Categorization


class TestCategorizationProperties:
    @pytest.fixture
    def categorization(self):
        args = Mock(spec=argparse.Namespace)
        args.input_path = Path('/test/annotated.vcf')
        args.output_dir_path = Path('/test/output')
        args.num_proc = 1
        with patch.object(Categorization, '_print_args'), \
             patch.object(Categorization, '_check_args_validity'):
            return Categorization(args)

    def test_num_proc_property(self, categorization):
        assert categorization.num_proc == 1

    def test_input_path_property(self, categorization):
        assert isinstance(categorization.input_path, Path)

    def test_output_dir_path_property(self, categorization):
        assert isinstance(categorization.output_dir_path, Path)


class TestResultPathGeneration:
    def test_result_path_conversion_vcf(self):
        args = Mock(spec=argparse.Namespace)
        args.input_path = Path('/test/annotated.vcf')
        args.output_dir_path = Path('/test/output')
        args.num_proc = 1
        with patch.object(Categorization, '_print_args'), \
             patch.object(Categorization, '_check_args_validity'):
            cat = Categorization(args)
        result = cat.result_path
        assert result.suffix == '.zarr'
        assert 'categorization_result' in result.name

    def test_result_path_conversion_vcf_gz(self):
        args = Mock(spec=argparse.Namespace)
        args.input_path = Path('/test/annotated.vcf.gz')
        args.output_dir_path = Path('/test/output')
        args.num_proc = 1
        with patch.object(Categorization, '_print_args'), \
             patch.object(Categorization, '_check_args_validity'):
            cat = Categorization(args)
        result = cat.result_path
        assert 'categorization_result.zarr' in result.name


class TestInitialization:
    def test_initialization_sets_attributes(self):
        args = Mock(spec=argparse.Namespace)
        args.input_path = Path('/test/annotated.vcf')
        args.output_dir_path = Path('/test/output')
        args.num_proc = 1
        with patch.object(Categorization, '_print_args'), \
             patch.object(Categorization, '_check_args_validity'):
            cat = Categorization(args)
        assert cat._annotated_vcf is None
        assert cat._gene_matrix is None
        assert cat._category_domain is None
        assert cat._redundant_categories is None
        assert cat._sample_ids is None
        assert cat._categories is None
        assert cat._result is None


class TestStaticMethods:
    def test_check_args_validity(self):
        args = Mock(spec=argparse.Namespace)
        args.input_path = Path('/valid/path')
        args.output_dir_path = Path('/valid/path')
        args.num_proc = 1
        with patch('cwas.categorization.check_num_proc'), \
             patch('cwas.categorization.check_is_file'), \
             patch('cwas.categorization.check_is_dir'):
            Categorization._check_args_validity(args)
