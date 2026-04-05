"""
Tests for cwas.burden_test module
"""
import pytest
import argparse
from pathlib import Path
from unittest.mock import Mock, patch
import pandas as pd
import numpy as np
from cwas.burden_test import BurdenTest, _contain_same_index, apply_region_mapping


class BurdenTestMock(BurdenTest):
    def run_burden_test(self):
        pass

    def run(self):
        pass


class TestContainSameIndex:
    def test_same_index(self):
        df1 = pd.DataFrame({'a': [1, 2, 3]}, index=['x', 'y', 'z'])
        df2 = pd.DataFrame({'b': [4, 5, 6]}, index=['x', 'y', 'z'])
        assert _contain_same_index(df1, df2) is True

    def test_different_index_order(self):
        df1 = pd.DataFrame({'a': [1, 2, 3]}, index=['x', 'y', 'z'])
        df2 = pd.DataFrame({'b': [4, 5, 6]}, index=['z', 'y', 'x'])
        assert _contain_same_index(df1, df2) is True

    def test_different_index_values(self):
        df1 = pd.DataFrame({'a': [1, 2, 3]}, index=['x', 'y', 'z'])
        df2 = pd.DataFrame({'b': [4, 5, 6]}, index=['a', 'b', 'c'])
        assert _contain_same_index(df1, df2) is False

    def test_subset_index(self):
        df1 = pd.DataFrame({'a': [1, 2, 3]}, index=['x', 'y', 'z'])
        df2 = pd.DataFrame({'b': [4, 5]}, index=['x', 'y'])
        assert _contain_same_index(df1, df2) is False


class TestApplyRegionMapping:
    def test_creates_binary_columns(self):
        df = pd.DataFrame({
            'gencode': ['CodingRegion', 'NoncodingRegion', 'PTVRegion'],
            'gene_set': ['geneset1', 'geneset2', 'geneset3'],
            'functional_score': ['score1', 'score2', 'score3'],
            'functional_annotation': ['annot1', 'annot2', 'annot3']
        })
        result = apply_region_mapping(df)
        assert 'is_coding' in result.columns
        assert 'is_noncoding' in result.columns
        assert 'is_PTV' in result.columns

    def test_coding_detection(self):
        df = pd.DataFrame({
            'gencode': ['CodingRegion', 'NoncodingRegion'],
            'gene_set': ['gs1', 'gs1'],
            'functional_score': ['fs1', 'fs1'],
            'functional_annotation': ['fa1', 'fa1']
        })
        result = apply_region_mapping(df)
        assert result['is_coding'].iloc[0] == 1
        assert result['is_coding'].iloc[1] == 0

    def test_ptv_detection(self):
        df = pd.DataFrame({
            'gencode': ['PTVRegion', 'MissenseRegion'],
            'gene_set': ['gs1', 'gs1'],
            'functional_score': ['fs1', 'fs1'],
            'functional_annotation': ['fa1', 'fa1']
        })
        result = apply_region_mapping(df)
        assert result['is_PTV'].iloc[0] == 1
        assert result['is_PTV'].iloc[1] == 0

    def test_noncoding_detection(self):
        df = pd.DataFrame({
            'gencode': ['IntronRegion', 'CodingRegion'],
            'gene_set': ['gs1', 'gs1'],
            'functional_score': ['fs1', 'fs1'],
            'functional_annotation': ['fa1', 'fa1']
        })
        result = apply_region_mapping(df)
        assert result['is_noncoding'].iloc[0] == 1
        assert result['is_noncoding'].iloc[1] == 0


class TestBurdenTestProperties:
    @pytest.fixture
    def burden_test(self):
        args = Mock(spec=argparse.Namespace)
        args.cat_path = Path('/test/cat_path')
        args.output_dir_path = Path('/test/output')
        args.sample_info_path = Path('/test/sample_info.txt')
        args.adj_factor_path = None
        args.use_n_carrier = False
        args.tag = 'test_tag'
        args.eff_test = 20
        args.marker_size = 15.0
        args.font_size = 12.0
        args.plot_size = 7.0
        args.plot_title = 'Test Plot'
        with patch.object(BurdenTestMock, '_print_args'), \
             patch.object(BurdenTestMock, '_check_args_validity'):
            return BurdenTestMock(args)

    def test_cat_path_property(self, burden_test):
        assert isinstance(burden_test.cat_path, Path)

    def test_output_dir_path_property(self, burden_test):
        assert isinstance(burden_test.output_dir_path, Path)

    def test_use_n_carrier_property(self, burden_test):
        assert burden_test.use_n_carrier is False

    def test_tag_property(self, burden_test):
        assert burden_test.tag == 'test_tag'

    def test_eff_test_property(self, burden_test):
        assert burden_test.eff_test == 20

    def test_adj_factor_path_none(self, burden_test):
        assert burden_test.adj_factor_path is None

    def test_result_path_conversion(self):
        args = Mock(spec=argparse.Namespace)
        args.cat_path = Path('/test/output/categorization_result.zarr')
        args.output_dir_path = Path('/test/output')
        args.sample_info_path = Path('/test/sample_info.txt')
        args.adj_factor_path = None
        args.use_n_carrier = False
        args.tag = None
        args.eff_test = 1
        args.marker_size = 15.0
        args.font_size = 12.0
        args.plot_size = 7.0
        args.plot_title = 'Test'
        with patch.object(BurdenTestMock, '_print_args'), \
             patch.object(BurdenTestMock, '_check_args_validity'):
            bt = BurdenTestMock(args)
        assert bt.result_path.name.endswith('burden_test.txt')

    def test_counts_path_conversion(self):
        args = Mock(spec=argparse.Namespace)
        args.cat_path = Path('/test/output/categorization_result.zarr')
        args.output_dir_path = Path('/test/output')
        args.sample_info_path = Path('/test/sample_info.txt')
        args.adj_factor_path = None
        args.use_n_carrier = False
        args.tag = None
        args.eff_test = 1
        args.marker_size = 15.0
        args.font_size = 12.0
        args.plot_size = 7.0
        args.plot_title = 'Test'
        with patch.object(BurdenTestMock, '_print_args'), \
             patch.object(BurdenTestMock, '_check_args_validity'):
            bt = BurdenTestMock(args)
        assert bt.counts_path.name.endswith('category_counts.txt')
