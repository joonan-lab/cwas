"""
Tests for cwas.extract_variant module
"""
import pytest
import argparse
import pandas as pd
import numpy as np
from pathlib import Path
from unittest.mock import patch
from cwas.extract_variant import ExtractVariant


def _make_ev(**kwargs):
    args = argparse.Namespace(**kwargs)
    with patch.object(ExtractVariant, '_print_args'), \
         patch.object(ExtractVariant, '_check_args_validity'):
        return ExtractVariant(args)


class TestExtractVariantStaticMethods:
    def test_check_args_validity_valid_args(self):
        with patch('cwas.extract_variant.check_is_file'), \
             patch('cwas.extract_variant.check_is_dir'):
            args = argparse.Namespace(
                input_path=Path('/test/input.vcf'),
                output_dir_path=Path('/test/output'),
                category_set_path=None,
            )
            ExtractVariant._check_args_validity(args)

    def test_check_args_validity_with_category_set(self):
        with patch('cwas.extract_variant.check_is_file'), \
             patch('cwas.extract_variant.check_is_dir'):
            args = argparse.Namespace(
                input_path=Path('/test/input.vcf'),
                output_dir_path=Path('/test/output'),
                category_set_path=Path('/test/categories.txt'),
            )
            ExtractVariant._check_args_validity(args)


class TestExtractVariantProperties:
    def test_input_path_property(self):
        ev = _make_ev(input_path=Path('/test/input.vcf'))
        assert isinstance(ev.input_path, Path)

    def test_output_dir_path_property(self):
        ev = _make_ev(output_dir_path=Path('/test/output'))
        assert isinstance(ev.output_dir_path, Path)

    def test_annotation_info_property(self):
        ev = _make_ev(annotation_info=True)
        assert ev.annotation_info is True

    def test_tag_property(self):
        ev = _make_ev(tag='mytag')
        assert ev.tag == 'mytag'

    def test_tag_property_none(self):
        ev = _make_ev(tag=None)
        assert ev.tag is None

    def test_category_set_path_property(self):
        ev = _make_ev(category_set_path=Path('/test/categories.txt'))
        assert isinstance(ev.category_set_path, Path)

    def test_category_set_path_property_none(self):
        ev = _make_ev(category_set_path=None)
        assert ev.category_set_path is None


class TestExtractIdxByInt:
    @pytest.fixture
    def ev(self):
        return _make_ev()

    def test_zero(self, ev):
        assert ev.extract_idx_by_int(0) == []

    def test_one(self, ev):
        assert ev.extract_idx_by_int(1) == [0]

    def test_two(self, ev):
        assert ev.extract_idx_by_int(2) == [1]

    def test_three(self, ev):
        assert ev.extract_idx_by_int(3) == [0, 1]

    def test_five(self, ev):
        assert ev.extract_idx_by_int(5) == [0, 2]

    def test_seven(self, ev):
        assert ev.extract_idx_by_int(7) == [0, 1, 2]

    def test_fifteen(self, ev):
        assert ev.extract_idx_by_int(15) == [0, 1, 2, 3]

    def test_twelve(self, ev):
        assert ev.extract_idx_by_int(12) == [2, 3]


class TestResultPath:
    def test_result_path_without_tag(self):
        ev = _make_ev(
            input_path=Path('/test/annotated.vcf.gz'),
            output_dir_path=Path('/output'),
            tag=None,
        )
        assert 'extracted_variants.txt.gz' in str(ev.result_path)

    def test_result_path_with_tag(self):
        ev = _make_ev(
            input_path=Path('/test/annotated.vcf.gz'),
            output_dir_path=Path('/output'),
            tag='mytag',
        )
        result = str(ev.result_path)
        assert 'mytag' in result
        assert 'extracted_variants.txt.gz' in result
