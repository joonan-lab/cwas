"""
Tests for cwas.argparser module
"""
import pytest
import argparse
from pathlib import Path
from unittest.mock import patch
from cwas import argparser


class TestStartParser:
    def test_start_returns_argument_parser(self):
        parser = argparser.start()
        assert isinstance(parser, argparse.ArgumentParser)

    def test_start_workspace_default_is_home_cwas(self):
        parser = argparser.start()
        args = parser.parse_args([])
        expected = Path.home() / ".cwas"
        assert args.workspace == expected

    def test_start_workspace_short_flag(self):
        parser = argparser.start()
        args = parser.parse_args(['-w', '/custom/path'])
        assert args.workspace == Path('/custom/path')

    def test_start_workspace_long_flag(self):
        parser = argparser.start()
        args = parser.parse_args(['--workspace', '/custom/path'])
        assert args.workspace == Path('/custom/path')


class TestConfigurationParser:
    def test_configuration_returns_argument_parser(self):
        parser = argparser.configuration()
        assert isinstance(parser, argparse.ArgumentParser)

    def test_configuration_force_overwrite_default(self):
        parser = argparser.configuration()
        args = parser.parse_args([])
        assert args.force_overwrite == 0

    def test_configuration_force_overwrite_flag(self):
        parser = argparser.configuration()
        args = parser.parse_args(['-f'])
        assert args.force_overwrite == 1

    def test_configuration_has_all_optional_arguments(self):
        parser = argparser.configuration()
        args = parser.parse_args([
            '-d', '/data',
            '-m', '/matrix',
            '-a', '/config.yaml',
            '-v', '/vep',
            '-f'
        ])
        assert args.data_dir == Path('/data')
        assert args.gene_matrix == Path('/matrix')
        assert args.annot_key_conf == Path('/config.yaml')
        assert args.vep == Path('/vep')
        assert args.force_overwrite == 1


class TestPreparationParser:
    def test_preparation_returns_argument_parser(self):
        parser = argparser.preparation()
        assert isinstance(parser, argparse.ArgumentParser)

    def test_preparation_num_proc_default(self):
        parser = argparser.preparation()
        args = parser.parse_args([])
        assert args.num_proc == 1

    def test_preparation_num_proc_custom(self):
        parser = argparser.preparation()
        args = parser.parse_args(['-p', '8'])
        assert args.num_proc == 8

    def test_preparation_force_overwrite_default(self):
        parser = argparser.preparation()
        args = parser.parse_args([])
        assert args.force_overwrite == 0


class TestAnnotationParser:
    @patch('dotenv.dotenv_values')
    def test_annotation_vcf_path_required(self, mock_dotenv):
        mock_dotenv.return_value = {"CWAS_WORKSPACE": "/default"}
        parser = argparser.annotation()
        with pytest.raises(SystemExit):
            parser.parse_args([])

    @patch('dotenv.dotenv_values')
    def test_annotation_vcf_path_provided(self, mock_dotenv):
        mock_dotenv.return_value = {"CWAS_WORKSPACE": "/default"}
        parser = argparser.annotation()
        args = parser.parse_args(['-v', '/path/to/file.vcf'])
        assert args.vcf_path == Path('/path/to/file.vcf')

    @patch('dotenv.dotenv_values')
    def test_annotation_num_proc_default(self, mock_dotenv):
        mock_dotenv.return_value = {"CWAS_WORKSPACE": "/default"}
        parser = argparser.annotation()
        args = parser.parse_args(['-v', '/file.vcf'])
        assert args.num_proc == 1


class TestCategorizationParser:
    @patch('dotenv.dotenv_values')
    def test_categorization_input_path_required(self, mock_dotenv):
        mock_dotenv.return_value = {"CWAS_WORKSPACE": "/default"}
        parser = argparser.categorization()
        with pytest.raises(SystemExit):
            parser.parse_args([])

    @patch('dotenv.dotenv_values')
    def test_categorization_input_path_provided(self, mock_dotenv):
        mock_dotenv.return_value = {"CWAS_WORKSPACE": "/default"}
        parser = argparser.categorization()
        args = parser.parse_args(['-i', '/annotated.vcf'])
        assert args.input_path == Path('/annotated.vcf')

    @patch('dotenv.dotenv_values')
    def test_categorization_num_proc_default(self, mock_dotenv):
        mock_dotenv.return_value = {"CWAS_WORKSPACE": "/default"}
        parser = argparser.categorization()
        args = parser.parse_args(['-i', '/file.vcf'])
        assert args.num_proc == 1


class TestBinomialTestParser:
    @patch('dotenv.dotenv_values')
    def test_binomial_test_required_args(self, mock_dotenv):
        mock_dotenv.return_value = {"CWAS_WORKSPACE": "/default"}
        parser = argparser.binomial_test()
        with pytest.raises(SystemExit):
            parser.parse_args([])

    @patch('dotenv.dotenv_values')
    def test_binomial_test_with_required_args(self, mock_dotenv):
        mock_dotenv.return_value = {"CWAS_WORKSPACE": "/default"}
        parser = argparser.binomial_test()
        args = parser.parse_args(['-i', '/cat.zarr', '-s', '/sample_info.txt'])
        assert args.cat_path == Path('/cat.zarr')
        assert args.sample_info_path == Path('/sample_info.txt')

    @patch('dotenv.dotenv_values')
    def test_binomial_test_use_n_carrier_default(self, mock_dotenv):
        mock_dotenv.return_value = {"CWAS_WORKSPACE": "/default"}
        parser = argparser.binomial_test()
        args = parser.parse_args(['-i', '/cat.zarr', '-s', '/sample_info.txt'])
        assert args.use_n_carrier is False

    @patch('dotenv.dotenv_values')
    def test_binomial_test_use_n_carrier_flag(self, mock_dotenv):
        mock_dotenv.return_value = {"CWAS_WORKSPACE": "/default"}
        parser = argparser.binomial_test()
        args = parser.parse_args(['-i', '/cat.zarr', '-s', '/sample_info.txt', '-u'])
        assert args.use_n_carrier is True

    @patch('dotenv.dotenv_values')
    def test_binomial_test_plot_title_default(self, mock_dotenv):
        mock_dotenv.return_value = {"CWAS_WORKSPACE": "/default"}
        parser = argparser.binomial_test()
        args = parser.parse_args(['-i', '/cat.zarr', '-s', '/sample_info.txt'])
        assert args.plot_title == "Binomial test result"


class TestPermutationTestParser:
    @patch('dotenv.dotenv_values')
    def test_permutation_test_required_args(self, mock_dotenv):
        mock_dotenv.return_value = {"CWAS_WORKSPACE": "/default"}
        parser = argparser.permutation_test()
        with pytest.raises(SystemExit):
            parser.parse_args([])

    @patch('dotenv.dotenv_values')
    def test_permutation_test_num_perm_default(self, mock_dotenv):
        mock_dotenv.return_value = {"CWAS_WORKSPACE": "/default"}
        parser = argparser.permutation_test()
        args = parser.parse_args(['-i', '/cat.zarr', '-s', '/sample_info.txt'])
        assert args.num_perm == 10000

    @patch('dotenv.dotenv_values')
    def test_permutation_test_burden_shift_flag(self, mock_dotenv):
        mock_dotenv.return_value = {"CWAS_WORKSPACE": "/default"}
        parser = argparser.permutation_test()
        args = parser.parse_args(['-i', '/cat.zarr', '-s', '/sample_info.txt', '-b'])
        assert args.burden_shift is True


class TestBurdenShiftParser:
    @patch('dotenv.dotenv_values')
    def test_burden_shift_required_args(self, mock_dotenv):
        mock_dotenv.return_value = {"CWAS_WORKSPACE": "/default"}
        parser = argparser.burden_shift()
        with pytest.raises(SystemExit):
            parser.parse_args([])

    @patch('dotenv.dotenv_values')
    def test_burden_shift_with_required_args(self, mock_dotenv):
        mock_dotenv.return_value = {"CWAS_WORKSPACE": "/default"}
        parser = argparser.burden_shift()
        args = parser.parse_args([
            '-i', '/burden_test.txt', '-b', '/burden_res.txt',
            '-c_info', '/cat_info.txt', '-c_count', '/cat_count.txt'
        ])
        assert args.input_path == Path('/burden_test.txt')
        assert args.burden_res == Path('/burden_res.txt')

    @patch('dotenv.dotenv_values')
    def test_burden_shift_count_cutoff_default(self, mock_dotenv):
        mock_dotenv.return_value = {"CWAS_WORKSPACE": "/default"}
        parser = argparser.burden_shift()
        args = parser.parse_args([
            '-i', '/a.txt', '-b', '/b.txt', '-c_info', '/c.txt', '-c_count', '/d.txt'
        ])
        assert args.count_cutoff == 7

    @patch('dotenv.dotenv_values')
    def test_burden_shift_pval_default(self, mock_dotenv):
        mock_dotenv.return_value = {"CWAS_WORKSPACE": "/default"}
        parser = argparser.burden_shift()
        args = parser.parse_args([
            '-i', '/a.txt', '-b', '/b.txt', '-c_info', '/c.txt', '-c_count', '/d.txt'
        ])
        assert args.pval == 0.05
