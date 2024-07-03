import pandas as pd
import pytest
from cwas.binomial_test import BinomialTest
import sys
import cwas.cli
from pathlib import Path
import zarr
import argparse

class BinomialTestMock(BinomialTest):
    def __init__(self, args):
        super().__init__(args)
        self._raw_counts = None

    def save_result(self):
        pass

    def update_env(self):
        pass

@pytest.fixture(scope="module")
def cat_path():
    return Path(__file__).parent / "test_file/test.categorization_result.zarr"

@pytest.fixture
def output_dir_path():
    return Path(__file__).parent / "test_file"

@pytest.fixture
def sample_info_path():
    return Path(__file__).parent / "test_file/sample.txt"

@pytest.fixture
def adj_factor_path():
    return Path(__file__).parent / "test_file/adj_factors.txt"

@pytest.fixture
def sample_info_other_sample_path():
    return Path(__file__).parent / "test_file/sample_other.txt"

@pytest.fixture
def adj_factor_other_sample_path():
    return Path(__file__).parent / "test_file/adj_factors_other.txt"

@pytest.fixture
def binomial_test(cat_path, output_dir_path, sample_info_path, adj_factor_path):
    # sys.argv = ['cwas', 'binomial_test', '-i', str(cat_path), '-s', str(sample_info_path), '-a', str(adj_factor_path)]
    # binom_inst = cwas.cli.main()
    args = argparse.Namespace(num_proc=1,
                              cat_path=cat_path,
                              output_dir_path = output_dir_path,
                              sample_info_path=str(sample_info_path),
                              adj_factor_path=str(adj_factor_path),
                              use_n_carrier=False,
                              eff_test = False,
                              tag = None,
                              plot_size = 7, plot_title = "Binomial test result", font_size = 15, marker_size = 15)
    inst = BinomialTestMock(args)
    root = zarr.open(str(cat_path), mode='r')
    inst._categorization_result = pd.DataFrame(data=root['data'],
                                               columns=root['metadata'].attrs['category'],
                                               index=root['metadata'].attrs['sample_id'])
    inst._sample_info = pd.read_csv(sample_info_path, sep='\t').set_index("SAMPLE")
    inst._adj_factor = pd.read_csv(adj_factor_path, sep='\t').set_index("SAMPLE")
    inst.save_counts_table('raw')  # Ensure _raw_counts is initialized
    return inst

@pytest.fixture
def binomial_test_with_inconsistent_sample(cat_path, output_dir_path, sample_info_other_sample_path, adj_factor_other_sample_path):
    # sys.argv = ['cwas', 'binomial_test', '-i', str(cat_path), '-s', str(sample_info_other_sample_path), '-a', str(adj_factor_other_sample_path)]
    # binom_inst = cwas.cli.main()
    args = argparse.Namespace(num_proc=1, cat_path=cat_path,
                              output_dir_path = output_dir_path,
                              sample_info_path=str(sample_info_other_sample_path),
                              adj_factor_path=str(adj_factor_other_sample_path),
                              use_n_carrier=False,
                              eff_test = False,
                              tag = None,
                              plot_size = 7, plot_title = "Binomial test result", font_size = 15, marker_size = 15)
    inst = BinomialTestMock(args)
    root = zarr.open(str(cat_path), mode='r')
    inst._categorization_result = pd.DataFrame(data=root['data'],
                                               columns=root['metadata'].attrs['category'],
                                               index=root['metadata'].attrs['sample_id'])
    inst._sample_info = pd.read_csv(sample_info_other_sample_path, sep='\t').set_index("SAMPLE")
    inst._adj_factor = pd.read_csv(adj_factor_other_sample_path, sep='\t').set_index("SAMPLE")
    inst.save_counts_table('raw')  # Ensure _raw_counts is initialized
    return inst

def test_adjust_categorization_result(binomial_test):
    assert binomial_test._categorization_result is not None
    assert not pd.DataFrame(binomial_test._categorization_result).empty

def test_adjust_categorization_with_inconsistent_sample(binomial_test_with_inconsistent_sample):
    with pytest.raises(ValueError):
        binomial_test_with_inconsistent_sample._adjust_categorization_result()

def test_run_with_inconsistent_sample(binomial_test_with_inconsistent_sample):
    with pytest.raises(ValueError):
        binomial_test_with_inconsistent_sample.run()

def test_run(binomial_test):
    binomial_test._adjust_categorization_result()
    binomial_test.run()
    assert binomial_test._result is not None
    assert binomial_test._result.index.name == "Category"
    expected_columns = [
        "variant_type",
        "gene_set",
        "functional_score",
        "gencode",
        "functional_annotation",
        "Case_DNV_Count",
        "Ctrl_DNV_Count",
        "Relative_Risk",
        "P",
        "P_1side",
        "Z_1side",
    ]
    expected_index = [
        'a_b_c_d_e',
        'b_c_d_e_f',
        'c_d_e_f_g',
        'd_e_f_g_h',
        'e_f_g_h_i',
        'f_g_h_i_j',
        'g_h_i_j_k',
        'h_i_j_k_l',
        'i_j_k_l_m',
        'j_k_l_m_n',
    ]
    assert list(binomial_test._result.columns.values) == expected_columns
    assert list(binomial_test._result.index.values) == expected_index

def test_binom_p(binomial_test):
    assert binomial_test.binom_p == 1 / 2  # A fraction of cases

def test_case_cnt(binomial_test):
    assert binomial_test.case_cnt == 3

def test_ctrl_cnt(binomial_test):
    assert binomial_test.ctrl_cnt == 3

def test_case_variant_cnt(binomial_test):
    binomial_test._adjust_categorization_result()
    assert list(binomial_test.case_variant_cnt) == [85.0, 29.8, 76.0, 63.800000000000004, 106.30000000000001, 24.8, 37.9, 86.9, 22.4, 72.6]

def test_ctrl_variant_cnt(binomial_test):
    binomial_test._adjust_categorization_result()
    assert list(binomial_test.ctrl_variant_cnt) == [59.8, 121.9, 166.3, 279.1, 249.7, 125.5, 88.3, 149.0, 140.3, 147.1]

def test_calculate_relative_risk(binomial_test):
    binomial_test._adjust_categorization_result()
    binomial_test.run()
    expected_relative_risk1 = (85.0 / 3) / (59.8 / 3)  # a_b_c_d_e
    expected_relative_risk2 = (29.8 / 3) / (121.9 / 3)  # b_c_d_e_f
    expected_relative_risk3 = (76.0 / 3) / (166.3 / 3)  # c_d_e_f_g
    expected_relative_risk4 = 0.22859190254389108 # (63.8 / 3) / (279.1 / 3)  # d_e_f_g_h
    expected_relative_risk5 = 0.4257108530236284 # (106.3 / 3) / (249.7 / 3)  # e_f_g_h_i
    expected_relative_risk6 = (24.8 / 3) / (125.5 / 3)  # f_g_h_i_j
    expected_relative_risk7 = 0.42921857304643257 # (37.9 / 3) / (88.8 / 3)  # g_h_i_j_k
    expected_relative_risk8 = (86.9 / 3) / (149.0 / 3)  # h_i_j_k_l
    expected_relative_risk9 = 0.15965787598004272 # (22.4 / 3) / (140.8 / 3)  # i_j_k_l_m
    expected_relative_risk10 = (72.6 / 3) / (147.1 / 3)  # j_k_l_m_n

    assert binomial_test._result["Relative_Risk"].to_list() == [
        expected_relative_risk1,
        expected_relative_risk2,
        expected_relative_risk3,
        expected_relative_risk4,
        expected_relative_risk5,
        expected_relative_risk6,
        expected_relative_risk7,
        expected_relative_risk8,
        expected_relative_risk9,
        expected_relative_risk10,
    ]
