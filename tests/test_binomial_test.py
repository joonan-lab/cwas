import pandas as pd
import pytest
from cwas.binomial_test import BinomialTest
import sys
import cwas.cli
import zarr
from pathlib import Path
import shutil

class BinomialTestMock(BinomialTest):
    """This class do not make outputs"""

    def save_result(self):
        pass

    def update_env(self):
        pass

# Fixture to create a temporary directory for the test
@pytest.fixture(scope="module")
def cat_path(cwas_workspace):
    cat_path = cwas_workspace / "categorized.zarr"
    yield cat_path  # fixture를 사용하고, 여기서 제어를 반환함

    # 여기에서 정리 작업을 수행
    if cat_path.exists():
        shutil.rmtree(cat_path)
        shutil.rmtree(cwas_workspace)

@pytest.fixture
def categorization_result():
    results = [
        {"SAMPLE": "Sample1", "A_B_C_D_E": 5, "a_b_c_d_e": 4},
        {"SAMPLE": "Sample2", "A_B_C_D_E": 10, "a_b_c_d_e": 16},
        {"SAMPLE": "Sample3", "A_B_C_D_E": 12, "a_b_c_d_e": 8},
        {"SAMPLE": "Sample4", "A_B_C_D_E": 7, "a_b_c_d_e": 11},
        {"SAMPLE": "Sample5", "A_B_C_D_E": 8, "a_b_c_d_e": 6},
        {"SAMPLE": "Sample6", "A_B_C_D_E": 15, "a_b_c_d_e": 10},
    ]
    return pd.DataFrame(results).set_index("SAMPLE")


@pytest.fixture
def create_zarr_group(cat_path, categorization_result):
    # TODO: Add VCF entries (variants)
    root = zarr.open(cat_path, mode='w')
    root.create_group('metadata')
    root['metadata'].attrs['sample_id'] = categorization_result.index.tolist()
    root['metadata'].attrs['category'] = categorization_result.columns.tolist()
    root.create_dataset('data', data=categorization_result.values, chunks=(1000, 1000), dtype='i4')
    return cat_path

# Test function that uses the created Zarr group
def test_zarr_data(create_zarr_group):
    cat_path = create_zarr_group
    # Your test logic that uses cat_path
    assert Path(cat_path).is_dir()  # Example test assertion

@pytest.fixture
def sample_info():
    samples = [
        {"SAMPLE": "Sample1", "PHENOTYPE": "case"},
        {"SAMPLE": "Sample2", "PHENOTYPE": "ctrl"},
        {"SAMPLE": "Sample3", "PHENOTYPE": "case"},
        {"SAMPLE": "Sample4", "PHENOTYPE": "ctrl"},
        {"SAMPLE": "Sample5", "PHENOTYPE": "ctrl"},
        {"SAMPLE": "Sample6", "PHENOTYPE": "ctrl"},
    ]
    return pd.DataFrame(samples).set_index("SAMPLE")


@pytest.fixture
def adjustment_factor():
    adj_factors = [
        {"SAMPLE": "Sample1", "AdjustFactor": 2.0},
        {"SAMPLE": "Sample2", "AdjustFactor": 0.5},
        {"SAMPLE": "Sample3", "AdjustFactor": 0.25},
        {"SAMPLE": "Sample4", "AdjustFactor": 1.0},
        {"SAMPLE": "Sample5", "AdjustFactor": 0.5},
        {"SAMPLE": "Sample6", "AdjustFactor": 0.2},
    ]
    return pd.DataFrame(adj_factors).set_index("SAMPLE")


@pytest.fixture
def sample_info_other_sample():
    samples = [
        {"SAMPLE": "SampleA", "PHENOTYPE": "case"},
        {"SAMPLE": "SampleB", "PHENOTYPE": "ctrl"},
    ]
    return pd.DataFrame(samples).set_index("SAMPLE")


@pytest.fixture
def adjustment_factor_other_sample():
    adj_factors = [
        {"SAMPLE": "SampleA", "AdjustFactor": 2.0},
        {"SAMPLE": "SampleB", "AdjustFactor": 0.5},
    ]
    return pd.DataFrame(adj_factors).set_index("SAMPLE")


@pytest.fixture
def binomial_test(
    cat_path, categorization_result, sample_info, adjustment_factor,
):
    # This is not an appropriate usage.
    factory_inst = cwas.factory.create("binomial_test")
    #sys.argv = ['cwas', 'binomial_test', None]
    #binom_inst = cwas.cli.main()
    #factory_inst.argparser().parse_args(['do'])
    #args = factory_inst.argparser().parse_args(['-i', str(cat_path)])
    args = factory_inst.argparser().parse_args([])
    
    args.cat_path = cat_path
    inst = BinomialTestMock(factory_inst.runnable(factory_inst.argparser()))
    inst._categorization_result = categorization_result
    inst._sample_info = sample_info
    inst._adj_factor = adjustment_factor
    return inst


@pytest.fixture
def binomial_test_with_inconsistent_sample(
    cat_path,
    categorization_result,
    sample_info_other_sample,
    adjustment_factor_other_sample,
):
    # This is not an appropriate usage.
    factory_inst = cwas.factory.create("binomial_test")
    #sys.argv = ['cwas', 'binomial_test', None]
    #binom_inst = cwas.cli.main()
    factory_inst._categorization_result = categorization_result
    factory_inst._cat_path = cat_path
    
    inst = BinomialTestMock(factory_inst.runnable(None))
    #inst = BinomialTestMock()
    inst._categorization_result = categorization_result
    inst._sample_info = sample_info_other_sample
    inst._adj_factor = adjustment_factor_other_sample
    return inst


def test_adjust_categorization_result(binomial_test):
    binomial_test._adjust_categorization_result()
    categorization_result = binomial_test.categorization_result
    assert categorization_result.loc["Sample1"].to_list() == [10, 8]
    assert categorization_result.loc["Sample2"].to_list() == [5, 8]
    assert categorization_result.loc["Sample3"].to_list() == [3, 2]
    assert categorization_result.loc["Sample4"].to_list() == [7, 11]
    assert categorization_result.loc["Sample5"].to_list() == [4, 3]
    assert categorization_result.loc["Sample6"].to_list() == [3, 2]


def test_adjust_categorization_with_inconsistent_sample(
    binomial_test_with_inconsistent_sample,
):
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
        "gene_list",
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
        "A_B_C_D_E",
        "a_b_c_d_e",
    ]
    assert list(binomial_test._result.columns.values) == expected_columns
    assert list(binomial_test._result.index.values) == expected_index


def test_binom_p(binomial_test):
    assert binomial_test.binom_p == 1 / 3  # A fraction of cases


def test_case_cnt(binomial_test):
    assert binomial_test.case_cnt == 2


def test_ctrl_cnt(binomial_test):
    assert binomial_test.ctrl_cnt == 4


def test_case_variant_cnt(binomial_test):
    binomial_test._adjust_categorization_result()
    assert list(binomial_test.case_variant_cnt) == [13, 10]


def test_ctrl_variant_cnt(binomial_test):
    binomial_test._adjust_categorization_result()
    assert list(binomial_test.ctrl_variant_cnt) == [19, 24]


def test_calculate_relative_risk(binomial_test):
    binomial_test._adjust_categorization_result()
    binomial_test.run()
    expected_relative_risk1 = (13 / 2) / (19 / 4)  # A_B_C_D_E
    expected_relative_risk2 = (10 / 2) / (24 / 4)  # a_b_c_d_e
    assert binomial_test._result["Relative_Risk"].to_list() == [
        expected_relative_risk1,
        expected_relative_risk2,
    ]

