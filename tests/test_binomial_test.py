import pandas as pd
import pytest
from cwas.binomial_test import BinomialTest


class BinomialTestMock(BinomialTest):
    """This class do not make outputs"""

    def save_result(self):
        pass

    def update_env(self):
        pass


@pytest.fixture(scope="module")
def sample_info_path(cwas_workspace):
    return cwas_workspace / "samples.txt"


@pytest.fixture(scope="module")
def adj_factor_path(cwas_workspace):
    return cwas_workspace / "adj_factors.txt"


@pytest.fixture(scope="module")
def args(sample_info_path, adj_factor_path):
    return ["-s", str(sample_info_path), "-a", str(adj_factor_path)]


@pytest.fixture(scope="module", autouse=True)
def setup_and_teardown(cwas_workspace, sample_info_path, adj_factor_path):
    cwas_workspace.mkdir()
    sample_info_path.touch()
    adj_factor_path.touch()
    yield
    sample_info_path.unlink()
    adj_factor_path.unlink()
    cwas_workspace.rmdir()


@pytest.fixture
def categorization_result():
    results = [
        {"SAMPLE": "Sample1", "A_B_C_D_E": 5, "a_b_c_d_e": 4},
        {"SAMPLE": "Sample2", "A_B_C_D_E": 10, "a_b_c_d_e": 16},
        {"SAMPLE": "Sample3", "A_B_C_D_E": 12, "a_b_c_d_e": 8},
        {"SAMPLE": "Sample4", "A_B_C_D_E": 7, "a_b_c_d_e": 11},
    ]
    return pd.DataFrame(results).set_index("SAMPLE")


@pytest.fixture
def sample_info():
    samples = [
        {"SAMPLE": "Sample1", "FAMILY": "F1", "PHENOTYPE": "case"},
        {"SAMPLE": "Sample2", "FAMILY": "F1", "PHENOTYPE": "ctrl"},
        {"SAMPLE": "Sample3", "FAMILY": "F2", "PHENOTYPE": "case"},
        {"SAMPLE": "Sample4", "FAMILY": "F2", "PHENOTYPE": "ctrl"},
    ]
    return pd.DataFrame(samples).set_index("SAMPLE")


@pytest.fixture
def adjustment_factor():
    adj_factors = [
        {"SAMPLE": "Sample1", "AdjustFactor": 2.0},
        {"SAMPLE": "Sample2", "AdjustFactor": 0.5},
        {"SAMPLE": "Sample3", "AdjustFactor": 0.25},
        {"SAMPLE": "Sample4", "AdjustFactor": 1.0},
    ]
    return pd.DataFrame(adj_factors).set_index("SAMPLE")


@pytest.fixture
def sample_info_other_sample():
    samples = [
        {"SAMPLE": "SampleA", "FAMILY": "F1", "PHENOTYPE": "case"},
        {"SAMPLE": "SampleB", "FAMILY": "F2", "PHENOTYPE": "ctrl"},
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
    args, categorization_result, sample_info, adjustment_factor,
):
    # This is not an appropriate usage.
    inst = BinomialTestMock.get_instance(args)
    inst._categorization_result = categorization_result
    inst._sample_info = sample_info
    inst._adj_factor = adjustment_factor
    return inst


@pytest.fixture
def binomial_test_with_inconsistent_sample(
    args,
    categorization_result,
    sample_info_other_sample,
    adjustment_factor_other_sample,
):
    # This is not an appropriate usage.
    inst = BinomialTestMock.get_instance(args)
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
        "conservation",
        "gencode",
        "region",
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
    assert binomial_test.binom_p == 0.5  # A fraction of cases


def test_case_cnt(binomial_test):
    assert binomial_test.case_cnt == 2


def test_ctrl_cnt(binomial_test):
    assert binomial_test.ctrl_cnt == 2


def test_case_variant_cnt(binomial_test):
    binomial_test._adjust_categorization_result()
    assert list(binomial_test.case_variant_cnt) == [13, 10]


def test_ctrl_variant_cnt(binomial_test):
    binomial_test._adjust_categorization_result()
    assert list(binomial_test.ctrl_variant_cnt) == [12, 19]
