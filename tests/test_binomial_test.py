import pandas as pd
import pytest
from cwas.binomial_test import BinomialTest


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
        {"SAMPLE": "Sample1", "Category1": 5, "Category2": 4},
        {"SAMPLE": "Sample2", "Category1": 10, "Category2": 16},
    ]
    return pd.DataFrame(results).set_index("SAMPLE")


@pytest.fixture
def sample_info():
    samples = [
        {"SAMPLE": "Sample1", "FAMILY": "F1", "PHENOTYPE": "case"},
        {"SAMPLE": "Sample2", "FAMILY": "F2", "PHENOTYPE": "ctrl"},
    ]
    return pd.DataFrame(samples).set_index("SAMPLE")


@pytest.fixture
def adjustment_factor():
    adj_factors = [
        {"SAMPLE": "Sample1", "AdjustFactor": 2.0},
        {"SAMPLE": "Sample2", "AdjustFactor": 0.5},
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
    inst = BinomialTest.get_instance(args)
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
    inst = BinomialTest.get_instance(args)
    inst._categorization_result = categorization_result
    inst._sample_info = sample_info_other_sample
    inst._adj_factor = adjustment_factor_other_sample
    return inst


def test_adjust_categorization_result(binomial_test):
    binomial_test._adjust_categorization_result()
    categorization_result = binomial_test.categorization_result
    assert categorization_result.loc["Sample1"].to_list() == [10, 8]
    assert categorization_result.loc["Sample2"].to_list() == [5, 8]


def test_adjust_categorization_with_inconsistent_sample(
    binomial_test_with_inconsistent_sample,
):
    with pytest.raises(ValueError):
        binomial_test_with_inconsistent_sample._adjust_categorization_result()


def test_run_with_inconsistent_sample(binomial_test_with_inconsistent_sample):
    with pytest.raises(ValueError):
        binomial_test_with_inconsistent_sample.run()


def test_run_burden_test(binomial_test):
    binomial_test._adjust_categorization_result()
    binomial_test.run_burden_test()
    assert binomial_test._result is not None
    assert binomial_test._result.index.name == "Category"
    expected_columns = [
        "Case_DNV_Count",
        "Ctrl_DNV_Count",
        "Relative_Risk",
        "P",
        "P_1side",
        "Z_1side",
    ]
    assert list(binomial_test._result.columns.values) == expected_columns


def test_binom_p(binomial_test):
    samples = [
        {"SAMPLE": "Sample1", "FAMILY": "F1", "PHENOTYPE": "case"},
        {"SAMPLE": "Sample2", "FAMILY": "F1", "PHENOTYPE": "case"},
        {"SAMPLE": "Sample3", "FAMILY": "F1", "PHENOTYPE": "case"},
        {"SAMPLE": "Sample4", "FAMILY": "F1", "PHENOTYPE": "case"},
        {"SAMPLE": "Sample5", "FAMILY": "F2", "PHENOTYPE": "ctrl"},
        {"SAMPLE": "Sample6", "FAMILY": "F2", "PHENOTYPE": "ctrl"},
        {"SAMPLE": "Sample7", "FAMILY": "F2", "PHENOTYPE": "ctrl"},
        {"SAMPLE": "Sample8", "FAMILY": "F2", "PHENOTYPE": "ctrl"},
        {"SAMPLE": "Sample9", "FAMILY": "F2", "PHENOTYPE": "ctrl"},
        {"SAMPLE": "Sample10", "FAMILY": "F2", "PHENOTYPE": "ctrl"},
    ]
    binomial_test._sample_info = pd.DataFrame(samples).set_index("SAMPLE")
    assert binomial_test.binom_p == 0.4  # A fraction of cases

