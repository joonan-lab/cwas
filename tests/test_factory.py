import argparse
import pytest
from pathlib import Path
import cwas.factory
from cwas.runnable import Runnable


def test_make_class_name():
    assert cwas.factory.make_class_name("categorization") == "Categorization"
    assert cwas.factory.make_class_name("annotation") == "Annotation"
    assert cwas.factory.make_class_name("burden_test") == "BurdenTest"
    assert cwas.factory.make_class_name("permutation_test") == "PermutationTest"
    assert cwas.factory.make_class_name("burden_shift") == "BurdenShift"
    assert cwas.factory.make_class_name("risk_score") == "RiskScore"
    assert cwas.factory.make_class_name("dawn") == "Dawn"


def test_create_factory():
    input_path = Path(__file__).parent / "test_file/test_annotated.vcf.gz"
    output_dir_path = Path(__file__).parent / "test_file"
    cat_path = Path(__file__).parent / "test_file/test.categorization_result.zarr"
    sample_info_path = Path(__file__).parent / "test_file/sample.txt"
    adj_factor_path = Path(__file__).parent / "test_file/adj_factors.txt"

    factory_inst = cwas.factory.create("categorization")
    args = argparse.Namespace(num_proc=1, input_path=str(input_path),
                              output_dir_path=str(output_dir_path))
    assert isinstance(factory_inst.argparser(), argparse.ArgumentParser)
    assert isinstance(factory_inst.runnable(args), Runnable)

    factory_inst = cwas.factory.create("binomial_test")
    args = argparse.Namespace(num_proc=1, cat_path=str(cat_path),
                              output_dir_path=str(output_dir_path),
                              sample_info_path = str(sample_info_path),
                              adj_factor_path = str(adj_factor_path),
                              use_n_carrier = False)
    assert isinstance(factory_inst.argparser(), argparse.ArgumentParser)
    assert isinstance(factory_inst.runnable(args), Runnable)

def test_create_factory_with_invalid_step():
    with pytest.raises(ValueError):
        cwas.factory.create("burden_test")