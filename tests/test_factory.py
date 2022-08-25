import argparse
import cwas.factory
import pytest

from cwas.runnable import Runnable


def test_make_class_name():
    assert cwas.factory.make_class_name("categorization") == "Categorization"
    assert cwas.factory.make_class_name("annotation") == "Annotation"
    assert cwas.factory.make_class_name("burden_test") == "BurdenTest"


def test_get_runnable():
    runnable = cwas.factory.get_runnable("categorization")
    assert runnable.__name__ == "Categorization"
    assert hasattr(runnable, "run")

    runnable = cwas.factory.get_runnable("burden_test")
    assert runnable.__name__ == "BurdenTest"
    assert hasattr(runnable, "run")

    with pytest.raises(ValueError):
        cwas.factory.get_runnable("not_runnable")


def test_create_factory():
    factory_inst = cwas.factory.create("categorization")
    assert isinstance(factory_inst.argparser(), argparse.ArgumentParser)
    assert isinstance(factory_inst.runnable(None), Runnable)

    factory_inst = cwas.factory.create("binomial_test")
    assert isinstance(factory_inst.argparser(), argparse.ArgumentParser)
    assert isinstance(factory_inst.runnable(None), Runnable)


def test_create_factory_with_invalid_step():
    with pytest.raises(ValueError):
        cwas.factory.create("burden_test")
