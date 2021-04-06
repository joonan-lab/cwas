import pytest

import cwas.factory as factory


def test_make_class_name():
    assert factory.make_class_name('categorization') == 'Categorization'
    assert factory.make_class_name('annotation') == 'Annotation'
    assert factory.make_class_name('burden_test') == 'BurdenTest'


def test_get_runnable():
    runnable = factory.get_runnable('categorization')
    assert runnable.__name__ == 'Categorization'
    assert hasattr(runnable, 'run')

    runnable = factory.get_runnable('burden_test')
    assert runnable.__name__ == 'BurdenTest'
    assert hasattr(runnable, 'run')

    with pytest.raises(ValueError):
        factory.get_runnable('not_runnable')
