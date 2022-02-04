import pytest
from cwas.core.categorization.category import Category


def test_category_hashable():
    category1 = Category("A1", "B1", "C1", "D1", "E1")
    category2 = Category("A2", "B2", "C2", "D2", "E2")
    category3 = Category("A1", "B1", "C1", "D1", "E1")

    test_dict = {}
    test_dict[category1] = 5
    test_dict[category2] = 2
    test_dict[category3] += 1

    assert test_dict[category1] == 6
    assert test_dict[category2] == 2


def test_category_hashable_keyerror():
    category1 = Category("A1", "B1", "C1", "D1", "E1")
    category2 = Category("A1", "B2", "C1", "D2", "E1")

    test_dict = {}
    test_dict[category1] = 1

    with pytest.raises(KeyError):
        test_dict[category2] += 1


def test_get_category_from_dict():
    category1 = Category("A1", "B1", "C1", "D1", "E1")
    test_dict = {}
    test_dict[category1] = 1

    category_from_dict = list(test_dict.keys())[0]
    assert category1 == category_from_dict

