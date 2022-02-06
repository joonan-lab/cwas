import cwas.core.categorization.utils as utils


def test_get_idx_dict():
    assert utils.get_idx_dict(["a", "b", "c"]) == {"a": 0, "b": 1, "c": 2}


def test_extract_sublist_by_int():
    test_list = ["a", "b", "c", "d"]
    assert utils.extract_sublist_by_int(test_list, 0) == []
    assert utils.extract_sublist_by_int(test_list, 1) == ["a"]
    assert utils.extract_sublist_by_int(test_list, 2) == ["b"]
    assert utils.extract_sublist_by_int(test_list, 4) == ["c"]
    assert utils.extract_sublist_by_int(test_list, 8) == ["d"]
    assert utils.extract_sublist_by_int(test_list, 3) == ["a", "b"]
    assert utils.extract_sublist_by_int(test_list, 11) == ["a", "b", "d"]
    assert utils.extract_sublist_by_int(test_list, 12) == ["c", "d"]
    assert utils.extract_sublist_by_int(test_list, 13) == ["a", "c", "d"]
