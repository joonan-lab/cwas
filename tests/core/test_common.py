"""
Test the methods in cwas.core.common
"""
import numpy as np
import pytest

import cwas.core.common as common


def test_cmp_two_arr():
    # Empty arrays
    arr1 = np.array([])
    arr2 = np.array([])
    assert common.cmp_two_arr(arr1, arr2)

    # Different lengths
    arr1 = np.array([1, 2, 3, 4, 5])
    arr2 = np.array([3, 4, 5])
    assert not common.cmp_two_arr(arr1, arr2)

    # Same arrays
    arr2 = np.array([1, 2, 3, 4, 5])
    assert common.cmp_two_arr(arr1, arr2)

    # Arrays with the same items but not in the same order
    np.random.shuffle(arr2)
    assert common.cmp_two_arr(arr1, arr2)

    # False if the items of each list are not the same.
    arr2 = np.array([1, 14, 3, 6, 5])
    assert not common.cmp_two_arr(arr1, arr2)


def test_div_dist_num():
    # Raise ValueError if one of any argument is not a positive integer.
    with pytest.raises(ValueError):
        _ = common.div_dist_num(-1, 5)

    with pytest.raises(ValueError):
        _ = common.div_dist_num(5, -1)

    # Raise ValueError if 'num' is less than 'num_group'.
    with pytest.raises(ValueError):
        _ = common.div_dist_num(4, 5)

    num = 5
    num_group = 5
    expected = [1, 1, 1, 1, 1]
    assert expected == common.div_dist_num(num, num_group)

    for i in range(5):
        num += 1
        expected[i] += 1
        assert expected == common.div_dist_num(num, num_group)


def test_chunk_list():
    # ValueError will be raised if the input list is empty.
    with pytest.raises(ValueError):
        _ = common.chunk_list([], 5)

    _list = [1, 2, 3, 4, 5]

    # ValueError will be raised
    # if the number of lists is not a positive integer.
    with pytest.raises(ValueError):
        _ = common.chunk_list(_list, -1)

    with pytest.raises(ValueError):
        _ = common.chunk_list(_list, 0)

    result = common.chunk_list(_list, 1)
    expected = [[1, 2, 3, 4, 5]]
    assert result == expected

    result = common.chunk_list(_list, 2)
    expected = [[1, 2, 3], [4, 5]]
    assert result == expected

    result = common.chunk_list(_list, 3)
    expected = [[1, 2], [3, 4], [5]]
    assert result == expected

    result = common.chunk_list(_list, 4)
    expected = [[1, 2], [3], [4], [5]]
    assert result == expected

    result = common.chunk_list(_list, 5)
    expected = [[1], [2], [3], [4], [5]]
    assert result == expected

    # ValueError is raised if the number of lists exceeds the length of the
    # input list.
    with pytest.raises(ValueError):
        _ = common.chunk_list(_list, 6)
