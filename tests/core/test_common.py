"""
Test the methods in cwas.core.common
"""
import numpy as np

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
