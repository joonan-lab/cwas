"""
Test the methods in cwas.core.common
"""
import cwas.core.common as common
import numpy as np
import pytest


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


def test_swap_label():
    # Make an input for testing
    label_group_pairs = []

    for n in range(1, 6):
        group_id = f"g{n}"
        label_group_pairs.append((f"{group_id}_l1", group_id))
        label_group_pairs.append((f"{group_id}_l2", group_id))

    np.random.shuffle(label_group_pairs)

    labels = np.array([pair[0] for pair in label_group_pairs])
    group_ids = np.array([pair[1] for pair in label_group_pairs])

    # Basic testing
    swap_labels = common.swap_label(labels, group_ids)
    assert len(swap_labels) == len(group_ids)
    assert np.any(labels != swap_labels)  # Swapping is succeeded.

    # Check if the labels are swapped within their groups.
    for swap_label, group_id in zip(swap_labels, group_ids):
        assert swap_label.startswith(group_id)

    # Test the case where more than 2 labels are in the same group.
    group_id = f"g{np.random.randint(1, 6)}"
    new_label = f"{group_id}_l3"
    labels = np.append(labels, new_label)
    group_ids = np.append(group_ids, group_id)
    np.random.shuffle(labels)
    np.random.shuffle(group_ids)

    with pytest.raises(AssertionError) as e_info:
        assert group_id in str(e_info.value)


def test_int_to_bit_arr():
    assert common.int_to_bit_arr(5, 0).size == 0
    assert common.int_to_bit_arr(5, 3).size == 3
    assert common.int_to_bit_arr(5, 5).size == 5
    assert (common.int_to_bit_arr(1, 3) == np.array([1, 0, 0])).all()
    assert (common.int_to_bit_arr(2, 3) == np.array([0, 1, 0])).all()
    assert (common.int_to_bit_arr(3, 3) == np.array([1, 1, 0])).all()
    assert (common.int_to_bit_arr(4, 3) == np.array([0, 0, 1])).all()
    assert (common.int_to_bit_arr(5, 3) == np.array([1, 0, 1])).all()
    assert (common.int_to_bit_arr(6, 3) == np.array([0, 1, 1])).all()
    assert (common.int_to_bit_arr(7, 3) == np.array([1, 1, 1])).all()


def test_int_to_bit_arr_invalid_args():
    # A negative integer is not allowed as an argument.
    with pytest.raises(ValueError):
        _ = common.int_to_bit_arr(-1, 5)

    with pytest.raises(ValueError):
        _ = common.int_to_bit_arr(5, -1)
