"""
Common algorithms for CWAS
"""
import numpy as np


def cmp_two_arr(array1: np.ndarray, array2: np.ndarray) -> bool:
    """ Return True if two arrays have the same items regardless of the order.
    Otherwise, it returns False.
    """
    if len(array1) != len(array2):
        return False

    array1_item_set = set(array1)

    for item in array2:
        if item not in array1_item_set:
            return False

    return True


def div_dist_num(num: int, num_group: int) -> list:
    """ Divide and distribute the number to each group almost equally. """
    if num <= 0 or num_group <= 0:
        raise ValueError(
            'Only positive integers are accepted as arguments.'
        )

    if num < num_group:
        raise ValueError(
            'The argument "num" must be larger than the argument "num_group".'
        )

    num_per_groups = []
    num_per_group = num // num_group
    remain_num = num % num_group

    for _ in range(num_group):
        if remain_num == 0:
            num_per_groups.append(num_per_group)
        else:
            num_per_groups.append(num_per_group + 1)
            remain_num -= 1

    return num_per_groups


def div_list(in_list: list, num_sub_list: int) -> list:
    """ Divide the input list into multiple sub-lists """
    if len(in_list) == 0:
        raise ValueError(
            'This function does not accept an empty list as an argument.'
        )

    if num_sub_list <= 0:
        raise ValueError(
            'The number of sub-lists must be a positive integer.'
        )

    if len(in_list) < num_sub_list:
        raise ValueError(
            'The length of the input list must be larger '
            'than the argument "num_sub_list".'
        )

    sub_lists = []
    sub_len = len(in_list) // num_sub_list

    for i in range(num_sub_list - 1):
        sub_list = in_list[sub_len * i: sub_len * (i + 1)]
        sub_lists.append(sub_list)

    sub_lists.append(in_list[sub_len * (num_sub_list - 1):])

    return sub_lists


def swap_label(labels: np.ndarray, group_ids: np.ndarray) -> np.ndarray:
    """ Randomly swap labels (case or control) in each group
    and return a list of swapped labels.

    :param labels: Array of labels
    :param group_ids: Array of group IDs corresponding to each label
    :return: Swapped labels
    """
    # Key: a group ID, Value: No. times the key is referred
    group_to_hit_cnt = {group_id: 0 for group_id in group_ids}
    # Key: A group, Value: The index of a label firstly matched with the group
    group_to_idx = {}
    swap_labels = np.copy(labels)

    # Make an array for random swapping
    num_group = len(group_to_hit_cnt.keys())
    do_swaps = np.random.binomial(1, 0.5, size=num_group)
    group_idx = 0

    for i, label in enumerate(labels):
        group_id = group_ids[i]
        group_hit_cnt = group_to_hit_cnt.get(group_id, 0)
        assert group_hit_cnt == 0 or group_hit_cnt == 1, \
            f'Too many labels (more than 2) in a group "{group_id}".'

        if group_to_hit_cnt == 0:
            group_to_hit_cnt[group_id] = 1
            group_to_idx[group_id] = i
        else:
            group_to_hit_cnt[group_id] += 1
            prev_idx = group_to_idx[group_id]

            # Swap
            do_swap = do_swaps[group_idx]
            if do_swap:
                swap_labels[i] = labels[prev_idx]
                swap_labels[prev_idx] = label

            group_idx += 1

    return swap_labels
