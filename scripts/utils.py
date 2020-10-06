"""
This script includes generalized functions used in analyses of this project.
"""
import os
import sys
from datetime import datetime

import numpy as np


def cmp_two_arr(array1: np.ndarray, array2: np.ndarray) -> bool:
    """ Return True if two arrays have the same items regardless of the order, else return False """
    if len(array1) != len(array2):
        return False

    array1_item_set = set(array1)

    for item in array2:
        if item not in array1_item_set:
            return False

    return True


def div_dist_num(num: int, num_group: int) -> list:
    """ Divide and distribute the number to each group almost equally. """
    num_per_groups = []
    num_per_group = num // num_group
    remain_num = num % num_group

    if num_per_group == 0:
        raise AssertionError(f'The number of groups ("{num_group:,d}") is larger then the input number("{num:,d}").')

    for _ in range(num_group):
        if remain_num == 0:
            num_per_groups.append(num_per_group)
        else:
            num_per_groups.append(num_per_group + 1)
            remain_num -= 1

    return num_per_groups


def div_list(in_list: list, num_sub_list: int) -> list:
    """ Divide the input list into multiple sub-lists """
    sub_lists = []
    sub_len = len(in_list) // num_sub_list

    if sub_len == 0:
        raise AssertionError(f'The number of sub-lists ("{num_sub_list:,d}") are larger than '
                             f'the length of the input list ("{len(in_list):,d}").')

    for i in range(num_sub_list - 1):
        sub_list = in_list[sub_len * i: sub_len * (i + 1)]
        sub_lists.append(sub_list)

    sub_lists.append(in_list[sub_len * (num_sub_list - 1):])

    return sub_lists


def get_curr_time() -> str:
    now = datetime.now()
    curr_time = now.strftime('%H:%M:%S %m/%d/%y')
    return curr_time


def swap_label(labels: np.ndarray, group_ids: np.ndarray) -> np.ndarray:
    """ Randomly swap labels (case or control) in each group and return a list of swapped labels.

    :param labels: Array of labels
    :param group_ids: Array of group IDs corresponding to each label
    :return: Swapped labels
    """
    group_to_hit_cnt = {group_id: 0 for group_id in group_ids}  # Key: A group, Value: No. times the key is referred
    group_to_idx = {}  # Key: A group, Value: The index of a label firstly matched with the group
    swap_labels = np.copy(labels)

    # Make an array for random swapping
    num_group = len(group_to_hit_cnt.keys())
    do_swaps = np.random.binomial(1, 0.5, size=num_group)
    group_idx = 0

    for i, label in enumerate(labels):
        group_id = group_ids[i]

        if group_to_hit_cnt.get(group_id, 0) == 0:
            group_to_hit_cnt[group_id] = 1
            group_to_idx[group_id] = i
        elif group_to_hit_cnt.get(group_id, 0) == 1:
            group_to_hit_cnt[group_id] += 1
            prev_idx = group_to_idx[group_id]

            # Swap
            do_swap = do_swaps[group_idx]
            if do_swap:
                swap_labels[i] = labels[prev_idx]
                swap_labels[prev_idx] = label

            group_idx += 1
        else:
            print(f'[ERROR] The group {group_id} has more than 2 labels.', file=sys.stderr)
            raise AssertionError('One group must have at most two labels.')

    return swap_labels


def execute_cmd(cmd: str):
    print(f'[{get_curr_time()}, CMD] {cmd}')
    exit_val = os.system(cmd)

    if exit_val != 0:
        print(f'[{get_curr_time()}, WARNING] This CMD is failed with this exit value {exit_val}.')
