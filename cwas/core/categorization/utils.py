def get_idx_dict(list_: list) -> dict:
    return {item: i for i, item in enumerate(list_)}


def extract_sublist_by_int(list_: list, n: int) -> list:
    """Get a sublist from the input list by using the input integer"""
    i = 0
    result = []

    while n != 0:
        if n % 2 == 1:
            result.append(list_[i])
        n >>= 1
        i += 1

    return result
