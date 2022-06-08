import numpy as np
from scipy.stats import fisher_exact

def make_table(n_case_carrier: int, n_case_non_carrier: int, n_ctrl_carrier: int, n_ctrl_non_carrier: int) -> np.ndarray:
    return np.array([[n_case_carrier, n_case_non_carrier], [n_ctrl_carrier, n_ctrl_non_carrier]])

def fisher_two_tail(n_case_carrier: int, n_case_non_carrier: int, n_ctrl_carrier: int, n_ctrl_non_carrier: int):
    table = make_table(n_case_carrier, n_case_non_carrier, n_ctrl_carrier, n_ctrl_non_carrier)
    return fisher_exact(table, alternative="two-sided")


def fisher_one_tail(n_case_carrier: int, n_case_non_carrier: int, n_ctrl_carrier: int, n_ctrl_non_carrier: int):
    table = make_table(n_case_carrier, n_case_non_carrier, n_ctrl_carrier, n_ctrl_non_carrier)
    return fisher_exact(table, alternative="greater")
