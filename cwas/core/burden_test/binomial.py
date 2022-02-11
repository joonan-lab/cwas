from scipy.stats import binom_test


def binom_two_tail(n1: int, n2: int, p: float):
    return binom_test(x=n1, n=n1 + n2, p=p, alternative="two-sided")


def binom_one_tail(n1: int, n2: int, p: float):
    return binom_test(x=n1, n=n1 + n2, p=p, alternative="greater")
