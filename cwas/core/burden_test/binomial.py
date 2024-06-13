from scipy.stats import binomtest


def binom_two_tail(n1: int, n2: int, p: float):
    n1 = int(n1)
    n2 = int(n2)
    return binomtest(k=n1, n=n1 + n2, p=p, alternative="two-sided").pvalue


def binom_one_tail(n1: int, n2: int, p: float):
    n1 = int(n1)
    n2 = int(n2)
    return binomtest(k=n1, n=n1 + n2, p=p, alternative="greater").pvalue
