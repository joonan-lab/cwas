"""
Tests for cwas.core.burden_test.binomial module
"""
import pytest
from scipy.stats import binomtest
from cwas.core.burden_test.binomial import binom_two_tail, binom_one_tail


class TestBinomTwoTail:
    def test_basic_two_tail(self):
        p_value = binom_two_tail(5, 5, 0.5)
        assert isinstance(p_value, float)
        assert 0 <= p_value <= 1

    def test_symmetric_case(self):
        p_value = binom_two_tail(50, 50, 0.5)
        assert p_value == 1.0

    def test_extreme_deviation(self):
        p_value = binom_two_tail(95, 5, 0.5)
        assert p_value < 0.001

    def test_zero_counts(self):
        p_value = binom_two_tail(0, 100, 0.5)
        assert 0 <= p_value <= 1

    def test_all_successes(self):
        p_value = binom_two_tail(100, 0, 0.5)
        assert 0 <= p_value <= 1

    def test_consistency_with_scipy(self):
        n1, n2, p = 20, 30, 0.4
        result = binom_two_tail(n1, n2, p)
        expected = binomtest(k=n1, n=n1 + n2, p=p, alternative="two-sided").pvalue
        assert result == expected

    def test_large_sample_size(self):
        p_value = binom_two_tail(5000, 5000, 0.5)
        assert 0 <= p_value <= 1


class TestBinomOneTail:
    def test_basic_one_tail(self):
        p_value = binom_one_tail(5, 5, 0.5)
        assert isinstance(p_value, float)
        assert 0 <= p_value <= 1

    def test_greater_than_null(self):
        p_value = binom_one_tail(70, 30, 0.5)
        assert p_value < 0.05

    def test_equal_to_expected(self):
        p_value = binom_one_tail(50, 50, 0.5)
        assert p_value >= 0.5

    def test_less_than_expected(self):
        p_value = binom_one_tail(30, 70, 0.5)
        assert p_value > 0.5

    def test_extreme_high_success(self):
        p_value = binom_one_tail(99, 1, 0.5)
        assert p_value < 0.001

    def test_consistency_with_scipy(self):
        n1, n2, p = 25, 75, 0.25
        result = binom_one_tail(n1, n2, p)
        expected = binomtest(k=n1, n=n1 + n2, p=p, alternative="greater").pvalue
        assert result == expected

    def test_comparison_two_vs_one_tail(self):
        n1, n2, p = 70, 30, 0.5
        one_tail = binom_one_tail(n1, n2, p)
        two_tail = binom_two_tail(n1, n2, p)
        assert one_tail <= two_tail
