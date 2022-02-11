import numpy as np
from scipy.stats import norm

from cwas.burden_test import BurdenTest
from cwas.core.burden_test.binomial import binom_one_tail, binom_two_tail


class BinomialTest(BurdenTest):
    @property
    def binom_p(self) -> float:
        return (self.phenotypes == "case").sum() / len(self.phenotypes)

    def run_burden_test(self):
        self._result["P"] = np.vectorize(binom_two_tail)(
            self.case_variant_cnt.round(),
            self.ctrl_variant_cnt.round(),
            self.binom_p,
        )

        # Add the pseudocount(1) in order to avoid p-values of one
        self._result["P_1side"] = np.vectorize(binom_one_tail)(
            self.case_variant_cnt.round() + 1,
            self.ctrl_variant_cnt.round() + 1,
            self.binom_p,
        )
        self._result["Z_1side"] = norm.ppf(1 - self._result["P_1side"].values)
