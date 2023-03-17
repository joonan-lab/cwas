import numpy as np
from scipy.stats import norm

from cwas.burden_test import BurdenTest
from cwas.core.burden_test.binomial import binom_one_tail, binom_two_tail
from cwas.utils.log import print_progress


class BinomialTest(BurdenTest):
    @property
    def binom_p(self) -> float:
        return (self.phenotypes == "case").sum() / len(self.phenotypes)

    def run_burden_test(self):
        print_progress("Run binomial test")
        if self.use_n_carrier_per_category:
            n1 = self.case_carrier_cnt
            n2 = self.ctrl_carrier_cnt
        else:
            n1 = self.case_variant_cnt.round()
            n2 = self.ctrl_variant_cnt.round()
        self._result["P"] = np.vectorize(binom_two_tail)(
            n1,
            n2,
            self.binom_p,
        )

        # Add the pseudocount(1) in order to avoid p-values of one
        self._result["P_1side"] = np.vectorize(binom_one_tail)(
            n1 + 1,
            n2 + 1,
            self.binom_p,
        )
        self._result["Z_1side"] = norm.ppf(1 - self._result["P_1side"].values)
