import numpy as np
import pandas as pd
from scipy.stats import norm

from cwas.burden_test import BurdenTest
from cwas.core.burden_test.binomial import binom_one_tail, binom_two_tail


class BinomialTest(BurdenTest):
    @property
    def binom_p(self) -> float:
        return (self.phenotypes == "case").sum() / len(self.phenotypes)

    def run_burden_test(self):
        variant_cnt_arr = np.concatenate(
            [
                self.case_variant_cnt[:, np.newaxis],
                self.ctrl_variant_cnt[:, np.newaxis],
            ],
            axis=1,
        )

        self._result = pd.DataFrame(
            variant_cnt_arr,
            index=self.categorization_result.columns.values,
            columns=["Case_DNV_Count", "Ctrl_DNV_Count"],
        )
        self._result.index.name = "Category"
        self._result["Relative_Risk"] = (
            self.case_variant_cnt / self.ctrl_variant_cnt
        )
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
