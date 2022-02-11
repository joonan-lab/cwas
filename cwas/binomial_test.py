import argparse

import numpy as np
import pandas as pd
from scipy.stats import norm

from cwas.burden_test import BurdenTest
from cwas.core.burden_test.binomial import binom_one_tail, binom_two_tail


class BinomialTest(BurdenTest):
    def __init__(self, args: argparse.Namespace):
        super().__init__(args)
        self._phenotypes = None
        self._case_variant_cnt = None
        self._ctrl_variant_cnt = None

    @property
    def phenotypes(self) -> np.ndarray:
        if self._phenotypes is None:
            self._phenotypes = np.vectorize(
                lambda sample_id: self.sample_info.to_dict()["PHENOTYPE"][
                    sample_id
                ]
            )(self.categorization_result.index.values)
        return self._phenotypes

    @property
    def case_variant_cnt(self) -> np.ndarray:
        if self._case_variant_cnt is None:
            self._case_variant_cnt = self.categorization_result.values[
                self.phenotypes == "case", :
            ].sum(axis=0)
        return self._case_variant_cnt

    @property
    def ctrl_variant_cnt(self) -> np.ndarray:
        if self._ctrl_variant_cnt is None:
            self._ctrl_variant_cnt = self.categorization_result.values[
                self.phenotypes == "ctrl", :
            ].sum(axis=0)
        return self._ctrl_variant_cnt

    @property
    def binom_p(self) -> float:
        return (self.phenotypes == "case").sum() / len(self.phenotypes)

    def run_burden_test(self):
        # Count the number of de novo variants (DNV) for cases and controls
        dnv_cnt_arr = np.concatenate(
            [
                self.case_variant_cnt[:, np.newaxis],
                self.ctrl_variant_cnt[:, np.newaxis],
            ],
            axis=1,
        )

        # Make a DataFrame for the results of binomial tests
        burden_df = pd.DataFrame(
            dnv_cnt_arr,
            index=self.categorization_result.columns.values,
            columns=["Case_DNV_Count", "Ctrl_DNV_Count"],
        )
        burden_df.index.name = "Category"
        burden_df["Relative_Risk"] = (
            self.case_variant_cnt / self.ctrl_variant_cnt
        )
        burden_df["P"] = np.vectorize(binom_two_tail)(
            self.case_variant_cnt.round(),
            self.ctrl_variant_cnt.round(),
            self.binom_p,
        )

        # Following metrics is for getting number of effective tests
        # Add a pseudo-count in order to avoid p-values of one
        burden_df["P_1side"] = np.vectorize(binom_one_tail)(
            self.case_variant_cnt.round() + 1,
            self.ctrl_variant_cnt.round() + 1,
            self.binom_p,
        )
        burden_df["Z_1side"] = norm.ppf(
            1 - burden_df["P_1side"].values
        )  # P-value -> Z-score

        self._result = burden_df
