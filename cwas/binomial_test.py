from typing import Tuple

import numpy as np
import pandas as pd
from scipy.stats import norm

from cwas.burden_test import BurdenTest
from cwas.core.burden_test.binomial import binom_one_tail, binom_two_tail


class BinomialTest(BurdenTest):
    def run_burden_test(self):
        # Count the number of de novo variants (DNV) for cases and controls
        cwas_cat_vals = self.categorization_result.values
        sample_info_dict = self.sample_info.to_dict()
        sample_ids = self.categorization_result.index.values
        sample_types = np.vectorize(
            lambda sample_id: sample_info_dict["PHENOTYPE"][sample_id]
        )(sample_ids)
        case_dnv_cnt, ctrl_dnv_cnt = cnt_case_ctrl_dnv(
            cwas_cat_vals, sample_types
        )
        dnv_cnt_arr = np.concatenate(
            [case_dnv_cnt[:, np.newaxis], ctrl_dnv_cnt[:, np.newaxis]], axis=1
        )

        # Make a DataFrame for the results of binomial tests
        burden_df = pd.DataFrame(
            dnv_cnt_arr,
            index=self.categorization_result.columns.values,
            columns=["Case_DNV_Count", "Ctrl_DNV_Count"],
        )
        burden_df.index.name = "Category"
        burden_df["Relative_Risk"] = case_dnv_cnt / ctrl_dnv_cnt
        burden_df["P"] = np.vectorize(binom_two_tail)(
            case_dnv_cnt.round(), ctrl_dnv_cnt.round(), 0.5
        )

        # Following metrics is for getting number of effective tests
        # Add a pseudo-count in order to avoid p-values of one
        burden_df["P_1side"] = np.vectorize(binom_one_tail)(
            case_dnv_cnt.round() + 1, ctrl_dnv_cnt.round() + 1, 0.5
        )
        burden_df["Z_1side"] = norm.ppf(
            1 - burden_df["P_1side"].values
        )  # P-value -> Z-score

        self._result = burden_df


def cnt_case_ctrl_dnv(
    sample_cat_vals: np.ndarray, sample_types: np.ndarray
) -> Tuple[float, float]:
    """ Count the number of the de novo variants for each phenotype, case and control.
    """
    are_case = sample_types == "case"
    case_dnv_cnt = sample_cat_vals[are_case, :].sum(axis=0)
    ctrl_dnv_cnt = sample_cat_vals[~are_case, :].sum(axis=0)

    return case_dnv_cnt, ctrl_dnv_cnt
