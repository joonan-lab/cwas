import numpy as np
from scipy.stats import norm

from cwas.burden_test import BurdenTest
from cwas.core.burden_test.fisher_exact import fisher_two_tail, fisher_one_tail
from cwas.utils.log import print_progress


class FisherExactTest(BurdenTest):
    def run_burden_test(self):
        print_progress("Run Fisher's exact test")
        
        if not self.use_n_carrier_per_category:
            raise RuntimeError(
                "This method is only for '-n' option."
            )
        
        self._result["OR"], self._result["P"] = np.vectorize(fisher_two_tail)(
            self.case_carrier_cnt,
            self.ctrl_carrier_cnt,
            self.case_cnt - self.case_carrier_cnt,
            self.ctrl_cnt - self.ctrl_carrier_cnt,
        )

        _, self._result["P_1side"] = np.vectorize(fisher_one_tail)(
            self.case_carrier_cnt,
            self.ctrl_carrier_cnt,
            self.case_cnt - self.case_carrier_cnt,
            self.ctrl_cnt - self.ctrl_carrier_cnt,
        )