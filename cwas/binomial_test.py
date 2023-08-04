import numpy as np
from scipy.stats import norm
import matplotlib.pyplot as plt
import os

from cwas.burden_test import BurdenTest
from cwas.core.burden_test.binomial import binom_one_tail, binom_two_tail
from cwas.utils.log import print_progress


class BinomialTest(BurdenTest):
    @property
    def binom_p(self) -> float:
        return self.case_cnt / (self.case_cnt + self.ctrl_cnt)

    def run_burden_test(self):
        print_progress("Run binomial test")
        if self.use_n_carrier:
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

        self._draw_volcano_plot()

    def _draw_volcano_plot(self):
        burden_res = self._result.copy()
        burden_res['log2_RR'] = burden_res['Relative_Risk'].apply(lambda x: np.log2(x))
        burden_res['-log_P'] = burden_res['P'].apply(lambda x: -np.log10(x))
        
        threshold = -np.log10(0.05)
        max_rr = max(burden_res.loc[burden_res.log2_RR!=np.inf, 'log2_RR'])
        min_rr = min(burden_res.loc[burden_res.log2_RR!=-np.inf, 'log2_RR'])
        max_val = max(abs(max_rr),abs(min_rr))
        max_x = np.trunc(max_val) + 2

        xticks = [int(x) for x in np.arange(-max_x, max_x+1, 2)]
        xlabels = xticks.copy()
        xlabels[0] = '-Inf'
        xlabels[-1] = 'Inf'
        yticks = [int(x) for x in np.arange(0, max(burden_res['-log_P']), 2)]
        ylabels = yticks.copy()

        def replace_inf(x, v):
            if x == np.inf:
                return v
            elif x == -np.inf:
                return -v
            else:
                return x

        if self.tag != None:
            tags = self.tag.strip().split(",")
            
            for t in tags:
                fig, axes = plt.subplots(figsize=(self.plot_size, self.plot_size))

                plt.title(self.plot_title, fontsize=self.font_size, loc='left', pad=5)
                axes.vlines(x=0, ymin=-0.5, ymax=max(burden_res['-log_P'])+0.5, linestyles='-', color='lightgray', linewidth=1.25, zorder=1)
                axes.scatter(x=burden_res['log2_RR'].apply(lambda x: replace_inf(x, max_x)), y=burden_res['-log_P'],
                            marker='o', color='silver', s=self.marker_size, label='Others', edgecolor='black', linewidth=0.5, zorder=2)
                axes.scatter(x=burden_res.loc[(burden_res.index.str.contains(t))&(burden_res['-log_P']>threshold)&(burden_res['log2_RR']>0), 'log2_RR'], 
                             y=burden_res.loc[(burden_res.index.str.contains(t))&(burden_res['-log_P']>threshold)&(burden_res['log2_RR']>0),'-log_P'],
                             marker='o', label=t, facecolor='#3d62a1', s=self.marker_size, alpha=.7, edgecolor='black', linewidth=0.5, zorder=3)
                axes.hlines(y=threshold, xmin=-(max_x+1), xmax=max_x+1, linestyles='--', linewidth=1.25, color='red')
                axes.text(-(max_x+1)+0.1, threshold+0.1, 'P=0.05', size=self.font_size*0.85, color='red')
                plt.ylim(-0.5, max(burden_res['-log_P'])+0.5)
                plt.xlim(-(max_x+1), (max_x+1))
                plt.xlabel('Relative Risk ($log_{2}$)', size=self.font_size)
                plt.ylabel("P ($-log_{10}$)", size=self.font_size)
                axes.set_xticks(xticks)
                axes.set_xticklabels(xlabels, fontsize=self.font_size)
                axes.set_yticks(yticks)
                axes.set_yticklabels(ylabels, fontsize=self.font_size)
                axes.legend(markerscale=self.font_size*0.15, fontsize=self.font_size*0.85)
                axes.spines['bottom'].set_linewidth(1.25)
                axes.spines['left'].set_linewidth(1.25)
                axes.tick_params(width=1.25)
                axes.spines['top'].set_visible(False)
                axes.spines['right'].set_visible(False)
                
                output = self.result_path.name.replace(".txt", f'.{t}.volcano_plot.pdf')
                
                plt.tight_layout()
                plt.savefig(os.path.join(self.output_dir_path, output), bbox_inches='tight')
                print_progress("Save the result to the volcano plot file {}".format(output))
        else:
            fig, axes = plt.subplots(figsize=(self.plot_size, self.plot_size))

            plt.title(self.plot_title, fontsize=self.font_size, loc='left', pad=5)
            axes.vlines(x=0, ymin=-0.5, ymax=max(burden_res['-log_P'])+0.5, linestyles='-', color='lightgray', linewidth=1.25, zorder=1)
            axes.scatter(x=burden_res['log2_RR'].apply(lambda x: replace_inf(x, max_x)), y=burden_res['-log_P'],
                        marker='o', color='silver', s=self.marker_size, edgecolor='black', linewidth=0.5, zorder=2)
            axes.hlines(y=threshold, xmin=-(max_x+1), xmax=max_x+1, linestyles='--', linewidth=1.25, color='red')
            axes.text(-(max_x+1)+0.1, threshold+0.1, 'P=0.05', size=self.font_size*0.85, color='red')
            plt.ylim(-0.5, max(burden_res['-log_P'])+0.5)
            plt.xlim(-(max_x+1), (max_x+1))
            plt.xlabel('Relative Risk ($log_{2}$)', size=self.font_size)
            plt.ylabel("P ($-log_{10}$)", size=self.font_size)
            axes.set_xticks(xticks)
            axes.set_xticklabels(xlabels, fontsize=self.font_size)
            axes.set_yticks(yticks)
            axes.set_yticklabels(ylabels, fontsize=self.font_size)
            axes.spines['bottom'].set_linewidth(1.25)
            axes.spines['left'].set_linewidth(1.25)
            axes.tick_params(width=1.25)
            axes.spines['top'].set_visible(False)
            axes.spines['right'].set_visible(False)

            output = self.result_path.name.replace(".txt", f'.volcano_plot.pdf')

            plt.tight_layout()
            plt.savefig(os.path.join(self.output_dir_path, output), bbox_inches='tight')
            print_progress("Save the result to the volcano plot file {}".format(output))