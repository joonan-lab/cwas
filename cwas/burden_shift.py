import argparse, os, sys
from pathlib import Path
import numpy as np
import pandas as pd
from tqdm import tqdm
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages
import time
import polars as pl


from cwas.utils.log import print_progress, print_arg
from cwas.runnable import Runnable
from scipy.stats import norm
from cwas.utils.check import check_is_file, check_is_dir

class BurdenShift(Runnable):
    def __init__(self, args: argparse.Namespace):
        super().__init__(args)
        self._burden_res = None
        self._burden_shift_res = None
        self._cat_sets = None
        self._cat_counts = None        
        
    @staticmethod
    def _print_args(args: argparse.Namespace):
        print_arg("Input file (burden test)", os.path.basename(args.input_path))
        print_arg("Burden shift result file", os.path.basename(args.burden_res))
        print_arg("Output directory", args.output_dir_path)
        print_arg("Category set file", os.path.basename(args.cat_set_file))
        print_arg("Category counts file", os.path.basename(args.cat_count_file))
        print_arg("Cutoff of p-value", args.pval)
        print_arg("Output tag (prefix of output files)", args.tag)

    @staticmethod
    def _check_args_validity(args: argparse.Namespace):
        check_is_file(args.input_path)
        check_is_file(args.burden_res)
        check_is_file(args.cat_set_file)
        check_is_file(args.cat_count_file)
        check_is_dir(args.output_dir_path)
        
    @property
    def input_file(self) -> Path:
        return self.args.input_path.resolve()
    
    @property
    def burden_res(self):
        if self._burden_res is None:
            self._burden_res = pl.read_csv(self.input_file, separator="\t")
        return self._burden_res
    
    @property
    def burden_shift_res(self) -> Path:
        if self._burden_shift_res is None:
            self._burden_shift_res = pl.read_csv(self.args.burden_res.resolve(), separator="\t")
            #self._burden_shift_res.index = self._burden_shift_res['Trial'].tolist()
        return self._burden_shift_res       
    
    @property
    def output_dir_path(self) -> Path:
        return self.args.output_dir_path.resolve()
    
    @property
    def cat_set_file(self) -> Path:
        return self.args.cat_set_file.resolve() 

    @property
    def cat_sets(self):
        if self._cat_sets is None:
            self._cat_sets = self._create_category_sets()
        return self._cat_sets
    
    @property
    def cat_counts(self) -> Path:
        if self._cat_counts is None:
            self._cat_counts = pd.read_csv(self.args.cat_count_file.resolve(), sep='\t', compression='gzip')
        return self._cat_counts
    
    @property
    def tag(self) -> str:
        return self.args.tag
    
    @property
    def c_cutoff(self) -> int:
        if self.args.count_cutoff < 0:
            raise ValueError(
                print("Count cutoff is must be positive value.")
            )        
        return self.args.count_cutoff
    
    @property
    def num_perms(self) -> int:
        return self.args.num_perms
    
    @property
    def pval(self) -> float:
        return self.args.pval
    

    def _create_category_sets(self):
        print_progress("Create category sets combined all of regions, biotypes, and gene list")
        
        catsets = pd.read_csv(self.cat_set_file, sep="\t", compression='gzip')
        catsets_dict = catsets.to_dict('list')
        
        genesets = sorted(list(set(catsets['gene_list'].unique()) - set(["Any"])))
        biotypes = ['coding','noncoding','promoter','UTR','intergenic','intron','lincRNA']
        ## for all genesets
        for g in genesets:
            col = f'is_{g}'
            catsets_dict[col] = list((catsets['gene_list'] == g).astype(int))
            
            # all genesets & all biotypes (except coding)
            for b in biotypes:
                if b != 'coding':
                    col = f'is_{b}_{g}'
                    catsets_dict[col] = list(((catsets['gene_list']==g)&(catsets[f'is_{b}']==1)).astype(int))
            
        ## for all regions
        regions = sorted(list(set(catsets['region'].unique()) - set(["Any"])))
        for r in regions:
            col = f'is_{r}'
            catsets_dict[col] = list((catsets['region'] == r).astype(int))
            
            ## all regions & all genesets
            for g in genesets:
                col = f'is_{r}_{g}'
                catsets_dict[col] = list(((catsets['region'] == r)&(catsets['gene_list'] == g)).astype(int))
            
            ## all regions & all biotypes (except lincRNA)
            for b in biotypes:
                if 'lincRNA' != b:
                    col = f'is_{b}_{r}'
                    catsets_dict[col] = list(((catsets['region'] == r)&(catsets[f'is_{b}'] == 1)).astype(int))
                    
                if 'coding' != b:
                    for g in genesets:
                        col = f'is_{b}_{r}_{g}'
                        catsets_dict[col] = list(((catsets['region'] == r)&(catsets[f'is_{b}'] == 1)&(catsets['gene_list'] == g)).astype(int))

        ## for all conserved
        conservation = sorted(list(set(catsets['conservation'].unique()) - set(["All"])))
        for c in conservation:
            col = f'is_{c}'
            catsets_dict[col] = list((catsets['conservation'] == c).astype(int))

        for key in ['variant_type', 'gene_list', 'conservation', 'gencode', 'region']:
            del catsets_dict[key]
            
        _cat_sets = pd.DataFrame(catsets_dict)
        return _cat_sets
    
    def _count_cats(self, pvals, pvalTrash):
        pvals = np.array(pvals)
        nCase = sum((pvals>0) & (pvals<=pvalTrash))
        nCtrl = sum((pvals<0) & (abs(pvals)<=pvalTrash))
        #nCase = len(pvals[(pvals>0) & (pvals<=pvalTrash)])
        #nCtrl = len(pvals[(pvals>0) & (abs(pvals)<=pvalTrash)])
        return (nCase, nCtrl)
    
    def _draw_shiftDistPlot(self, df, setName, nObsCase, nObsCtrl, pCase, pCtrl):
        ggpermCounts = pd.DataFrame(df['case'].tolist()+df['control'].tolist(), columns=['N_signif_tests'])

        ft_size = 10
        larger = max(nObsCase, nObsCtrl)
        max_cnts = round(ggpermCounts['N_signif_tests'].max())
        xticks = np.arange(max_cnts/4, (max_cnts/4)*4+1, max_cnts/4) if max_cnts > 0 else [max_cnts]
        
        fig1 = plt.figure(figsize=(5,4))
        plt.title(setName, fontsize=ft_size, loc='left', weight='bold')
        ax = sns.kdeplot(ggpermCounts['N_signif_tests'], edgecolor='black', color='#EBEBEB', alpha=1, fill=True, linewidth=1)
        ymax = ax.get_ylim()[-1]
        plt.vlines(ymin=0, ymax=ymax, x=nObsCase, color='red', linestyles='--', linewidth=1)
        plt.text(nObsCase+larger*0.02, ymax*0.8, pCase, color='red', fontsize=ft_size)
        plt.vlines(ymin=0, ymax=ymax, x=nObsCtrl, color='blue', linestyles='--', linewidth=1)
        plt.text(nObsCtrl+larger*0.02, ymax*0.9, pCtrl, color='blue', fontsize=ft_size)
        plt.xlabel('Number of significant tests', fontsize=ft_size, labelpad=5, weight='bold')
        plt.ylabel('Density', fontsize=ft_size, labelpad=5, weight='bold')
        plt.xticks(xticks, fontsize=ft_size)
        plt.yticks(fontsize=ft_size)
        plt.ylim(0, ymax)
        plt.tight_layout()
        plt.close()

        df2 = pd.melt(df, var_name='type',value_name='value')
        fig2 = plt.figure(figsize=(5,4))
        plt.title(setName, fontsize=ft_size, loc='left', weight='bold')
        ax = sns.kdeplot(df2.loc[df2['type']=='case','value'], edgecolor='black', color='#ff8a89', 
                        alpha=.6, fill=True, linewidth=1, label='Case')
        ax = sns.kdeplot(df2.loc[df2['type']=='control','value'], edgecolor='black', color='#8b8aff', 
                        alpha=.6, fill=True, linewidth=1, label='Control')
        ymax = ax.get_ylim()[-1]
        plt.vlines(ymin=0, ymax=ymax, x=nObsCase, color='red', linestyles='--', linewidth=1)
        plt.text(nObsCase+larger*0.02, ymax*0.7, pCase, color='red', fontsize=ft_size)
        plt.vlines(ymin=0, ymax=ymax, x=nObsCtrl, color='blue', linestyles='--', linewidth=1)
        plt.text(nObsCtrl+larger*0.02, ymax*0.8, pCtrl, color='blue', fontsize=ft_size)
        plt.xlabel('Number of significant tests', fontsize=ft_size, labelpad=5, weight='bold')
        plt.ylabel('Density', fontsize=ft_size, labelpad=5, weight='bold')
        plt.xticks(xticks, fontsize=ft_size)
        plt.yticks(fontsize=ft_size)
        plt.ylim(0, ymax)
        plt.legend(fontsize=7, ncol=2)
        plt.tight_layout()
        plt.close()

        return fig1, fig2        
        
    def run(self):
        self.burden_shift()
        print_progress("Done")
        
    def burden_shift(self):
        print(self.input_file, type(self.input_file))
        filt_cats = self.cat_counts.loc[self.cat_counts['Raw_counts'] >= self.c_cutoff, "Category"].tolist()
        filt_cat_sets = self.cat_sets.loc[self.cat_sets['Category'].isin(filt_cats)]
        print_progress("Number of category sets above category counts cutoff: {}".format(len(filt_cat_sets)))
            
        ## Loop through categories specified in the catsets object
        print_progress("Compare burden test and permutation test")
        
        obsTab = pd.DataFrame()
        plot_output = f'plotDistr_p{self.pval}_cutoff{self.c_cutoff}.{self.tag}.pdf'
        pdfsave = PdfPages(os.path.join(self.output_dir_path, plot_output))
        for i in tqdm(range(len(filt_cat_sets.columns))):
            setName = None
            testCats = None
            ## Define the name of the set, and the cateories within it
            if i == 0:
                setName = 'All'
                testCats = filt_cat_sets['Category']
            else:
                setName = filt_cat_sets.columns[i]
                testCats = filt_cat_sets.loc[filt_cat_sets.iloc[:,i]==1, 'Category'].tolist()    
            
            if len(testCats) > 0:
                # Subset true results to only the categories in the set
                burdenResTrim = self.burden_res.filter(self.burden_res['Category'].is_in(testCats))

                # Find case and control counts from this subset of categories (Two-sided)
                nObsCase = len(burdenResTrim.filter((burdenResTrim['P']<=self.pval)&(burdenResTrim['Relative_Risk']>1)))
                nObsCtrl = len(burdenResTrim.filter((burdenResTrim['P']<=self.pval)&(burdenResTrim['Relative_Risk']<1)))
                
                # Subset burden shift results to only the categories in the set
                burdenShitTrim = self.burden_shift_res.select(list(set(self.burden_shift_res.columns)&set(testCats))).to_pandas()
                permCounts = pd.DataFrame(burdenShitTrim.apply(lambda x: self._count_cats(x, self.pval), axis=1).tolist(), columns=['case','control'])
                                
                # Compare the observed counts to the permuted counts to calculate shift p-values
                nPermCase = len(np.where(permCounts['case']>=nObsCase)[0])
                pCase = nPermCase / len(permCounts)
                nPermCtrl = len(np.where(permCounts['control']>=nObsCtrl)[0])
                pCtrl = nPermCtrl / len(permCounts)
                                
                # Plot the results
                fig1, fig2 = self._draw_shiftDistPlot(permCounts, setName, nObsCase, nObsCtrl, pCase, pCtrl)
                pdfsave.savefig(fig1)
                pdfsave.savefig(fig2)
                
                # Add data to output object
                tmp_df = pd.DataFrame([{'Category_set':setName,
                                       "N_cats_cse":nObsCase,
                                       'N_cats_control':nObsCtrl,
                                       'P_case':pCase,
                                       'P_control':pCtrl}])
                obsTab = pd.concat([obsTab, tmp_df])
                    
        pdfsave.close()
        
        output_name = os.path.basename(self.input_file).replace('burden_test.txt.gz','')
        obsFile = output_name + f".nCats_obs_p_{self.pval}_{self.c_cutoff}.{self.tag}.txt"
        
        obsTab.to_csv(os.path.join(self.output_dir_path, obsFile), sep="\t", index=False)