import argparse, os
from pathlib import Path
import numpy as np
import pandas as pd
from tqdm import tqdm
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.ticker as mticker
import polars as pl
import re

from cwas.utils.log import print_progress, print_arg
from cwas.runnable import Runnable
from cwas.utils.check import check_is_file, check_is_dir

pd.set_option('mode.chained_assignment',  None)

class BurdenShift(Runnable):
    def __init__(self, args: argparse.Namespace):
        super().__init__(args)
        self._burden_res = None
        self._burden_shift_res = None
        self._cat_sets = None
        self._cat_counts = None
        self._cat_set_list = None
        
    @staticmethod
    def _print_args(args: argparse.Namespace):
        print_arg("Input file (burden test)", os.path.basename(args.input_path))
        print_arg("Burden shift result file", os.path.basename(args.burden_res))
        print_arg("Output directory", args.output_dir_path)
        print_arg("Category set file", os.path.basename(args.cat_set_file))
        print_arg("Category counts file", os.path.basename(args.cat_count_file))
        print_arg("Cutoff of category sets and p-value", [args.count_cutoff, args.pval])
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
    def plot_title(self) -> float:
        return self.args.plot_title
    
    @property
    def burden_res(self):
        if self._burden_res is None:
            self._burden_res = pl.read_csv(self.input_file, separator="\t")
        return self._burden_res
    
    @property
    def burden_shift_res(self) -> Path:
        if self._burden_shift_res is None:
            self._burden_shift_res = pl.read_csv(self.args.burden_res.resolve(), separator="\t")
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
            self._cat_counts = pd.read_csv(self.args.cat_count_file.resolve(), sep='\t')
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
    def pval(self) -> float:
        return self.args.pval
    
    @property
    def cat_set_list(self):
        if not self.args.cat_set_list is None:
            if self._cat_set_list is None:
                with open(self.args.cat_set_list, 'r') as f:
                    self._cat_set_list = f.read().splitlines()
        return self._cat_set_list

    @property
    def n_cat_sets(self) -> int:
        return self.args.n_cat_sets
    
    @property
    def fontsize(self) -> int:
        return self.args.fontsize
    

    def _create_category_sets(self):
        print_progress("Create category sets combined all of GENCODE, regions, and gene list")
        catsets = pd.read_csv(self.cat_set_file, sep="\t")
        catsets_dict = catsets.to_dict('list')        
        genesets = sorted(list(set(catsets['gene_set'].unique()) - set(["Any"])))
        gencodes = ['coding','noncoding','promoter','UTR','intergenic','intron','lincRNA']

        ## all gencodes
        for b in gencodes:
            col = f'is_{b}'
            catsets_dict[col] = list((catsets[f'is_{b}'] == 1).astype(int))

        ## for all genesets
        for g in genesets:
            col = f'is_{g}'
            catsets_dict[col] = list((catsets['gene_set'] == g).astype(int))
            
            # all genesets & all gencodes
            for b in gencodes:
                #if b != 'coding':
                col = f'is_{b}_{g}'
                catsets_dict[col] = list(((catsets['gene_set']==g)&(catsets[f'is_{b}']==1)).astype(int))
            
        ## for all functional_annotations
        regions = sorted(list(set(catsets['functional_annotation'].unique()) - set(["Any"])))
        for r in regions:
            col = f'is_{r}'
            catsets_dict[col] = list((catsets['functional_annotation'] == r).astype(int))
            
            ## all regions
            for b in gencodes:
                #if 'lincRNA' != b:
                col = f'is_{b}_{r}'
                catsets_dict[col] = list(((catsets['functional_annotation'] == r)&(catsets[f'is_{b}'] == 1)).astype(int))

            ## all regions & all genesets
            for g in genesets:
                col = f'is_{r}_{g}'
                catsets_dict[col] = list(((catsets['functional_annotation'] == r)&(catsets['gene_set'] == g)).astype(int))

        ## for all conserved
        functional_score = sorted(list(set(catsets['functional_score'].unique()) - set(["All"])))
        for c in functional_score:
            col = f'is_{c}'
            catsets_dict[col] = list((catsets['functional_score'] == c).astype(int))

            ## all conserved & all gencodes
            for b in gencodes:
                #if 'lincRNA' != b:
                col = f'is_{b}_{c}'
                catsets_dict[col] = list(((catsets['functional_score'] == c)&(catsets[f'is_{b}'] == 1)).astype(int))

            ## all conserved & all genesets
            for g in genesets:
                col = f'is_{r}_{g}'
                catsets_dict[col] = list(((catsets['functional_score'] == c)&(catsets['gene_set'] == g)).astype(int))

        for key in ['variant_type', 'gene_set', 'functional_score', 'gencode', 'functional_annotation']:
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

        ft_size = self.fontsize
        larger = max(nObsCase, nObsCtrl)
        max_cnts = round(ggpermCounts['N_signif_tests'].max())
        xticks = np.arange(max_cnts/4, (max_cnts/4)*4+1, max_cnts/4) if max_cnts > 0 else [max_cnts]
        
        plot_title = setName.replace('is_', '').capitalize()
        subtitle = f'No. cats in Cases: {nObsCase}, No. cats in Controls: {nObsCtrl}'
        
        fig1 = plt.figure(figsize=(5,4))
        plt.title(plot_title, fontsize=ft_size, loc='left', weight='bold')
        plt.suptitle(subtitle, fontsize=ft_size, y=0.92)
        ax = sns.kdeplot(ggpermCounts['N_signif_tests'], edgecolor='black', color='#EBEBEB',
                         alpha=1, fill=True, linewidth=1, warn_singular=False)
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

        df2 = pd.melt(df, var_name='type')
        fig2 = plt.figure(figsize=(5,4))
        plt.title(plot_title, fontsize=ft_size, loc='left', weight='bold')
        plt.suptitle(subtitle, fontsize=ft_size, y=0.92)
        ax = sns.kdeplot(df2.loc[df2['type']=='case','value'], edgecolor='black', color='#ff8a89', 
                        alpha=.6, fill=True, linewidth=1, label='Case', warn_singular=False)
        ax = sns.kdeplot(df2.loc[df2['type']=='control','value'], edgecolor='black', color='#8b8aff', 
                        alpha=.6, fill=True, linewidth=1, label='Control', warn_singular=False)
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
    
    def _burden_shift_size(self, x, bins):
        if x <= bins[0]:
            return 1
        elif bins[0] < x <= bins[1]:
            return 3
        elif bins[1] < x <= bins[2]:
            return 5
        elif bins[2] < x <= bins[3]:
            return 7
        elif bins[3] < x <= bins[4]:
            return 9
        elif bins[4] < x <= bins[5]:
            return 11
        else:
            return 13

    
    def _change_cre_name(self, x):
        pat = re.search(r'([a-zA-Z]*CRE[a-zA-Z0-9]*)', x).group(1)
        return x.replace(pat, pat.replace('CRE', " ")+'-CRE')
    
    def _match_cat_sets(self, pre, x):
        if len(pre.split("&")) == len(x.replace('is_','').split("_")):
            res = set(pre.split("&"))&set(x.split("_"))
            
            if len(res) == len(pre.split("&")):
                return x
            
    def _create_shiftResPlot_df(self, df, isin_inf=False):
        main_domain = ['Coding', 'PTV', 'Missense', 'Damaging', 'Noncoding', 'Promoter', 'Intron', 'Intergenic', 'UTR', 'LincRNA']
        domain_order = ['Coding (All)','PTV','Missense','Coding w/o PTV','Noncoding (All)','Promoter','Intron','Intergenic','UTR','LincRNA','CRE','Others']
        
        ## modify dataframe as validated form
        df = df.loc[(df['Category_set'] != 'All')&(df['Category_set']!='is_lincRNA_lincRNA')]
        df['Category_set'] = df['Category_set'].str.replace('ptv','PTV')
        df['Category_set'] = df['Category_set'].str.replace('_no_','_w/o_')
        df['Category_set'] = df['Category_set'].str.replace('_wo_','_w/o_')
        df['Domain'] = df['Category_set'].str.split('_').str[1]
        df['Domain'] = df['Domain'].str[0].str.upper() + df['Domain'].str[1:]

        df_case = df.loc[:, ['Category_set','N_cats_case','P_case','Domain']]
        df_case.columns = ['Category_set','N_cats','P','Domain']
        df_case['Phenotype'] = 'Case'
        df_ctrl = df.loc[:, ['Category_set','N_cats_control','P_control','Domain']]
        df_ctrl.columns = ['Category_set','N_cats','P','Domain']
        df_ctrl['Phenotype'] = 'Control'

        df2 = pd.concat([df_case, df_ctrl])#.sort_values(['Category_set','Phenotype'])

        maxs = df2['N_cats'].max()
        mins = df2['N_cats'].min()
        if maxs >= 12:
            num_bins = 6
        elif 6 <= maxs < 12:
            num_bins = 3
        elif 1 < maxs < 6:
            num_bins = 2
        else:
            num_bins = 1
        bins = np.arange(mins, maxs + 1, (maxs - mins) // num_bins)
        df2['Size'] = df2['N_cats'].apply(lambda x: self._burden_shift_size(x, bins))

        df2.loc[~df2.Domain.isin(main_domain), 'Domain'] = 'Others'
        df2["Domain2"] = df2["Domain"]
        df2.loc[df2['Domain2'].isin(["Coding","Noncoding"]), "Domain2"] = df2.loc[df2['Domain2'].isin(["Coding","Noncoding"]), "Domain2"] + " (All)"
        df2.loc[df2['Domain2']=="Damaging", "Domain2"] = "Missense"
        df2.loc[df2['Category_set']=="is_coding_w/o_promoter", "Domain2"] = "Coding w/o PTV"
        df2['Domain_order'] = df2.Domain2.apply(lambda x: domain_order.index(x))

        df2["Category_term"] = df2["Category_set"].str.replace("is_","")
        df2.loc[df2["Category_term"]=='coding_w/o_PTV', "Category_term"] = "coding w/o PTV"
        df2.loc[df2["Category_term"].str.contains("CRE"), "Category_term"] = df2.loc[df2["Category_term"].str.contains("CRE"), "Category_term"].apply(lambda x: self._change_cre_name(x))

        df2.loc[df2['Domain2']=='Others', 'Category_term'] = df2.loc[df2['Domain2']=='Others', 'Category_term'].apply(lambda x: x.replace('_', "\n"))
        df2['Category_term'] = df2["Category_term"].apply(lambda x: " ".join(x.split('_')[1:]) if len(x.split('_'))>2 else x.split('_')[-1])

        df2["Category_term"] = df2["Category_term"].str.replace("CHD8Common", "CHD8 targets")
        df2["Category_term"] = df2["Category_term"].str.replace("FMRPDarnel", "FMRP targets")
        df2["Category_term"] = df2["Category_term"].str.replace("L23", "L2/3")
        df2["Category_term"] = df2["Category_term"].str.replace("L56", "L5/6")
        df2["Category_term"] = df2["Category_term"].str.replace("WillseyUnion", "ASD coexpression")
        df2["Category_term"] = df2["Category_term"].str.replace("ASDTADAFDR03", "ASD risk")
        df2["Category_term"] = df2["Category_term"].str.replace("LOEUF37", "Constrained genes")
        df2["Category_term"] = df2["Category_term"].str.replace("Micro", "Microglia")
        df2["Category_term"] = df2["Category_term"].str.replace("Astro", "Astrocytes")
        df2['Category_term'] = df2['Category_term'].str[0].str.upper() + df2['Category_term'].str[1:]

        df2["new_name"] = df2["Category_term"] + "::" + df2["Domain2"]
        df2["-log10P"] = df2["P"].apply(lambda x: -np.log10(x))
        
        ## extract queried data frame for drawing plot
        all_cat_sets = self.cat_sets.columns.tolist()[1:]
        query_df = pd.DataFrame()
        
        if not self.cat_set_list is None:
            query_cat_sets = []
            for x in self.cat_set_list:
                if len(x.split('&'))>1:
                    query = list(set([x for x in list(map(lambda y: self._match_cat_sets(x, y), all_cat_sets)) if not x is None]))
                    query_cat_sets.extend(query)
                else:
                    query_cat_sets.append('is_'+x)
            query_df = df2.loc[df2["Category_set"].isin(query_cat_sets)]    
        else:
            topN_cats = df2.loc[df2.Phenotype=='Case'].sort_values(["-log10P","N_cats"], ascending=False)["Category_set"].tolist()[:self.n_cat_sets]
            query_df = df2.loc[df2["Category_set"].isin(topN_cats)]

        isin_inf = True if sum(np.isinf(query_df["-log10P"]))>0 else False
        max_lim = df2.loc[~np.isinf(df2["-log10P"]), "-log10P"].max() + 0.5
        query_df.loc[np.isinf(query_df["-log10P"]), "-log10P"] = max_lim
        query_df.sort_values(["Domain_order","Category_set"], ascending=False, inplace=True)   
        
        return (query_df, isin_inf)
    

    def burden_shift(self):
        filt_cats = self.cat_counts.loc[self.cat_counts['Raw_counts'] >= self.c_cutoff, "Category"].tolist()
        filt_cat_sets = self.cat_sets.loc[self.cat_sets['Category'].isin(filt_cats)]
        print_progress("Number of category sets above category counts cutoff: {}".format(len(filt_cat_sets)))
            
        ## Loop through categories specified in the catsets object
        print_progress("Compare burden test and permutation test")
        obsTab = pd.DataFrame()
        output_name = re.sub(r'burden_test\.txt\.gz|burden_test\.txt', '', os.path.basename(self.input_file))
        plot_output = output_name + f'burdenshift_p{self.pval}_cutoff{self.c_cutoff}'
        if self.tag is not None:
            plot_output += f".{self.tag}"
        plot_output += ".dist_plot.pdf"
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
                                       "N_cats_case":nObsCase,
                                       'N_cats_control':nObsCtrl,
                                       'P_case':pCase,
                                       'P_control':pCtrl}])
                obsTab = pd.concat([obsTab, tmp_df])
                    
        pdfsave.close()
        
        self._obsTab = obsTab
        
        output_name = re.sub(r'burden_test\.txt\.gz|burden_test\.txt', '', os.path.basename(self.input_file))
        obsFile = output_name + f"burdenshift_p{self.pval}_cutoff{self.c_cutoff}"
        if self.tag is not None:
            obsFile += f".{self.tag}"
        obsFile += ".txt"
        
        obsTab.to_csv(os.path.join(self.output_dir_path, obsFile), sep="\t", index=False)
        
    def draw_shiftResPlot(self):
        output_name = re.sub(r'burden_test\.txt\.gz|burden_test\.txt', '', os.path.basename(self.input_file))
        plot_output = output_name + f'burdenshift_p{self.pval}_cutoff{self.c_cutoff}'
        if self.tag is not None:
            plot_output += f".{self.tag}"
        plot_output += ".result_plot.pdf"
        isin_inf = None
        
        plot_df, isin_inf = self._create_shiftResPlot_df(self._obsTab)
        case_df = plot_df.loc[plot_df.Phenotype=='Case'].reset_index(drop=True)
        ctrl_df = plot_df.loc[plot_df.Phenotype=='Control'].reset_index(drop=True)

        ## Draw plot
        plt.rcParams['font.size'] = self.fontsize
        h = len(plot_df)/2 * 0.7
        fig, ax = plt.subplots(ncols=2, figsize=(14, h),
                               width_ratios=[.5,13.5])

        ## main plot
        ax[1].set_title(self.plot_title, weight='bold', loc='left', pad=5)
        ax[1].scatter(case_df['-log10P'], case_df['new_name'],
                      s=case_df['Size']*20, label='case', color='#ff8a89', edgecolor='black', linewidth=0.5)
        ax[1].scatter(ctrl_df['-log10P'], ctrl_df['new_name'],
                      s=ctrl_df['Size']*20, label='contorl', color='#8b8aff', edgecolor='black', linewidth=0.5)

        y_min, y_max = ax[1].get_ylim()
        ax[1].vlines(ymin=y_min, ymax=y_max, x=-np.log10(0.05), color='red', linestyle='--', linewidth=1.5)

        x_max = ax[1].get_xlim()[-1]
        xticks = np.arange(0,x_max+.5, 0.5)
        xticklabels = [str(x) for x in xticks]
        xticklabels[-2] = 'inf' if isin_inf else xticklabels[-2]
        ax[1].set_xticks(np.arange(0, x_max+.5, 0.5))
        ax[1].set_xticklabels(xticklabels)
        ax[1].set_xlim(0-0.2, round(x_max,1)+0.2)
        ax[1].set_xlabel('P(-log10)')

        ## create and draw grouped y-label
        cat_pos = ax[1].get_yticks()
        cat_labels = [x.get_text() for x in ax[1].get_yticklabels()]
        tmp = pd.DataFrame({'y_pos':cat_pos, 'new_name':cat_labels}, index=range(len(cat_pos)))
        y_axis_df = pd.merge(tmp, plot_df.loc[:,["Category_set","Domain_order","Domain2","new_name"]], on="new_name").drop_duplicates()
        minor_ypos = pd.DataFrame(y_axis_df.groupby("Domain2").mean('y_pos'))["y_pos"].reset_index().rename({"y_pos":"minor_y_pos"}, axis=1)

        ### create vertical line for category domain name
        ax[0].set_yticks(cat_pos)
        ax[0].set_ylim(y_min, y_max)
        for d in y_axis_df["Domain2"].unique():
            y_min_val = y_axis_df.loc[y_axis_df["Domain2"]==d, "y_pos"].min()
            y_max_val = y_axis_df.loc[y_axis_df["Domain2"]==d, "y_pos"].max()
            ax[0].vlines(0, y_min_val-.25, y_max_val+0.25, color='black', linewidth=1.5, ls='-' , clip_on=False)
            
        ### add domain name (grouping y-label)
        for y_pos, y_label in zip(minor_ypos["minor_y_pos"], minor_ypos['Domain2']):
            ax[0].text(-.015, y_pos, y_label, color='black', ha="right", va='center')
        ax[0].axis('off')

        ## set yaxis of main plot
        yticklabels = [text.split("::")[0] for text in cat_labels]
        ax[1].yaxis.set_major_locator(mticker.FixedLocator(cat_pos))
        ax[1].set_yticklabels(yticklabels)
        ax[1].set_ylim(y_min, y_max)

        ## create multiple legends
        ### legend1 - phenotype
        pheno_lds = [Line2D([0],[0], markerfacecolor='#ff8a89', marker='o', markersize=self.fontsize*1.15, markeredgewidth=0.5, linewidth=0, markeredgecolor='black', label='Case'),
                     Line2D([0],[0], markerfacecolor='#8b8aff', marker='o', markersize=self.fontsize*1.15, markeredgewidth=0.5, linewidth=0, markeredgecolor='black', label='Control')]
        legend1 = ax[1].legend(handles=pheno_lds, loc='center left', bbox_to_anchor=(1,0.25), labelspacing=.7, title='Phenotype', frameon=False, borderaxespad=0.)
        legend1._legend_box.align = 'left'

        ### legend2 - size
        size_lds = []
        sizes = [1, 3, 5, 7, 9, 11, 13]
        size_labels = ['0','<10','<50','<100','<150','<200',"≥200"]
        for s, sl in zip(sizes, size_labels):
            tmp_lds = Line2D([0],[0], markerfacecolor='white', marker='o', markersize=s, markeredgewidth=0.5, markeredgecolor='black', linewidth=0, label=sl)
            size_lds.append(tmp_lds)
        legend2 = ax[1].legend(handles=size_lds, loc='center left', bbox_to_anchor=(1,0.6), labelspacing=.7, title='Number of\nsignificant\ncategories', frameon=False, borderaxespad=0.)
        legend2._legend_box.align = 'left'

        ax[1].add_artist(legend1)
        ax[1].add_artist(legend2)
        ax[1].legend(handles=[], loc='center left', bbox_to_anchor=(1.1,0.5), frameon=False)
        
        ## set ticks params
        ax[1].spines['top'].set_visible(False)
        ax[1].spines['right'].set_visible(False)
        ax[1].spines['left'].set_linewidth(1.5)
        ax[1].spines['bottom'].set_linewidth(1.5)
        ax[1].tick_params(width=1.5, size=7)

        plt.tight_layout()
        plt.savefig(os.path.join(self.output_dir_path, plot_output), bbox_inches='tight')    

    def run(self):
        self.burden_shift()
        self.draw_shiftResPlot()
        print_progress("Done")