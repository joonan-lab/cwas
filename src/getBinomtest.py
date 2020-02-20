#!/usr/bin/env python
__version__ = 0.2
__author__ = 'Joon An'
__date__ = 'October 5th, 2018'

description = '''
				Script for burden test.
			'''

import os,sys,argparse,gzip
from os.path import expanduser
home = expanduser("~")

import pandas as pd
import numpy as np

import multiprocessing as mp
from multiprocessing import Pool, cpu_count
from functools import partial
import pyximport; pyximport.install(language_level=3, setup_args={'include_dirs': np.get_include()})
import perm as ctest


def run_binom(df_raw0, df_adj0, adj0, output_tag0, cwas):
    ## Calculate binom test for Adjusted variants.
    if adj0 == 'no':
        print('[Progress] Delete df_adj0 to keep memory')
        # this option does not need a data frame with adjustment
        del df_adj0
    else:
        pass

    ## Raw rate
    ## Collapse the number of variant by phenotype
    print('[Progress] Calculate binomial statistics for the unadjusted DNV rate')
    df_raw0 = df_raw0.drop('Fam', 1)
    df_raw0 = df_raw0.drop('SampleID', 1)
    cols = df_raw0.columns.tolist()[1:]

    df_pro = pd.DataFrame(df_raw0[df_raw0['Role'].str.contains('p1')][cols].sum())
    df_pro = df_pro.rename(index=str, columns={0: 'Case_count_raw'})
    df_sib = pd.DataFrame(df_raw0[df_raw0['Role'].str.contains('s1')][cols].sum())
    df_sib = df_sib.rename(index=str, columns={0: 'Con_count_raw'})

    burdenMat_raw = pd.concat([df_pro, df_sib], axis=1, join='inner')
    burdenMat_raw['Annotation_combo'] = burdenMat_raw.index
    burdenMat_raw['Total_count_raw'] = burdenMat_raw['Case_count_raw'] + burdenMat_raw['Con_count_raw']
    burdenMat_raw = burdenMat_raw[['Annotation_combo', 'Total_count_raw', 'Case_count_raw', 'Con_count_raw']]
    burdenMat_raw['Unadjusted_relative_risk'] = burdenMat_raw['Case_count_raw']/burdenMat_raw['Con_count_raw']

    ## Calculate binom test
    burdenMat_raw['Binom_p_raw'] = pd.Series(ctest.binomTest(burdenMat_raw.Case_count_raw.astype('int64').values,
                                                             burdenMat_raw.Con_count_raw.astype('int64').values), index=burdenMat_raw.index)
    burdenMat_raw['Binom_p_1side_raw'] = pd.Series(ctest.binomTest_onesided(burdenMat_raw.Case_count_raw.astype('int64').values,
                                                                            burdenMat_raw.Con_count_raw.astype('int64').values), index=burdenMat_raw.index)

    ## (Optional) Convert to z-score
    if cwas == 'Yes':
        print('[Progress] Convert p-values to z scores.')
        burdenMat_raw['Binom_Z_raw'] = ctest.check_pToZ(burdenMat_raw.Binom_p_1side_raw.values)

    else:
        print('[Progress] No option given for convert p-values.')

    ## Calculate binom test for Adjusted variants.
    if adj0 == 'no':
        print('[Progress] No option given to calculate binomial statistics for the adjusted DNV rate')
    else:
        print('[Progress] Calculate binomial statistics for the adjusted DNV rate')
        ## Collapse the number of variant by phenotype
        df_adj0 = df_adj0.drop('Fam', 1)
        df_adj0 = df_adj0.drop('SampleID', 1)
        cols = df_adj0.columns.tolist()[1:]
        df_pro = pd.DataFrame(df_adj0[df_adj0['Role'].str.contains('p1')][cols].sum())
        df_pro = df_pro.rename(index=str, columns={0: 'Case_count_adj'})
        df_sib = pd.DataFrame(df_adj0[df_adj0['Role'].str.contains('s1')][cols].sum())
        df_sib = df_sib.rename(index=str, columns={0: 'Con_count_adj'})

        burdenMat_adj = pd.concat([df_pro, df_sib], axis=1, join='inner')
        burdenMat_adj['Annotation_combo'] = burdenMat_adj.index
        burdenMat_adj['Total_count_adj'] = burdenMat_adj['Case_count_adj'] + burdenMat_adj['Con_count_adj']
        burdenMat_adj = burdenMat_adj[['Annotation_combo', 'Total_count_adj', 'Case_count_adj', 'Con_count_adj']]
        burdenMat_adj['Adjusted_relative_risk'] = burdenMat_adj['Case_count_adj']/burdenMat_adj['Con_count_adj']

        ## Calculate binom test
        burdenMat_adj['Binom_p'] = pd.Series(ctest.binomTest(burdenMat_adj['Case_count_adj'].round(decimals=0).astype('int64').values,
                                                             burdenMat_adj['Con_count_adj'].round(decimals=0).astype('int64').values),
                                             index=burdenMat_adj.index)
        burdenMat_adj['Binom_p_1side'] = pd.Series(ctest.binomTest_onesided(burdenMat_adj['Case_count_adj'].round(decimals=0).astype('int64').values,
                                                                            burdenMat_adj['Con_count_adj'].round(decimals=0).astype('int64').values),
                                                   index=burdenMat_adj.index)

        ## Rounding numbers
        burdenMat_adj['Total_count_adj'] = burdenMat_adj['Total_count_adj'].round(decimals=2)
        burdenMat_adj['Case_count_adj'] = burdenMat_adj['Case_count_adj'].round(decimals=2)
        burdenMat_adj['Con_count_adj'] = burdenMat_adj['Con_count_adj'].round(decimals=2)


    ## Out to CSV
    print('[Progress] Save the result to csv')
    ## Merge into the cats dataframe
    outfile = '.'.join(['result','burden', output_tag0, 'txt'])
    if adj0 == 'no':
        burdenMat_raw.to_csv(outfile, sep='\t', index=False)
    else:
        burdenMat = pd.merge(burdenMat_raw, burdenMat_adj, how='inner', on='Annotation_combo')
        burdenMat.to_csv(outfile, sep='\t', index=False)

def main(infile, adj, trim_file, output_tag, cwas):
    ## Load data
    if '.gz' in infile:
        ## Set dtypes
        with gzip.open(infile, 'r') as f:
            header = f.readline()
        dtypes = {'SampleID': 'str'}
        for h in header.split('\t')[1:]:
            dtypes[h] = 'int'
        ## Open the infile
        print('[Progress] Loading data: gzip compressed file')
        df_raw = pd.io.parsers.read_csv(infile, sep='\t', index_col=False, compression='gzip', dtype=dtypes)
    else:
        ## Set dtypes
        with open(infile, 'r') as f:
            header = f.readline()
        dtypes = {'SampleID': 'str'}
        for h in header.split('\t')[1:]:
            dtypes[h] = 'int'
        ## Open the infile
        print('[Progress] Loading data: text file')
        df_raw = pd.io.parsers.read_csv(infile, sep='\t', index_col=False, dtype=dtypes)

    print('[Progress] DataFrame of input is %s' % (df_raw.shape,))

    ## Remove redundant categories
    if trim_file == 'no':
        '[Progress] The option for trimming redundant categories is not given'
    else:
        rdd_cats = open(trim_file).read().splitlines()
        df_raw = df_raw[df_raw.columns[~df_raw.columns.isin(rdd_cats)]]

    ## Adjust the rate of de novo variants by covariates
    nSamples = len(df_raw.SampleID.unique())
    print('[Progress] Total %s samples in this dataset' % str(nSamples))
    if adj == 'no':
        df_adj = 'to_be_deleted'

        ## Add Fam and Role
        df_raw[['Fam', 'Role']] = df_raw['SampleID'].str.split('_', expand=True)

        df_raw.loc[ df_raw.Role.isin(['s2','s3']), 'Role'] = 's1'

        cols = ['Role'] + df_raw.columns.tolist()[:-1]
        df_raw = df_raw.loc[:, cols]

    else:
        print('[Progress] Adjust the rate of de novo variants based on %s' % adj)
        ## Get the number of categories
        ncols = len(df_raw.columns)
        ## Load information for adjustment
        adj_info = pd.io.parsers.read_csv(adj, sep='\t', index_col=False)
        ## Check all samples in an adjustment file
        nOverlap = len(list(set(adj_info['SampleID']).intersection(set(df_raw['SampleID']))))
        if nOverlap == nSamples:
            print('[Progress] Adjustment information is given for all samples')
        else:
            nMissing = nSamples - nOverlap
            print('[Progress] Adjustment information is missing in %s sample(s)' % str(nMissing))
            sys.exit(0)
        ## Merge into the cats dataframe
        df_adj = pd.merge(df_raw, adj_info, how='inner', on='SampleID')
        ## Multiple by the rate adjustment
        df_adj = df_adj.iloc[:,1:ncols].multiply(df_adj['AdjustFactor'], axis="index")
        df_adj['SampleID'] = df_raw['SampleID']
        cols = df_adj.columns.tolist()
        cols = cols[-1:] + cols[:-1]
        df_adj = df_adj.loc[:, cols]
        # print df_adj['All_ASDTADAFDR03_All_Any_Any'].round(decimals=0).astype('int64').values

        ## Add Fam and Role
        df_raw[['Fam', 'Role']] = df_raw['SampleID'].str.split('_', expand=True)
        df_adj[['Fam', 'Role']] = df_adj['SampleID'].str.split('_', expand=True)
        df_raw.loc[ df_raw.Role.isin(['s2','s3']), 'Role'] = 's1'
        df_adj.loc[ df_adj.Role.isin(['s2','s3']), 'Role'] = 's1'
        cols = ['Role'] + df_raw.columns.tolist()[:-1]
        df_raw = df_raw.loc[:, cols]
        df_adj = df_adj.loc[:, cols]

    ## Random label assign for case/control
    if cwas == 'Yes':
        print('[Progress] CWAS option is given')
        fams_to_be_random = df_raw.Fam.unique()
        no_fams = len(fams_to_be_random) # 1902
        fams_swap = list(fams_to_be_random[np.random.binomial(1, 0.5, size=no_fams)==1])
        df_raw['Role2'] = np.where(df_raw.Fam.isin(fams_swap), np.where(df_raw.Role == 'p1', 's1', 'p1'), np.where(df_raw.Role == 'p1', 'p1', 's1'))
        df_raw['Role'] = df_raw['Role2']
        df_raw = df_raw.drop('Role2',1)

        # Save swap_families
        o_swap = open('.'.join(['swap',output_tag,'txt']), 'w')
        o_swap.write(','.join(fams_swap) + '\n')
        o_swap.close()
    else:
        pass

    ## Run burden test
    output_tag = output_tag + '.noPerm'
    run_binom(df_raw, df_adj, adj, output_tag, cwas)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('-i','--infile', required=True, type=str, help='Input File')
    parser.add_argument('-a','--adj', required=False, type=str, help='File to adjust the DNV rate for covariates', default='no')
    parser.add_argument('-r','--trim_file', required=False, type=str, help='File to remove redundant categories', default='no')
    parser.add_argument('-o','--output_tag', required=False, type=str, help='Output tag', default='output')
    parser.add_argument('-cwas','--cwas', required=False, type=str, help='Option for cwas to get zscores and random label assign for case/control', default='no')
    args = parser.parse_args()
    main(args.infile, args.adj, args.trim_file, args.output_tag, args.cwas)



