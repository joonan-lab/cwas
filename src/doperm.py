#!/usr/bin/env python
__version__ = 0.2
__author__ = 'Joon An'
__date__ = 'October 5th, 2018'

description = '''
                Script for burden permutation.
            '''

import os,sys,argparse,glob,gzip
import numpy as np
import pandas as pd
import pickle

import multiprocessing as mp
from multiprocessing import Pool
from functools import partial
import pyximport; pyximport.install(language_level=3, setup_args={'include_dirs': np.get_include()})
import perm as ctest


def main(mode, infile, burden_file, adj_file, trim_file, swap_file, output_tag, number_threads, cats_start, cats_end, s3_path, family_number):
    ## Get swap index information
    if mode == 'index':
        print('[Progress] The option for creating family swap index is given')
        list_idx = ctest.create_index(family_number)

        ## pickle
        pickle.dump(list_idx, open(swap_file, 'wb'))
        print('[Progress] Done creating family swap index')
        sys.exit(0)

    elif mode == 'trim':
        ## Load sumvar data
        option_compression = 'gzip' if '.gz' in infile else None
        df_raw = pd.io.parsers.read_csv(infile, sep='\t', index_col=False, compression=option_compression)

        ## Remove redundant categories
        rdd_cats = open(trim_file).read().splitlines()
        df_raw = df_raw[df_raw.columns[~df_raw.columns.isin(rdd_cats)]]
        print(df_raw.shape)

        ## Write an output
        outfile = infile.replace('txt','trimmed.txt')
        df_raw.to_csv(outfile, sep='\t', index=False)
        print('[Progress] Done trimming redundant categories')
        sys.exit(0)

    elif mode == 'perm':
        ## Load family swap index
        print('[Progress] Loading a file for family swap index. File: %s' % swap_file)
        list_idx = pickle.load(open(swap_file, 'rb'))
        print('[Progress] Loaded family swap index. Index contains %s lists and %s families ' % (str(len(list_idx)), str(len(list_idx[0]))))

    elif mode == 'merge':
        ## merge p-values
        outfile_perm_p = '.'.join(['result','perm_p', output_tag, 'txt.gz'])
        outfile_burdenshift = '.'.join(['result','perm_burdenshift', output_tag, 'txt.gz'])
        o = gzip.open(outfile_burdenshift, 'w')
        fs = sorted(glob.glob('perm_p.*.gz'))
        list_perm_p1 = [['Annotation_combo', 'Perm_p']]
        header = ['Annotation_combo'] + [str(n) for n in range(1,10001)]
        o.write('\t'.join(header) + '\n')
        for f in fs:
            cat = f.split('.')[1]
            fh = gzip.open(f).read().splitlines()
            perm_p1 = [cat, fh[0]] # perm p values / permutations
            list_perm_p1.append(perm_p1)
            perm_p2 = [cat] + fh[1:10001] # all permutation p values
            o.write('\t'.join(perm_p2) + '\n')

        o.close()

        ## concat the burden and perm matrix
        df_perm = pd.DataFrame(list_perm_p1[1:], columns=list_perm_p1[0])
        option_compression = 'gzip' if '.gz' in burden_file else None
        df_burden = pd.io.parsers.read_csv(burden_file, sep='\t', index_col=False, compression=option_compression)
        df_burden = df_burden[df_burden['Annotation_combo'].isin(df_perm['Annotation_combo'].tolist())]
        pd.merge(df_burden, df_perm, how='inner', on='Annotation_combo').to_csv(outfile_perm_p, sep='\t', index=False, compression='gzip')

        ## merge rr
        outfile_perm_rr = '.'.join(['result','perm_rr', output_tag, 'txt.gz'])
        o = gzip.open(outfile_perm_rr, 'w')
        o.write('\t'.join(header) + '\n')

        fs = sorted(glob.glob('perm_rr.*.gz'))
        for f in fs:
            print(fs.index(f))
            cat = f.split('.')[1]
            fh = gzip.open(f).read().splitlines()
            perm_rr = [cat] + fh[0:10000] # all permutation rr
            o.write('\t'.join(perm_rr) + '\n')
        o.close()

        ## merge counts
        outfile_perm_count_pro = '.'.join(['result','perm_count_pro', output_tag, 'txt.gz'])
        outfile_perm_count_sib = '.'.join(['result','perm_count_sib', output_tag, 'txt.gz'])
        o_pro = gzip.open(outfile_perm_count_pro, 'w')
        o_sib = gzip.open(outfile_perm_count_sib, 'w')
        o_pro.write('\t'.join(header) + '\n')
        o_sib.write('\t'.join(header) + '\n')

        fs = sorted(glob.glob('perm_count.*.gz'))
        for f in fs:
            print(fs.index(f))
            cat = f.split('.')[1]
            fh = gzip.open(f).read().splitlines()
            perm_pro = [cat]
            perm_sib = [cat]
            for l in fh[0:10000]:
                perm_pro.append(l.split('\t')[0])
                perm_sib.append(l.split('\t')[0])
            o_pro.write('\t'.join(perm_pro) + '\n')
            o_sib.write('\t'.join(perm_pro) + '\n')
        o_pro.close()
        o_sib.close()

        ## Send file to s3
        if s3_path != 'no':
            print('[Progress] Copy results to s3 %s' % s3_path)
            cmd = ' '.join(['for file in result.perm_*gz; do aws s3 cp $file', s3_path, '; done'])
            os.system(cmd)
        else:
            print('[Progress] No option is given for s3')

        sys.exit(0)

    else:
        print('[Error] Please specify your mode (-m)')
        sys.exit(0)

    ## Load sumvar data
    option_usecols = ''
    if cats_start == 'no':
        option_usecols = 'None'
    else:
        option_usecols = [0] + list(range(cats_start, cats_end))
        print(option_usecols)
    option_compression = 'gzip' if '.gz' in infile else None
    df_raw = pd.io.parsers.read_csv(infile, sep='\t', index_col=False, compression=option_compression, usecols=option_usecols)

    ## Load burden data (without permutation)
    option_compression = 'gzip' if '.gz' in burden_file else None
    df_burden = pd.io.parsers.read_csv(burden_file, sep='\t', index_col=False, compression=option_compression)

    ## Check the number of families between sumvar and family swap index
    if len(df_raw.SampleID.unique()) == (len(list_idx[1]) * 2):
        print('[Progress] The number of families is matching between sumvar and family swap index')
    else:
        print('[Progress] The number of families is not matching between sumvar and family swap index')
        sys.exit(0)

    ## Adjust the rate of de novo variants by covariates
    nSamples = len(df_raw.SampleID.unique())
    print('[Progress] Total %s samples in this dataset' % str(nSamples))
    print('[Progress] Adjust the rate of de novo variants based on %s' % adj_file)

    ## Get the number of categories
    ncols = len(df_raw.columns)

    ## Load information for adjustment
    adj_info = pd.io.parsers.read_csv(adj_file, sep='\t', index_col=False)

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
    del df_raw
    cols = df_adj.columns.tolist()
    cols = cols[-1:] + cols[:-1]
    df_adj = df_adj.loc[:, cols]
    print('[Progress] Adjustment for the DNV rate')

    ## Add Fam and Role
    print('[Progress] Update information for family and role')
    df_adj[['Fam', 'Role']] = df_adj['SampleID'].str.split('_', expand=True)
    df_adj = df_adj.drop('SampleID', 1)
    df_adj.loc[ df_adj.Role.isin(['s2','s3']), 'Role'] = 's1'

    ## Do Permutation!
    print('[Progress] Start permutation')
    cats = df_adj.columns.tolist()
    cats.remove('Fam')
    cats.remove('Role')
    df_cats = [df_adj[[c,'Fam','Role']] for c in df_adj[cats].columns]
    print('[Progress] Total %s categories to be permuted' % str(len(cats)))
    pool = mp.Pool(number_threads)
    pool.map_async(partial(ctest.doperm, df_burden=df_burden, swap_index=list_idx), df_cats)
    pool.close()
    pool.join()
    print('[Progress] Completed permutation')

    ## Send file to s3
    if s3_path != 'no':
        print('[Progress] Copy results to s3 %s' % s3_path)
        comp = '.'.join(['set_perm',str(cats_start),'tar.gz'])
        os.system(' '.join(['tar', '-czvf', comp, 'perm*gz']))
        cmd = ' '.join(['for file in set_perm*gz; do aws s3 cp $file', s3_path, '; done'])
        os.system(cmd)
        os.system('rm perm*gz set_perm*')
    else:
        print('[Progress] No option is given for s3')


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=description)

    ## Arguments for perm
    parser.add_argument('-m','--mode', required=True, type=str, choices=['perm', 'index', 'trim', 'merge'],
                        help='Please choose a mode to do permutations (perm) or create family swap index (index)')
    parser.add_argument('-i','--infile', type=str, help='Input File', default='no')
    parser.add_argument('-b','--burden_file', type=str, help='Non-permutaiton burden matrix file', default='no')
    parser.add_argument('-a','--adj_file', required=False, type=str, help='File to adjust the DNV rate for covariates', default='no')
    parser.add_argument('-r','--trim_file', required=False, type=str, help='File to remove redundant categories', default='no')
    parser.add_argument('-o','--output_tag', required=False, type=str, help='Output tag', default='output')
    parser.add_argument('-t','--number_threads', required=False, type=int, help='Number of threads', default=4)
    parser.add_argument('-s','--swap_file', type=str, help='File for family swap index', default='no')
    parser.add_argument('-cats_start','--cats_start', type=int, help='Start position of categories', default=0)
    parser.add_argument('-cats_end','--cats_end', type=int, help='End position of categories', default=1)
    parser.add_argument('-s3_path','--s3_path', type=str, help='Copy path for s3', default='no')

    ## Arguments for index
    parser.add_argument('-n','--family_number', type=int, help='File to save family swap index', default=0)

    args = parser.parse_args()

    ## parse argument lists
    main(args.mode,
         args.infile, args.burden_file,
         args.adj_file, args.trim_file, args.swap_file,
         args.output_tag, args.number_threads,
         args.cats_start, args.cats_end, args.s3_path,
         args.family_number)
