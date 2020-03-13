#!/usr/bin/env python
__version__ = 0.2
__author__ = 'Joon An'
__date__ = 'October 5th, 2018'

description = '''
				Script for categorization and summerization.
			'''

import argparse
import glob
import multiprocessing as mp
import os
import re
from functools import partial

import pandas as pd
import pyximport

pyximport.install(language_level=3, reload_support=True)
from categorization import get_col_index, parCat


def main(vep_vcf_path, gene_mat_path, num_threads, output_tag, af_known):
    # Print the run setting
    print('[Setting] Input file: %s' % vep_vcf_path)  # This file contains the list of the variants annotated by VEP.
    print('[Setting] Gene matrix file: %s' % gene_mat_path)
    print('[Setting] Number of threads: %s' % num_threads)
    print('[Setting] Output tag: %s' % output_tag)
    print('[Progress] Loading the input file into the data frame')

    # Make the DataFrame of the annotated variants from the VCF file
    rdd_colnames = ["CHROM", "POS", "QUAL", "FILTER", "INFO", "Allele", "Allele_Rm", "IMPACT", "Gene", "Feature_type",
                    "Feature", "EXON", "INTRON", "HGVSc", "HGVSp", "cDNA_position", "CDS_position", "Protein_position",
                    "Amino_acids", "Codons", "Existing_variation", "STRAND", "FLAGS", "SYMBOL_SOURCE", "HGNC_ID",
                    "CANONICAL", "TSL", "APPRIS", "CCDS", "SOURCE", "gnomADg"]  # The list of redundant columns
    variant_df = parse_vep_vcf(vep_vcf_path, rdd_colnames)

    print(f'[Progress] The number of the input variants: {len(variant_df.index)}')

    # (Optional) Filter variants by whether allele frequency is known or not in gnomAD
    if af_known == 'no':
        variant_df = variant_df[variant_df['gnomADg_AF'] == '']
        print(f'[Progress] Remove AF-known variants (The number of the remained variants: {len(variant_df.index)})')
    elif af_known == 'only':
        variant_df = variant_df[variant_df['gnomADg_AF'] != '']
        print(f'[Progress] Remove AF-unknown variants (The number of the remained variants: {len(variant_df.index)})')
    else:
        print('[Progress] Keep all the variants')

    # Get the sample information
    sample_ids = variant_df['SampleID'].unique()
    print(f'[Progress] Total {len(sample_ids)} samples are ready for analysis.')

    # Creating the header information
    header_index = get_col_index(list(variant_df.columns), gene_mat_path)
    print('[Progress] Start processing' + '\n')

    # Split dataframes by samples
    s = variant_df.groupby('SampleID')
    inputs = [s.get_group(x) for x in s.groups]
    print('[Progress] Split the dataframe by samples. Total %s dataframes' % len(inputs))

    # Creating a pool for parallel processing
    pool = mp.Pool(num_threads)
    pool.imap_unordered(partial(parCat, header_index=header_index), inputs)
    pool.close()
    pool.join()
    print('[Progress] Calculation for each samples are done' + '\n')

    # Merging files after run
    fs = sorted(glob.glob('tmp_catego*'))
    print('[Progress] Start merging %s files' % str(len(fs)))
    f = fs[0]
    fh = open(f).read().splitlines()
    out = []
    sample_id = f.replace('tmp_catego.', '').replace('.txt', '').replace('.', '_')
    match = ['SampleID'] + [a.split(';')[0] for a in fh]
    out.append(match)
    match = [sample_id] + [int(a.split(';')[1]) for a in fh]
    out.append(match)

    for f in fs[1:]:
        match = [f.replace('tmp_catego.', '').replace('.txt', '').replace('.', '_')]
        print(fs.index(f))
        with open(f) as fh:
            for l in fh:
                match.append(int(l.rstrip('\n').split(';')[1]))
        out.append(match)

    df = pd.DataFrame(out[1:], columns=out[0])

    # Writing out cat var matrix
    outfile_sumvar = '.'.join(['result', 'sumVar', output_tag, 'txt'])
    file_cats_zero = '.'.join(['result', 'cats_zero', output_tag, 'txt'])
    print('[Progress] Remove columns with zero variant')
    is_zero_df = pd.DataFrame(df == 0)
    df_zero = df.loc[:, is_zero_df.any(axis=0)]
    o = open(file_cats_zero, 'w')
    for l in list(df_zero.columns):
        o.write(l + '\n')
    o.close()
    is_nonzero_df = pd.DataFrame(df != 0)
    df = df.loc[:, is_nonzero_df.any(axis=0)]
    print('[Progress] Writing the sumvar matrix to csv files')
    df.to_csv(outfile_sumvar, sep='\t', index=False)

    # Clean up
    os.system('rm tmp*')


def parse_vep_vcf(vep_vcf_path: str, rm_colnames: list = None) -> pd.DataFrame:
    """ Parse the VCF file from VEP and make a pandas.DataFrame object for the list of the annotated variants.

    :param vep_vcf_path: The path of the VCF file that contains the list of the annotated variants by VEP
    :param rm_colnames: The list of column names that should be removed from the DataFrame.
                        (Note: the column names must exist in the original DataFrame.)
    :return: The DataFrame object that consists of the annotated variants
    """
    variant_df_rows = []
    variant_df_colnames = []
    info_field_names = []  # The list of the field names that make up the INFO field of the VCF file

    # Parse the VCF file
    with open(vep_vcf_path, 'r') as vep_vcf_file:
        for line in vep_vcf_file:
            if line.startswith('#'):  # The comments
                if line.startswith('#CHROM'):  # The header
                    variant_df_colnames = line[1:].rstrip('\n').split('\t')
                elif line.startswith('##INFO=<ID=CSQ'):
                    csq_line = line.rstrip('">\n')
                    info_format_start_idx = re.search(r'Format: ', csq_line).span()[1]
                    info_field_names = csq_line[info_format_start_idx:].split('|')
            else:
                variant_df_row = line.rstrip('\n').split('\t')
                variant_df_rows.append(variant_df_row)

    assert variant_df_rows
    assert variant_df_colnames
    assert info_field_names

    variant_df = pd.DataFrame(variant_df_rows, columns=variant_df_colnames)

    # Expand the INFO column
    variant_df[info_field_names] = variant_df['INFO'].str.split('|', expand=True)

    # Expand the Allele column (which is one field of the INFO)
    allele_field_names = ['SampleID', 'Batch', 'Allele_Rm']
    variant_df[allele_field_names] = variant_df['Allele'].str.split(';', expand=True)

    if rm_colnames is not None:
        variant_df.drop(columns=rm_colnames, inplace=True)

    return variant_df


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('-i', '--infile', required=True, type=str, help='Input File')
    parser.add_argument('-g', '--gene_matrix', required=False, type=str, help='Gene matrix File',
                        default='geneMatrixV38_v1.txt')
    parser.add_argument('-t', '--number_threads', required=False, type=int, help='Number of threads', default=1)
    parser.add_argument('-o', '--output_tag', required=False, type=str, help='Output tag', default='output')
    parser.add_argument('-a', '--af_known', required=False, type=str, help='Keep known variants', default='yes')

    # Arguments that are not used yet
    parser.add_argument('-lof', '--lof', required=False, type=str, help='Keep lof variants', default='No')

    args = parser.parse_args()
    main(args.infile, args.gene_matrix, args.number_threads, args.output_tag, args.af_known)
