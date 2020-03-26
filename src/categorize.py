#!/usr/bin/env python
"""
Script for categorizing variants annotated by VEP into CWAS categories, which are combinations of annotation terms,
and counting the number of variants for each category.

For more detailed information, please refer to An et al., 2018 (PMID 30545852).

"""
import argparse
import multiprocessing as mp
import os
import re
import sys
from datetime import datetime
from functools import partial

import numpy as np
import pandas as pd
import pyximport

pyximport.install(language_level=3, reload_support=True, setup_args={'include_dirs': np.get_include()})
from categorization import cwas_cat


def main():
    # Create the argument parser
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-i', '--infile', dest='in_vcf_path', required=True, type=str,
                        help='Input VCF file from VEP')
    parser.add_argument('-g', '--gene_matrix', dest='gene_mat_path', required=True, type=str,
                        help='Gene matrix file')
    parser.add_argument('-r', '--rdd_cat_file', dest='rdd_cat_path', required=False, type=str,
                        help='File that contains a list of redundant CWAS categories', default='')
    parser.add_argument('-o', '--outfile', dest='outfile_path', required=False, type=str,
                        help='Path of the output', default='cwas_cat_result.txt')
    parser.add_argument('-p', '--num_proc', dest='num_proc', required=False, type=int,
                        help='Number of processes for this script', default=1)
    parser.add_argument('-a', '--af_known', dest='af_known', required=False, type=str, choices=['yes', 'no', 'only'],
                        help='Keep the variants with known allele frequencies', default='yes')

    # Parse the arguments
    args = parser.parse_args()

    # Print the description
    print(__doc__)

    # Print and check the validity of the settings
    print(f'[Setting] The input VCF file: {args.in_vcf_path}')  # VCF from VEP
    print(f'[Setting] The gene matrix file: {args.gene_mat_path}')
    print(f'[Setting] The list of redundant CWAS categories: {args.rdd_cat_path}')
    print(f'[Setting] The output path: {args.outfile_path}')
    print(f'[Setting] No. processes for this script: {args.num_proc:,d}')
    assert os.path.isfile(args.in_vcf_path), f'The input VCF file "{args.in_vcf_path}" cannot be found.'
    assert os.path.isfile(args.gene_mat_path), f'The gene matrix file: "{args.gene_mat_path}" cannot be found.'
    assert args.rdd_cat_path == '' or os.path.isfile(args.rdd_cat_path), \
        f'The list of redundant CWAS categories "{args.rdd_cat_path}" cannot be found.'
    outfile_dir = os.path.dirname(args.outfile_path)
    assert outfile_dir == '' or os.path.isdir(outfile_dir), f'The outfile directory "{outfile_dir}" cannot be found.'
    assert 1 <= args.num_proc <= mp.cpu_count(), \
        f'Invalid number of processes "{args.num_proc:,d}". It must be in the range [1, {mp.cpu_count()}].'
    print()

    print(f'[{get_curr_time()}, Progress] Load the input VCF file into the DataFrame')
    # Make the DataFrame of the annotated variants from the VCF file
    rdd_colnames = ["CHROM", "POS", "QUAL", "FILTER", "INFO", "Allele", "Allele_Rm", "IMPACT", "Gene", "Feature_type",
                    "Feature", "EXON", "INTRON", "HGVSc", "HGVSp", "cDNA_position", "CDS_position", "Protein_position",
                    "Amino_acids", "Codons", "Existing_variation", "STRAND", "FLAGS", "SYMBOL_SOURCE", "HGNC_ID",
                    "CANONICAL", "TSL", "APPRIS", "CCDS", "SOURCE", "gnomADg"]  # The list of redundant columns
    variant_df = parse_vep_vcf(args.in_vcf_path, rdd_colnames)
    print(f'[{get_curr_time()}, Progress] No. the input variants: {len(variant_df.index):,d}')

    # Create the information for the 'gene_list' annotation terms for each gene symbol
    gene_list_set_dict = parse_gene_mat(args.gene_mat_path)

    # (Optional) Filter variants by whether allele frequency is known or not in gnomAD
    if args.af_known == 'no':
        variant_df = variant_df[variant_df['gnomADg_AF'] == '']
        print(f'[{get_curr_time()}, Progress] Remove AF-known variants '
              f'(No. the remained variants: {len(variant_df.index):,d})')
    elif args.af_known == 'only':
        variant_df = variant_df[variant_df['gnomADg_AF'] != '']
        print(f'[{get_curr_time()}, Progress] Remove AF-unknown variants '
              f'(No. the remained variants: {len(variant_df.index):,d})')
    else:
        print(f'[{get_curr_time()}, Progress] Keep all the variants')

    # Split the DataFrame by SampleIDs
    print(f'[{get_curr_time()}, Progress] Split the DataFrame by SampleIDs')
    groupby_sample = variant_df.groupby('SampleID')
    sample_ids = groupby_sample.groups
    sample_var_dfs = [groupby_sample.get_group(sample_id) for sample_id in sample_ids]
    print(f'[{get_curr_time()}, Progress] No. the DataFrames (No. of the samples): {len(sample_var_dfs):,d}')

    # Categorize the variants in each sample
    print(f'[{get_curr_time()}, Progress] Categorize the variants in each samples')
    if args.num_proc == 1:
        cat_results, sample_ids = cwas_cat_samples(sample_var_dfs, gene_list_set_dict)
    else:
        try:
            var_df_sub_lists = div_list(sample_var_dfs, args.num_proc)
        except AssertionError:
            print(f'[{get_curr_time()}, ERROR] Too many number of processes "{args.num_proc:,d}". '
                  f'Please make it less than the number of samples.', file=sys.stderr)
            raise

        pool = mp.Pool(args.num_proc)
        proc_outputs = pool.map(partial(cwas_cat_samples, gene_list_set_dict=gene_list_set_dict), var_df_sub_lists)
        pool.close()
        pool.join()

        cat_results = []
        sample_ids = []

        for proc_cat_results, proc_sample_ids in proc_outputs:
            cat_results += proc_cat_results
            sample_ids += proc_sample_ids

    # Create the DataFrame for the result of the categorization
    print(f'[{get_curr_time()}, Progress] Create the DataFrame for the result of the categorization')
    cat_result_df = pd.DataFrame(cat_results).fillna(0)
    cat_result_df = cat_result_df.astype(int)
    cat_result_df['SampleID'] = sample_ids
    cat_result_df = cat_result_df.set_index('SampleID')
    print(f'[{get_curr_time()}, Progress] No. CWAS categories with at least 1 variant: '
          f'{len(cat_result_df.columns) - 1:,d}')

    # Remove redundant categories
    if args.rdd_cat_path is None:
        print(f'[{get_curr_time()}, Progress] Keep redundant categories')
    else:
        with open(args.rdd_cat_path, 'r') as rdd_cat_file:
            rdd_cats = rdd_cat_file.read().splitlines()

        cat_result_df.drop(rdd_cats, axis='columns', inplace=True, errors='ignore')  # Remove only existing columns
        print(f'[{get_curr_time()}, Progress] No. non-redundant CWAS categories with at least 1 variant: '
              f'{len(cat_result_df.columns) - 1:,d}')

    # Write the result of the categorization
    print(f'[{get_curr_time()}, Progress] Write the result of the categorization')
    cat_result_df.to_csv(args.outfile_path, sep='\t')

    print(f'[{get_curr_time()}, Progress] Done')


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


def parse_gene_mat(gene_mat_path: str) -> dict:
    """ Parse the gene matrix file and make a dictionary which key and value are a gene symbol and the set of the names
    of the gene lists where this gene is involved, respectively.

    :param gene_mat_path: The path of the gene matrix file
    :return: The dictionary that is mentioned above
    """
    gene_list_set_dict = {}

    with open(gene_mat_path, 'r') as gene_mat_file:
        header = gene_mat_file.readline()
        all_gene_list_names = np.array(header.rstrip('\n').split('\t')[1:])

        for line in gene_mat_file:
            fields = line.rstrip('\n').split('\t')
            gene_symbol = fields[0]
            in_gene_list_arr = (np.array(fields[1:]) == '1')  # Convert to the boolean array
            gene_list_names = all_gene_list_names[in_gene_list_arr]
            gene_list_set_dict[gene_symbol] = set(gene_list_names)

    return gene_list_set_dict


def cwas_cat_samples(sample_var_dfs: list, gene_list_set_dict: dict) -> (list, list):
    """ This is a wrapper function to execute 'cwas_cat' for multiple samples

    :param sample_var_dfs: The list of pd.DataFrame objects for each sample's variants
    :param gene_list_set_dict: The dictionary from 'parse_gene_mat' function
    :returns:
        1. The list of the pd.Series objects for 'cwas_cat' results of the samples
        2. The list of the sample IDs
    """
    cat_results = []  # Item: pd.Series object
    sample_ids = []

    for sample_var_df in sample_var_dfs:
        sample_id = sample_var_df['SampleID'].values[0]
        sample_ids.append(sample_id)
        cat_result_dict = cwas_cat(sample_var_df, gene_list_set_dict)
        cat_results.append(pd.Series(cat_result_dict))

    return cat_results, sample_ids


def div_list(in_list: list, n_sub: int) -> list:
    """ Divide the input list into multiple sub-lists """
    sub_lists = []
    sub_len = len(in_list) // n_sub

    if sub_len == 0:
        raise AssertionError(f'The number of sub-lists ("{n_sub:,d}") are larger than '
                             f'the length of the input list ("{len(in_list):,d}").')

    for i in range(n_sub - 1):
        sub_list = in_list[sub_len * i : sub_len * (i + 1)]
        sub_lists.append(sub_list)

    sub_lists.append(in_list[sub_len * (n_sub - 1) :])

    return sub_lists


def get_curr_time() -> str:
    now = datetime.now()
    curr_time = now.strftime('%H:%M:%S %m/%d/%y')
    return curr_time


if __name__ == "__main__":
    main()
