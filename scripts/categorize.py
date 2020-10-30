#!/usr/bin/env python3
"""
This is a script in purpose of categorizing de novo variants (DNVs)
annotated by Variant Effect Predictor (VEP). The categories are combinations
of annotation terms, and counting the number of variants.

The output is a matrix which consists of the numbers of variants
for each category per sample.

For more detailed information, please refer to An et al., 2018 (PMID 30545852).

"""
import argparse
import multiprocessing as mp
import os
import re
import sys
from functools import partial

import numpy as np
import pandas as pd
import pyximport
import yaml

pyximport.install(
    language_level=3,
    reload_support=True,
    setup_args={'include_dirs': np.get_include()}
)
from categorization import categorize_variant as _categorize_variant
from utils import div_list, get_curr_time


def main():
    # Print the description
    print(__doc__)

    # Paths to essential configuration files
    curr_dir = os.path.dirname(os.path.abspath(__file__))
    project_dir = os.path.dirname(curr_dir)
    gene_mat_path = os.path.join(project_dir, 'conf', 'gene_matrix.txt')
    cat_conf_path = os.path.join(project_dir, 'conf', 'categories.yaml')
    rdd_cat_path = \
        os.path.join(project_dir, 'conf', 'redundant_categories.yaml')

    # Configuration file validity check
    try:
        assert os.path.isfile(gene_mat_path), f'The gene matrix file "{gene_mat_path}" cannot be found.'
        assert os.path.isfile(cat_conf_path), f'The category configuration file "{cat_conf_path}" cannot be found.'
        assert os.path.isfile(rdd_cat_path), f'The file listing redundant categories "{cat_conf_path}" cannot be found.'
    except AssertionError:
        print('[ERROR] One of default configuration files cannot be found. '
              'Please do not remove any files in the conf directory.', file=sys.stderr)
        raise

    # Parse the arguments
    parser = create_arg_parser()
    args = parser.parse_args()
    print_args(args)
    check_args_validity(args)
    print()

    # Make the DataFrame of the annotated variants from the VCF file
    print(f'[{get_curr_time()}, Progress] Load the input VCF file into a DataFrame')
    rdd_colnames = ["CHROM", "POS", "QUAL", "FILTER", "INFO", "Allele", "IMPACT", "Gene", "Feature_type",
                    "Feature", "EXON", "INTRON", "HGVSc", "HGVSp", "cDNA_position", "CDS_position", "Protein_position",
                    "Amino_acids", "Codons", "Existing_variation", "STRAND", "FLAGS", "SYMBOL_SOURCE", "HGNC_ID",
                    "CANONICAL", "TSL", "APPRIS", "CCDS", "SOURCE", "gnomADg"]  # The list of redundant columns
    variant_df = parse_vep_vcf(args.in_vcf_path, rdd_colnames)
    print(f'[{get_curr_time()}, Progress] No. input DNVs: {len(variant_df.index):,d}')

    # Create the information of the 'gene_list' annotation terms for each gene symbol
    gene_list_dict = parse_gene_mat(gene_mat_path)

    # (Optional) Filter the DNVs by whether allele frequency is known or not in gnomAD
    if args.af_known == 'no':
        variant_df = variant_df[variant_df['gnomADg_AF'] == '']
        print(f'[{get_curr_time()}, Progress] Remove AF-known variants '
              f'(No. the remained variants: {len(variant_df.index):,d})')
    elif args.af_known == 'only':
        variant_df = variant_df[variant_df['gnomADg_AF'] != '']
        print(f'[{get_curr_time()}, Progress] Remove AF-unknown variants '
              f'(No. the remained variants: {len(variant_df.index):,d})')
    else:
        print(f'[{get_curr_time()}, Progress] Keep all variants')

    # Categorize the DNVs
    print(f'[{get_curr_time()}, Progress] Categorize DNVs of each sample')
    with open(cat_conf_path, 'r') as cat_conf_file:
        category_dict = yaml.safe_load(cat_conf_file)
    try:
        cat_result_df = categorize_variant(variant_df, category_dict, gene_list_dict, args.num_proc)
    except AssertionError:
        print(f'[{get_curr_time()}, ERROR] Too many number of processes "{args.num_proc:,d}". '
              f'This number must be lower than the number of the samples.', file=sys.stderr)
        raise
    print(f'[{get_curr_time()}, Progress] No. samples: {len(cat_result_df.index.values):,d}')
    print(f'[{get_curr_time()}, Progress] No. CWAS categories with at least 1 DNV: '
          f'{len(cat_result_df.columns):,d}')

    # Remove redundant categories
    with open(rdd_cat_path, 'r') as rdd_cat_file:
        rdd_cats = yaml.safe_load(rdd_cat_file)

    cat_result_df.drop(rdd_cats, axis='columns', inplace=True, errors='ignore')  # Remove only existing columns
    print(f'[{get_curr_time()}, Progress] No. non-redundant CWAS categories with at least 1 DNV: '
          f'{len(cat_result_df.columns):,d}')

    # Write the result of the categorization
    print(f'[{get_curr_time()}, Progress] Write the result of the categorization')
    cat_result_df.to_csv(args.outfile_path, sep='\t')

    print(f'[{get_curr_time()}, Progress] Done')


def create_arg_parser() -> argparse.ArgumentParser:
    """ Create an argument parser for this script and return it """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-i', '--infile', dest='in_vcf_path', required=True, type=str,
                        help='Input VCF file from VEP')
    parser.add_argument('-o', '--outfile', dest='outfile_path', required=False, type=str,
                        help='Path of the output', default='cwas_cat_result.txt')
    parser.add_argument('-p', '--num_proc', dest='num_proc', required=False, type=int,
                        help='Number of processes for this script', default=1)
    parser.add_argument('-a', '--af_known', dest='af_known', required=False, type=str, choices=['yes', 'no', 'only'],
                        help='Keep the variants with known allele frequencies', default='yes')

    return parser


def print_args(args: argparse.Namespace):
    """ Print the settings (arguments) """
    print(f'[Setting] The input VCF file: {args.in_vcf_path}')  # VCF from VEP
    print(f'[Setting] The output path: {args.outfile_path}')
    print(f'[Setting] No. processes for this script: {args.num_proc:,d}')
    print(f'[Setting] Keep the variants with known allele frequencies: {args.af_known}')


def check_args_validity(args: argparse.Namespace):
    assert os.path.isfile(args.in_vcf_path), f'The input VCF file "{args.in_vcf_path}" cannot be found.'
    outfile_dir = os.path.dirname(args.outfile_path)
    assert outfile_dir == '' or os.path.isdir(outfile_dir), f'The outfile directory "{outfile_dir}" cannot be found.'
    assert 1 <= args.num_proc <= mp.cpu_count(), \
        f'Invalid number of processes "{args.num_proc:,d}". It must be in the range [1, {mp.cpu_count()}].'


def parse_vep_vcf(vep_vcf_path: str, rdd_colnames: list = None) -> pd.DataFrame:
    """ Parse the VCF file from VEP and make a pandas.DataFrame object listing the annotated variants.

    :param vep_vcf_path: The path of the VCF file listing annotated variants by VEP
    :param rdd_colnames: The list of column names redundant for CWAS
                         (Warning: Unavailable column names will be ignored.)
    :return: The DataFrame object listing annotated variants
    """
    variant_df_rows = []
    variant_df_colnames = []
    csq_field_names = []  # The list of the field names that make up the CSQ information (the VEP result)
    annot_field_names = []

    # Parse the VCF file
    with open(vep_vcf_path, 'r') as vep_vcf_file:
        for line in vep_vcf_file:
            if line.startswith('#'):  # The comments
                if line.startswith('#CHROM'):  # The header
                    variant_df_colnames = line[1:].rstrip('\n').split('\t')
                elif line.startswith('##INFO=<ID=CSQ'):  # A VCF from VEP must contain this line.
                    csq_line = line.rstrip('">\n')
                    info_format_start_idx = re.search(r'Format: ', csq_line).span()[1]
                    csq_field_names = csq_line[info_format_start_idx:].split('|')
                elif line.startswith('##INFO=<ID=ANNOT'):
                    annot_line = line.rstrip('">\n')
                    annot_field_str_idx = re.search(r'Key=', annot_line).span()[1]
                    annot_field_names = annot_line[annot_field_str_idx:].split('|')
            else:
                variant_df_row = line.rstrip('\n').split('\t')
                variant_df_rows.append(variant_df_row)

    vep_vcf_df = pd.DataFrame(variant_df_rows, columns=variant_df_colnames)

    # Parse the INFO field
    info_strs = vep_vcf_df['INFO'].values
    info_dicts = list(map(parse_info_str, info_strs))
    info_df = pd.DataFrame(info_dicts)

    # Parse the CSQ strings (VEP results)
    csq_strs = info_df['CSQ'].values
    csq_records = list(map(lambda csq_str: csq_str.split('|'), csq_strs))
    csq_df = pd.DataFrame(csq_records, columns=csq_field_names)

    # Parse the annotation integers
    annot_ints = info_df['ANNOT'].values.astype(int)
    annot_records = list(map(lambda annot_int: int_to_one_hot(annot_int, len(annot_field_names)), annot_ints))
    annot_df = pd.DataFrame(annot_records, columns=annot_field_names)

    # Concatenate those DataFrames
    variant_df = pd.concat([vep_vcf_df.drop(columns='INFO'), info_df.drop(columns=['CSQ', 'ANNOT']), csq_df, annot_df],
                           axis='columns')

    # Trim the columns redundant for CWAS
    if rdd_colnames is not None:
        variant_df.drop(columns=rdd_colnames, inplace=True, errors='ignore')

    return variant_df


def parse_info_str(info_str: str) -> dict:
    """ Parse the string in the INFO field of the VCF file from VEP and make a dictionary """
    info_dict = {}
    key_value_pairs = info_str.split(';')

    for key_value_pair in key_value_pairs:
        key, value = key_value_pair.split('=', 1)
        info_dict[key] = value

    return info_dict


def parse_gene_mat(gene_mat_path: str) -> dict:
    """ Parse the gene matrix file and make a dictionary which key and value are a gene symbol and the set of the names
    of the gene lists where this gene is involved, respectively.

    :param gene_mat_path: The path of the gene matrix file
    :return: The dictionary that is mentioned above
    """
    gene_list_dict = {}

    with open(gene_mat_path, 'r') as gene_mat_file:
        header = gene_mat_file.readline()
        all_gene_list_names = np.array(header.rstrip('\n').split('\t')[1:])

        for line in gene_mat_file:
            fields = line.rstrip('\n').split('\t')
            gene_symbol = fields[0]
            in_gene_list_arr = (np.array(fields[1:]) == '1')  # Convert to the boolean array
            gene_list_names = all_gene_list_names[in_gene_list_arr]
            gene_list_dict[gene_symbol] = set(gene_list_names)

    return gene_list_dict


def categorize_variant(variant_df: pd.DataFrame, category_dict: dict, gene_list_dict: dict, num_proc: int) \
        -> pd.DataFrame:
    """ Categorize the variants in the input DataFrame into CWAS categories and return DataFrame that contains
    No. variants of each CWAS category for each sample.

    :param variant_df: The DataFrame that contains a list of variants annotated by VEP
    :param category_dict: The dictionary from parsing the category configuration file
    :param gene_list_dict: The dictionary from the 'parse_gene_mat' function
    :param num_proc: No. processes used for the categorization
    :return: The DataFrame that contains No. variants of each CWAS category for each sample (Sample IDs are its indices)
    """
    # Split the DataFrame by SampleIDs
    groupby_sample = variant_df.groupby('SAMPLE')
    sample_ids = list(groupby_sample.groups)
    sample_var_dfs = [groupby_sample.get_group(sample_id) for sample_id in sample_ids]

    # Categorize the variants in each sample
    if num_proc == 1:
        cat_result_dicts = categorize_each_sample(sample_var_dfs, category_dict, gene_list_dict)
    else:
        var_df_sub_lists = div_list(sample_var_dfs, num_proc)  # It can raise AssertionError.
        pool = mp.Pool(num_proc)
        proc_outputs = \
            pool.map(partial(categorize_each_sample, category_dict=category_dict, gene_list_dict=gene_list_dict),
                     var_df_sub_lists)
        pool.close()
        pool.join()

        cat_result_dicts = []

        for proc_output in proc_outputs:
            cat_result_dicts += proc_output

    # Create the DataFrame for the result of the categorization
    cat_result_df = pd.DataFrame(cat_result_dicts).fillna(0)
    cat_result_df = cat_result_df.astype(int)
    cat_result_df['SAMPLE'] = sample_ids
    cat_result_df = cat_result_df.set_index('SAMPLE')

    return cat_result_df


def categorize_each_sample(sample_var_dfs: list, category_dict: dict, gene_list_dict: dict) -> list:
    """ This is a wrapper function to execute 'cwas_cat' for multiple samples

    :param sample_var_dfs: The list of pd.DataFrame objects listing each sample's variants
    :param category_dict: The dictionary from parsing the category configuration file
    :param gene_list_dict: The dictionary from the 'parse_gene_mat' function
    :return: The list of dictionaries for each sample's 'cwas_cat' result
    """
    cat_result_dicts = []

    for sample_var_df in sample_var_dfs:
        cat_result_dict = _categorize_variant(sample_var_df, category_dict, gene_list_dict)
        cat_result_dicts.append(cat_result_dict)

    return cat_result_dicts


def int_to_one_hot(n, one_hot_len):
    one_hot = np.zeros(one_hot_len)

    for i in range(one_hot_len):
        bit = n % 2
        one_hot[i] += bit
        n >>= 1

        if n == 0:
            break

    return one_hot


if __name__ == "__main__":
    main()
