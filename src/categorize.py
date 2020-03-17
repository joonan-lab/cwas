#!/usr/bin/env python
__version__ = 0.2
__author__ = 'Joon An'
__date__ = 'October 5th, 2018'

description = '''
				Script for categorization and summerization.
			'''

import argparse
import re

import numpy as np
import pandas as pd
import pyximport

pyximport.install(language_level=3, reload_support=True, setup_args={'include_dirs': np.get_include()})
from categorization import cwas_cat


def main(vep_vcf_path, gene_mat_path, num_threads, output_tag, af_known):
    # Print the run setting
    print('[Setting] Input VCF file: %s' % vep_vcf_path)  # This contains the list of the variants annotated by VEP.
    print('[Setting] Gene matrix file: %s' % gene_mat_path)
    print('[Setting] Number of threads: %s' % num_threads)
    print('[Setting] Output tag: %s' % output_tag)
    print('[Progress] Load the input VCF file into the DataFrame')

    # Make the DataFrame of the annotated variants from the VCF file
    rdd_colnames = ["CHROM", "POS", "QUAL", "FILTER", "INFO", "Allele", "Allele_Rm", "IMPACT", "Gene", "Feature_type",
                    "Feature", "EXON", "INTRON", "HGVSc", "HGVSp", "cDNA_position", "CDS_position", "Protein_position",
                    "Amino_acids", "Codons", "Existing_variation", "STRAND", "FLAGS", "SYMBOL_SOURCE", "HGNC_ID",
                    "CANONICAL", "TSL", "APPRIS", "CCDS", "SOURCE", "gnomADg"]  # The list of redundant columns
    variant_df = parse_vep_vcf(vep_vcf_path, rdd_colnames)
    print(f'[Progress] No. the input variants: {len(variant_df.index)}')

    # Create the information for the 'gene_list' annotation terms for each gene symbol
    gene_list_set_dict = parse_gene_mat(gene_mat_path)

    # (Optional) Filter variants by whether allele frequency is known or not in gnomAD
    if af_known == 'no':
        variant_df = variant_df[variant_df['gnomADg_AF'] == '']
        print(f'[Progress] Remove AF-known variants (The number of the remained variants: {len(variant_df.index)})')
    elif af_known == 'only':
        variant_df = variant_df[variant_df['gnomADg_AF'] != '']
        print(f'[Progress] Remove AF-unknown variants (The number of the remained variants: {len(variant_df.index)})')
    else:
        print('[Progress] Keep all the variants')

    # Split the DataFrame by SampleIDs
    print('[Progress] Split the DataFrame by SampleIDs')
    groupby_sample = variant_df.groupby('SampleID')
    sample_var_dfs = [groupby_sample.get_group(sample_id) for sample_id in groupby_sample.groups]
    print(f'[Progress] Total No. the DataFrames: {len(sample_var_dfs):,d}')

    # Categorize the variants in each sample
    print('[Progress] Categorize the variants in each samples' + '\n')
    cat_results = []  # Item: pd.Series object

    for sample_var_df in sample_var_dfs:
        cat_result_dict = cwas_cat(sample_var_df, gene_list_set_dict)
        cat_results.append(pd.Series(cat_result_dict))


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
