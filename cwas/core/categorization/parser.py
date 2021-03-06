"""
This module includes functions for parsing files used in the CWAS
categorization step. By parsing those files, these functions make the
pandas.DataFrame objects that can be directly used in the categorization
algorithm.
"""
import pathlib
import re

import pandas as pd

from cwas.core.common import int_to_one_hot


def parse_vep_vcf(vep_vcf_path: pathlib.Path) -> pd.DataFrame:
    """ Parse a Variant Calling File (VCF) from Variant Effect Predictor (VEP)
    and make a pandas.DataFrame object listing annotated variants.
    """
    variant_col_names = []
    variant_rows = []  # Item: a list of values of each column
    csq_field_names = []  # CSQ is information from VEP
    annot_field_names = []  # Custom annotation field names

    # TODO: check the format of the VCF, the input VCF must have such lines.
    # Read and parse the input VCF
    with vep_vcf_path.open('r') as vep_vcf_file:
        for line in vep_vcf_file:
            if line.startswith('#'):  # Comments
                if line.startswith('#CHROM'):
                    variant_col_names = line[1:].rstrip('\n').split('\t')
                elif line.startswith('##INFO=<ID=CSQ'):
                    csq_line = line.rstrip('">\n')
                    info_format_start_idx = \
                        re.search(r'Format: ', csq_line).span()[1]
                    csq_field_names = \
                        csq_line[info_format_start_idx:].split('|')
                elif line.startswith('##INFO=<ID=ANNOT'):
                    annot_line = line.rstrip('">\n')
                    annot_field_str_idx = \
                        re.search(r'Key=', annot_line).span()[1]
                    annot_field_names = \
                        annot_line[annot_field_str_idx:].split('|')
            else:
                variant_row = line.rstrip('\n').split('\t')
                variant_rows.append(variant_row)

    vep_vcf_df = pd.DataFrame(variant_rows, columns=variant_col_names)

    # Parse the INFO field
    info_strs = vep_vcf_df['INFO'].values
    info_dicts = list(map(_parse_info_str, info_strs))
    info_df = pd.DataFrame(info_dicts)

    # Parse the CSQ strings (VEP results)
    csq_strs = info_df['CSQ'].values
    csq_records = list(map(lambda csq_str: csq_str.split('|'), csq_strs))
    csq_df = pd.DataFrame(csq_records, columns=csq_field_names)

    # Parse the annotation integers
    annot_ints = info_df['ANNOT'].values.astype(int)
    annot_field_cnt = len(annot_field_names)
    annot_records = \
        list(map(lambda annot_int: int_to_one_hot(annot_int, annot_field_cnt),
                 annot_ints))
    annot_df = pd.DataFrame(annot_records, columns=annot_field_names)

    # Concatenate those DataFrames
    variant_df = pd.concat([
        vep_vcf_df.drop(columns='INFO'),
        info_df.drop(columns=['CSQ', 'ANNOT']),
        csq_df,
        annot_df
    ], axis='columns')

    return variant_df


def _parse_info_str(info_str: str) -> dict:
    """ Parse the string of the INFO field to make a dictionary """
    info_dict = {}
    key_value_pairs = info_str.split(';')

    for key_value_pair in key_value_pairs:
        key, value = key_value_pair.split('=', 1)
        info_dict[key] = value

    return info_dict
