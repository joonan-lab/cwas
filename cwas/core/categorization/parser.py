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
    info_df = _parse_info_column(vep_vcf_df['INFO'], csq_field_names,
                                 annot_field_names)
    vep_vcf_df.drop(columns='INFO', inplace=True)
    vep_vcf_df = pd.concat([vep_vcf_df, info_df], axis='columns')

    return vep_vcf_df


def _parse_info_column(info_column: pd.Series, csq_field_names: list,
                       annot_field_names: list) -> pd.DataFrame:
    """ Parse the INFO column and make a pd.DataFrame object """
    info_values = info_column.values
    info_dicts = list(map(_parse_info_str, info_values))
    info_df = pd.DataFrame(info_dicts)
    csq_df = _parse_csq_column(info_df['CSQ'], csq_field_names)
    annot_df = _parse_annot_column(info_df['ANNOT'], annot_field_names)
    info_df.drop(columns=['CSQ', 'ANNOT'], inplace=True)
    info_df = pd.concat([info_df, csq_df, annot_df], axis='columns')

    return info_df


def _parse_info_str(info_str: str) -> dict:
    """ Parse the string of the INFO field to make a dictionary """
    info_dict = {}
    key_value_pairs = info_str.split(';')

    for key_value_pair in key_value_pairs:
        key, value = key_value_pair.split('=', 1)
        info_dict[key] = value

    return info_dict


def _parse_csq_column(csq_column: pd.Series, csq_field_names: list) -> \
        pd.DataFrame:
    """ Parse the CSQ strings in the CSQ column and make a pd.DataFrame
    object
    """
    csq_values = csq_column.values
    csq_records = \
        list(map(lambda csq_str: csq_str.split('|'), csq_values))
    csq_df = pd.DataFrame(csq_records, columns=csq_field_names)

    return csq_df


def _parse_annot_column(annot_column: pd.Series, annot_field_names: list) -> \
        pd.DataFrame:
    """ Parse the annotation integer in the ANNOT column and make a
    pd.DataFrame object
    """
    annot_ints = annot_column.values.astype(int)
    annot_field_cnt = len(annot_field_names)
    annot_records = \
        list(map(lambda annot_int: int_to_one_hot(annot_int, annot_field_cnt),
                 annot_ints))
    annot_df = pd.DataFrame(annot_records, columns=annot_field_names)

    return annot_df
