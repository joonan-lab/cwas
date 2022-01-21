"""
This module includes functions for parsing files used in the CWAS
categorization step. By parsing those files, these functions make the
pandas.DataFrame objects that can be directly used in the categorization
algorithm.
"""
import pathlib
import re

import numpy as np
import pandas as pd
from cwas.core.common import int_to_bit_arr
from cwas.utils.log import print_err

# TODO: Make the code much clearer
def parse_annotated_vcf(vcf_path: pathlib.Path) -> pd.DataFrame:
    """ Parse a Variant Calling File (VCF) that includes Variant Effect
    Predictor (VEP) and CWAS annotation information and make a
    pandas.DataFrame object listing annotated variants.
    """
    variant_col_names = []
    variant_rows = []  # Item: a list of values of each column
    csq_field_names = []  # CSQ is information from VEP
    annot_field_names = []  # Custom annotation field names

    # Read and parse the input VCF
    with vcf_path.open("r") as vep_vcf_file:
        has_col_name = False
        has_csq_info = False
        has_annot_info = False
        for line in vep_vcf_file:
            if line.startswith("#"):  # Comments
                if line.startswith("#CHROM"):
                    has_col_name = True
                    variant_col_names = line[1:].rstrip("\n").split("\t")
                elif line.startswith("##INFO=<ID=CSQ"):
                    has_csq_info = True
                    csq_line = line.rstrip('">\n')
                    info_format_start_idx = re.search(
                        r"Format: ", csq_line
                    ).span()[1]
                    csq_field_names = csq_line[info_format_start_idx:].split(
                        "|"
                    )
                elif line.startswith("##INFO=<ID=ANNOT"):
                    has_annot_info = True
                    annot_line = line.rstrip('">\n')
                    annot_field_str_idx = re.search(r"Key=", annot_line).span()[
                        1
                    ]
                    annot_field_names = annot_line[annot_field_str_idx:].split(
                        "|"
                    )
            else:  # Rows of variant information follow the comments.
                assert has_col_name, "The VCF does not have column names."
                assert has_csq_info, "The VCF does not have CSQ information."
                assert has_annot_info, (
                    "The VCF does not have annotation " "information."
                )
                variant_row = line.rstrip("\n").split("\t")
                variant_rows.append(variant_row)

    result = pd.DataFrame(variant_rows, columns=variant_col_names)
    try:
        info_df = _parse_info_column(
            result["INFO"], csq_field_names, annot_field_names
        )
    except KeyError:
        print_err(
            "The VCF does not have INFO column or "
            "the INFO values do not have expected field keys."
        )
        raise

    result.drop(columns="INFO", inplace=True)
    result = pd.concat([result, info_df], axis="columns")

    return result


def _parse_info_column(
    info_column: pd.Series, csq_field_names: list, annot_field_names: list
) -> pd.DataFrame:
    """ Parse the INFO column and make a pd.DataFrame object """
    info_values = info_column.values
    info_dicts = list(map(_parse_info_str, info_values))
    info_df = pd.DataFrame(info_dicts)
    csq_df = _parse_csq_column(info_df["CSQ"], csq_field_names)
    annot_df = _parse_annot_column(info_df["ANNOT"], annot_field_names)
    info_df.drop(columns=["CSQ", "ANNOT"], inplace=True)
    info_df = pd.concat([info_df, csq_df, annot_df], axis="columns")

    return info_df


def _parse_info_str(info_str: str) -> dict:
    """ Parse the string of the INFO field to make a dictionary """
    info_dict = {}
    key_value_pairs = info_str.split(";")

    for key_value_pair in key_value_pairs:
        key, value = key_value_pair.split("=", 1)
        info_dict[key] = value

    return info_dict


def _parse_csq_column(
    csq_column: pd.Series, csq_field_names: list
) -> pd.DataFrame:
    """ Parse the CSQ strings in the CSQ column and make a pd.DataFrame
    object
    """
    csq_values = csq_column.values
    csq_records = list(map(lambda csq_str: csq_str.split("|"), csq_values))
    csq_df = pd.DataFrame(csq_records, columns=csq_field_names)

    return csq_df


def _parse_annot_column(
    annot_column: pd.Series, annot_field_names: list
) -> pd.DataFrame:
    """ Parse the annotation integer in the ANNOT column and make a
    pd.DataFrame object
    """
    annot_ints = annot_column.values.astype(int)
    annot_field_cnt = len(annot_field_names)
    annot_records = list(
        map(
            lambda annot_int: int_to_bit_arr(annot_int, annot_field_cnt),
            annot_ints,
        )
    )
    annot_df = pd.DataFrame(annot_records, columns=annot_field_names)

    return annot_df


def parse_gene_matrix(gene_matrix_path: pathlib.Path) -> dict:
    """ Parse the gene matrix file and make a dictionary.
    The keys and values of the dictionary are gene symbols
    and a set of type names where the gene is associated,
    respectively.
    """
    gene_matrix_dict = {}

    with gene_matrix_path.open("r") as gene_matrix_file:
        header = gene_matrix_file.readline()
        all_gene_type_names = np.array(header.rstrip("\n").split("\t")[1:])

        for line in gene_matrix_file:
            fields = line.rstrip("\n").split("\t")
            gene_symbol = fields[0]
            gene_matrix_values = (
                np.array(fields[1:]) == "1"
            )  # Convert to the boolean array
            gene_type_names = all_gene_type_names[gene_matrix_values]
            gene_matrix_dict[gene_symbol] = set(gene_type_names)

    return gene_matrix_dict
