from pathlib import Path
from typing import Optional
import argparse

import numpy as np
import pandas as pd
import re

from collections import defaultdict
from itertools import product

from cwas.runnable import Runnable
from cwas.utils.log import print_progress
from cwas.core.categorization.parser import (
    parse_annotated_vcf)
from cwas.utils.check import check_is_file, check_is_dir

class ExtractVariant(Runnable):
    def __init__(self, args: Optional[argparse.Namespace] = None):
        super().__init__(args)
        self._annotated_vcf = None
        self._gene_matrix = None
        self._tag = None
        self._category_set_path = None
        self._category_set = None

    @staticmethod
    def _check_args_validity(args: argparse.Namespace):
        check_is_file(args.input_path)
        check_is_dir(args.output_dir_path)
        if args.category_set_path :
            check_is_file(args.category_set_path)

    @property
    def input_path(self):
        return self.args.input_path.resolve()

    @property
    def output_dir_path(self):
        return self.args.output_dir_path.resolve()

    @property
    def annotation_info(self) -> bool:
        return self.args.annotation_info

    @property
    def tag(self) -> str:
        return self.args.tag

    @property
    def mis_info_key(self) -> str:
        return self.get_env("VEP_MIS_INFO_KEY")

    @property
    def mis_thres(self) -> float:
        return float(self.get_env("VEP_MIS_THRES"))

    @property
    def annotated_vcf(self) -> pd.DataFrame:
        if self._annotated_vcf is None:
            print_progress("Parse the annotated VCF")
            self._annotated_vcf = parse_annotated_vcf(
                Path(self.input_path)
            )
        return self._annotated_vcf

    @property
    def gene_matrix(self) -> pd.DataFrame:
        if self._gene_matrix is None:
            self._gene_matrix = pd.read_csv(self.get_env("GENE_MATRIX"), sep='\t')
            # Keep the last gene in duplicates
            self._gene_matrix = self._gene_matrix[~self._gene_matrix.duplicated(subset=['gene_name'], keep='last')]
        return self._gene_matrix

    @property
    def category_set_path(self) -> Optional[Path]:
        return (
            self.args.category_set_path.resolve()
            if self.args.category_set_path
            else None
        )
    @property
    def category_set(self) -> pd.DataFrame:
        if self._category_set is None and self.category_set_path:
            self._category_set = pd.read_csv(self.category_set_path, sep='\t')
        return self._category_set
    
    @property
    def result_path(self) -> Path:
        if self.tag is None:
            save_name = 'extracted_variants.txt.gz'
        else:
            save_name = '.'.join([self.tag, 'extracted_variants.txt.gz'])
        f_name = re.sub(r'annotated\.vcf\.gz|annotated\.vcf', save_name, self.input_path.name)
        return Path(
            f"{self.output_dir_path}/" +
            f"{f_name}"
        )
    
    def annotate_variants(self):
        print_progress("Annotate variants with annotation dataset")
        self.annotated_vcf['CLASS'] = np.where(self.annotated_vcf['REF'].str.len() > 1, 'Deletion',
                                               np.where((self.annotated_vcf['REF'].str.len() == 1) & (self.annotated_vcf['ALT'].str.len() == 1), 'SNV', 'Insertion'))
        self.annotated_vcf['DEF.GENE'] = self.annotated_vcf.apply(lambda x: x['NEAREST']
                                                                    if "downstream_gene_variant" in x['Consequence']
                                                                    or "intergenic_variant" in x['Consequence']
                                                                    else x['SYMBOL'], axis=1)
        merged_df = pd.merge(self.annotated_vcf, self.gene_matrix, left_on='DEF.GENE', right_on='gene_name', how='left')
        cols_to_fillna = self.gene_matrix.columns.tolist()
        merged_df[cols_to_fillna] = merged_df[cols_to_fillna].fillna(0)
        merged_df = merged_df.drop(columns=['gene_id', 'gene_name'])
        ## Coding
        # Define the list of string patterns to search for
        patterns = ['stop_gained', 'splice_donor', 'splice_acceptor', 'frameshift_variant', 'missense_variant', 'protein_altering_variant', 'start_lost', 'stop_lost', 'inframe_deletion', 'inframe_insertion', 'synonymous_variant', 'stop_retained_variant', 'incomplete_terminal_codon_variant', 'protein_altering_variant', 'coding_sequence_variant']
        # Use str.contains() to search for each pattern in the Consequence column
        merged_df['is_CodingRegion'] = merged_df['Consequence'].str.contains('|'.join(patterns)).astype(int)
        merged_df['is_CodingRegion'] = np.where(merged_df['ProteinCoding'] == 1, merged_df['is_CodingRegion'], 0)
        ## Noncoding
        merged_df['is_NoncodingRegion'] = np.where(merged_df['is_CodingRegion'] == 1, 0, 1)
        ## PTV
        # Define the list of string patterns to search for
        patterns = ['stop_gained', 'splice_donor', 'splice_acceptor', 'frameshift_variant']
        # Use str.contains() to search for each pattern in the Consequence column
        merged_df['is_PTVRegion'] = merged_df['Consequence'].str.contains('|'.join(patterns)).astype(int)
        merged_df['is_PTVRegion'] = np.where((merged_df['is_CodingRegion'] == 1) & (merged_df['LoF'] == 'HC') & ((merged_df['LoF_flags'] == 'SINGLE_EXON') | (merged_df['LoF_flags'] == '')), merged_df['is_PTVRegion'], 0)
        ## Frameshift
        # Define the list of string patterns to search for
        patterns = ['frameshift_variant']
        # Define the list of string patterns to exclude
        exclude_patterns = ['stop_gained', 'splice_donor', 'splice_acceptor']
        merged_df['is_FrameshiftRegion'] = ((merged_df['Consequence'].str.contains('|'.join(patterns))) 
                                            & (~merged_df['Consequence'].str.contains('|'.join(exclude_patterns)))).astype(int)
        merged_df['is_FrameshiftRegion'] = np.where(merged_df['is_PTVRegion'] == 1, merged_df['is_FrameshiftRegion'], 0)
        ## Missense
        patterns = ['missense_variant', 'protein_altering_variant', 'start_lost', 'stop_lost']
        merged_df['is_MissenseRegion'] = merged_df['Consequence'].str.contains('|'.join(patterns)).astype(int)
        merged_df['is_MissenseRegion'] = np.where((merged_df['is_CodingRegion'] == 1) & (merged_df['is_PTVRegion'] == 0), merged_df['is_MissenseRegion'], 0)
        ## Damaging missense
        merged_df['is_DamagingMissenseRegion'] = np.where((merged_df['is_MissenseRegion'] == 1) & (pd.to_numeric(merged_df["MisDb_" + self.mis_info_key], errors='coerce').fillna(0) >= self.mis_thres), 1, 0)
        # Define the list of string patterns to search for
        patterns = ['inframe_deletion', 'inframe_insertion']
        merged_df['is_InFrameRegion'] = merged_df['Consequence'].str.contains('|'.join(patterns)).astype(int)
        merged_df['is_InFrameRegion'] = np.where((merged_df['is_CodingRegion'] == 1) & (merged_df['is_PTVRegion'] == 0) & (merged_df['is_MissenseRegion'] == 0), merged_df['is_InFrameRegion'], 0)
        ## Silent
        patterns = ['synonymous_variant']
        merged_df['is_SilentRegion'] = merged_df['Consequence'].str.contains('|'.join(patterns)).astype(int)
        merged_df['is_SilentRegion'] = np.where((merged_df['is_CodingRegion'] == 1) &
                                                (merged_df['is_PTVRegion'] == 0) &
                                                (merged_df['is_MissenseRegion'] == 0) &
                                                (merged_df['is_InFrameRegion'] == 0),
                                                merged_df['is_SilentRegion'], 0)
        ## UTRs
        # Define the list of string patterns to search for
        patterns = ['_UTR_']
        merged_df['is_UTRsRegion'] = merged_df['Consequence'].str.contains('|'.join(patterns)).astype(int)
        merged_df['is_UTRsRegion'] = np.where((merged_df['is_NoncodingRegion'] == 1),
                                              merged_df['is_UTRsRegion'], 0)
        ## Promoter
        patterns = ['upstream_gene_variant']
        merged_df['is_PromoterRegion'] = merged_df['Consequence'].str.contains('|'.join(patterns)).astype(int)
        merged_df['is_PromoterRegion'] = np.where((merged_df['is_NoncodingRegion'] == 1) &
                                                  (merged_df['is_UTRsRegion'] == 0),
                                                  merged_df['is_PromoterRegion'], 0)
        ## Intron
        patterns = ['intron_variant']
        merged_df['is_IntronRegion'] = merged_df['Consequence'].str.contains('|'.join(patterns)).astype(int)
        merged_df['is_IntronRegion'] = np.where((merged_df['is_NoncodingRegion'] == 1) &
                                                  (merged_df['is_UTRsRegion'] == 0) &
                                                  (merged_df['is_PromoterRegion'] == 0),
                                                  merged_df['is_IntronRegion'], 0)
        ## SpliceSite non-canonical
        patterns = ['splice_region_variant']
        merged_df['is_SpliceSiteNoncanonRegion'] = merged_df['Consequence'].str.contains('|'.join(patterns)).astype(int)
        merged_df['is_SpliceSiteNoncanonRegion'] = np.where((merged_df['is_NoncodingRegion'] == 1) &
                                                  (merged_df['is_UTRsRegion'] == 0) &
                                                  (merged_df['is_PromoterRegion'] == 0) &
                                                  (merged_df['is_IntronRegion'] == 0),
                                                  merged_df['is_SpliceSiteNoncanonRegion'], 0)
        ## Intergenic
        patterns = ['downstream_gene_variant', 'intergenic_variant']
        merged_df['is_IntergenicRegion'] = merged_df['Consequence'].str.contains('|'.join(patterns)).astype(int)
        merged_df['is_IntergenicRegion'] = np.where((merged_df['is_NoncodingRegion'] == 1) &
                                                  (merged_df['is_UTRsRegion'] == 0) &
                                                  (merged_df['is_PromoterRegion'] == 0) &
                                                  (merged_df['is_IntronRegion'] == 0) &
                                                  (merged_df['is_SpliceSiteNoncanonRegion'] == 0),
                                                  merged_df['is_IntergenicRegion'], 0)
        ## lincRNA
        merged_df['is_lincRnaRegion'] = np.where((merged_df['is_NoncodingRegion'] == 1) &
                                                  (merged_df['is_UTRsRegion'] == 0) &
                                                  (merged_df['is_PromoterRegion'] == 0) &
                                                  (merged_df['is_IntronRegion'] == 0) &
                                                  (merged_df['is_SpliceSiteNoncanonRegion'] == 0) &
                                                  (merged_df['is_IntergenicRegion'] == 0) &
                                                  (merged_df['lincRNA'] == 1) &
                                                  (merged_df['ProteinCoding'] == 0),
                                                  1, 0)
        ## Others
        merged_df['is_OtherTranscriptRegion'] = np.where((merged_df['is_NoncodingRegion'] == 1) &
                                                  (merged_df['is_UTRsRegion'] == 0) &
                                                  (merged_df['is_PromoterRegion'] == 0) &
                                                  (merged_df['is_IntronRegion'] == 0) &
                                                  (merged_df['is_SpliceSiteNoncanonRegion'] == 0) &
                                                  (merged_df['is_IntergenicRegion'] == 0) &
                                                  (merged_df['is_lincRnaRegion'] == 0) &
                                                  (merged_df['ProteinCoding'] == 0),
                                                  1, 0)
        self._result = merged_df

    def allocate_variants(self, category: pd.DataFrame):
        if category['variant_type'] == 'All':
            filtered_result = self._result[self._result['CLASS'].isin(['SNV', 'Deletion', 'Insertion'])]
        elif category['variant_type'] == 'SNV':
            filtered_result = self._result.query('CLASS == "SNV"')
        else:
            filtered_result = self._result[self._result['CLASS'].isin(['Deletion', 'Insertion'])]

        if category['gene_set'] != 'Any':
            filtered_result = filtered_result[filtered_result[category['gene_set']] == 1]
        if category['functional_score'] != 'All':
            filtered_result = filtered_result[filtered_result[category['functional_score']] == 1]
        if category['gencode'] != 'Any':
            filtered_result = filtered_result[filtered_result['_'.join(['is', category['gencode']])] == 1]
        if category['functional_annotation'] != 'Any':
            filtered_result = filtered_result[filtered_result[category['functional_annotation']] == 1]
        filtered_result['CATEGORY'] = category['Category']
    
        return(filtered_result)

    def filter_variants(self):
        print_progress(f"Filter variants in {self.category_set.shape[0]} categories")
        self.category_set[['variant_type', 'gene_set', 'functional_score', 'gencode', 'functional_annotation']] = self.category_set['Category'].str.split('_', expand=True)
        # Filter variants by categories and concatenate them vertically
        self._result = pd.concat(self.category_set.apply(lambda x: self.allocate_variants(category = x), axis=1).tolist(), axis=0)
    
    def remove_annotation_info(self):
        print_progress("No annotation information attached")
        if self.category_set_path :
            self._result = self._result.loc[:, ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'SAMPLE', 'CATEGORY']]
        else :
            self._result = self._result.loc[:, ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'SAMPLE']]
    
    def save_result(self):
        print_progress(f"Save the result to the file {self.result_path}")
        self._result.to_csv(self.result_path, sep='\t', compression='gzip', index=False)
    
    def run(self):
        self.annotate_variants()
        if self.category_set_path :
            self.filter_variants()
        if self.annotation_info is None :
            self.remove_annotation_info()
        else:
            print_progress("Annotation information attached")
        self.save_result()
        print_progress("Done")
        
        