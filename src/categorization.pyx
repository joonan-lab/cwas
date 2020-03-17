#!/usr/bin/env python
cimport numpy as cnp
import numpy as np
import pandas as pd


cpdef dict cwas_cat(variant_df: pd.DataFrame, gene_list_set_dict: dict):
    """ Categorize the variants in the input data frame into the combinations of the annotation terms in the 5 groups
    and return the dictionary that contains the distribution of the variants for the combinations. These combinations
    are defined as CWAS categories.

    :param variant_df: The DataFrame object that contains the annotated variants from VEP 
    :param gene_list_set_dict: The dictionary which key and value are a gene symbol and the set of the names of 
                               the gene lists where this gene is involved, respectively.
    :return: The dictionary which key and value are a CWAS category and the number of the variants in this category, 
             respectively
    """
    cdef dict cwas_annot_dict, cwas_annot_idx_dict, cwas_annot_term_conv, cat_result_dict
    cdef cnp.ndarray var_type_annot_ints, cons_annot_ints, gene_list_annot_ints, effect_annot_ints, region_annot_ints
    cdef int var_type_int, cons_int, gene_list_int, effect_int, region_int
    cdef list var_type_annots, cons_annots, gene_list_annots, effect_annots, region_annots
    cdef str var_type_annot, cons_annot, gene_list_annot, effect_annot, region_annot

    cwas_annot_dict, cwas_annot_idx_dict = get_cwas_annot()
    cwas_annot_term_conv = get_annot_term_conv()

    # For annotating each variant with multiple annotation terms from each group efficiently,
    # "Annotation integer" (annot_int) is used.
    # Annotation integer: A bitwise representation of the annotation of each variant where each bit means
    # each annotation term
    # e.g. If a list of annotation terms is ['A', 'B', 'C'] and the annotation integer is 0b101, it means that
    # the variant is annotated as 'A' and 'B'.
    var_type_annot_ints = annot_var_type(variant_df, cwas_annot_idx_dict['var_type'])
    cons_annot_ints = annot_cons(variant_df, cwas_annot_idx_dict['cons'])
    gene_list_annot_ints = annot_gene_list(variant_df, cwas_annot_idx_dict['gene_list'], gene_list_set_dict)
    effect_annot_ints = annot_effect(variant_df, cwas_annot_idx_dict['effect'], gene_list_set_dict)
    region_annot_ints = annot_region(variant_df, cwas_annot_idx_dict['region'])

    # Categorize by the annotation terms for each variant
    cat_result_dict = {}

    for var_type_annot_int, cons_annot_int, gene_list_annot_int, effect_annot_int, region_annot_int in \
            zip(var_type_annot_ints, cons_annot_ints, gene_list_annot_ints, effect_annot_ints, region_annot_ints):
        var_type_annots = ['All'] + parse_annot_int(var_type_annot_int, cwas_annot_dict['var_type'])
        cons_annots = ['All'] + parse_annot_int(cons_annot_int, cwas_annot_dict['cons'])
        gene_list_annots = ['Any'] + parse_annot_int(gene_list_annot_int, cwas_annot_dict['gene_list'])
        effect_annots = ['Any'] + parse_annot_int(effect_annot_int, cwas_annot_dict['effect'])
        region_annots = ['Any'] + parse_annot_int(region_annot_int, cwas_annot_dict['region'])

        # Make combinations using the annotation terms
        for var_type_annot in var_type_annots:
            for cons_annot in cons_annots:
                for gene_list_annot in gene_list_annots:
                    for effect_annot in effect_annots:
                        for region_annot in region_annots:
                            cwas_cat = f'{var_type_annot}' \
                                       f'_{cwas_annot_term_conv["gene_list"][gene_list_annot]}' \
                                       f'_{cwas_annot_term_conv["cons"][cons_annot]}' \
                                       f'_{effect_annot}' \
                                       f'_{cwas_annot_term_conv["region"][region_annot]}'

                            if cat_result_dict.get(cwas_cat) is None:
                                cat_result_dict[cwas_cat] = 1
                            else:
                                cat_result_dict[cwas_cat] += 1

    return cat_result_dict


# Functions for the information of annotation terms for CWAS
cdef tuple get_cwas_annot():
    """ Return the dictionary for a list of annotation terms for CWAS and the dictionary 
    for indices of the annotation terms. 

    --- The groups of the annotation terms for CWAS ---
    1. Variant types (var_type)
    2. Conservation (cons)
    3. Gene lists (gene_list)
    4. GENCODE annotation categories (effect)
    5. Functional annotation categories (region)

    :returns
    1. A dictionary which key is a group name and value is the list of the annotation terms of this group.
    2. A dictionary which key is a group name and value is the dictionary
       which key is an annotation term of this group and value is the index of the annotation term.
    """
    cdef list var_type_annots, cons_annots, gene_list_annots, effect_annots, region_annots
    cdef dict cwas_annot_dict, cwas_annot_idx_dict
    cdef tuple returns

    # The lists of CWAS categories
    var_type_annots = [
        'SNV',
        'Indel'
    ]
    cons_annots = [
        'phyloP46wayVt',
        'phastCons46wayVt',
    ]
    gene_list_annots = [
        'ASD_TADA_FDR03',
        'Willsey_Union',
        'geneSet_PLI90Score',
        'geneSet_PSD',
        'geneSet_DDD',
        'geneSet_BE',
        'geneSet_CHD8_Common',
        'geneSet_FMRP_Darnell',
        'geneSet_Protein_Coding',
        'geneSet_Pseudogene',
        'geneSet_lincRNA',
        'geneSet_Antisense',
        'geneSet_Processed_Transcript',
    ]
    effect_annots = [
        'CodingRegion',
        'FrameshiftRegion',
        'InFrameRegion',
        'SilentRegion',
        'LoFRegion',
        'MissenseHVARDRegionSimple',
        'MissenseRegion',
        'NoncodingRegion',
        'SpliceSiteNoncanonRegion',
        'IntronRegion',
        'PromoterRegion',
        'IntergenicRegion',
        'UTRsRegion',
        'AntisenseRegion',
        'lincRnaRegion',
        'OtherTranscriptRegion',
    ]
    region_annots = [
        'ChmmState15_E1_Brain',
        'ChmmState15_E2_Brain',
        'ChmmState15_E3_Brain',
        'ChmmState15_E4_Brain',
        'ChmmState15_E5_Brain',
        'ChmmState15_E6_Brain',
        'ChmmState15_E7_Brain',
        'ChmmState15_E8_Brain',
        'ChmmState15_E9_Brain',
        'ChmmState15_E10_Brain',
        'ChmmState15_E11_Brain',
        'ChmmState15_E12_Brain',
        'ChmmState15_E13_Brain',
        'ChmmState15_E14_Brain',
        'ChmmState15_E15_Brain',
        'EpigenomeByGroup4_DNaseFDR001_Brain',
        'EpigenomeByGroup4_H3K27ac_Brain',
        'EpigenomeByGroup4_H3K27me3_Brain',
        'EpigenomeByGroup4_H3K36me3_Brain',
        'EpigenomeByGroup4_H3K4me1_Brain',
        'EpigenomeByGroup4_H3K4me3_Brain',
        'EpigenomeByGroup4_H3K9ac_Brain',
        'EpigenomeByGroup4_H3K9me3_Brain',
        'H3K27ac_160407_multiInt_filtBy2_merge_3col',
        'atac_norep_160407_multiInt_filtBy2_merge_3col',
        'HARs_Doan2016',
        'fantom5_enhancer_robust',
        'EncodeDNaseClustersUCSC',
        'EncodeTfbsClusterV2UCSC',
        'vistaEnhancerUCSC',
        'Yale_H3K27ac_CBC',
        'Yale_H3K27ac_DFC',
    ]

    cwas_annot_dict = {
        'var_type': var_type_annots,
        'cons': cons_annots,
        'gene_list': gene_list_annots,
        'effect': effect_annots,
        'region': region_annots,
    }
    cwas_annot_idx_dict = {
        'var_type': {annot: i for i, annot in enumerate(var_type_annots)},
        'cons': {annot: i for i, annot in enumerate(cons_annots)},
        'gene_list': {annot: i for i, annot in enumerate(gene_list_annots)},
        'effect': {annot: i for i, annot in enumerate(effect_annots)},
        'region': {annot: i for i, annot in enumerate(region_annots)},
    }
    returns = (cwas_annot_dict, cwas_annot_idx_dict)

    return returns


cdef dict get_annot_term_conv():
    """ Return the dictionary to convert the annotation terms in 'cons', 'gene_list', and 'region' groups to the 
    previous annotation terms. This is for backward compatibility with the previous CWAS.
    """
    cdef dict cons_term_conv, gene_list_term_conv, region_term_conv, cwas_annot_term_conv

    cons_term_conv = {
        'All': 'All',
        'phyloP46wayVt': 'phyloP46way',
        'phastCons46wayVt': 'phastCons46way',
    }

    gene_list_term_conv = {
        'Any': 'Any',
        'ASD_TADA_FDR03': 'ASDTADAFDR03',
        'Willsey_Union': 'WillseyUnion',
        'geneSet_PLI90Score': 'PLI90Score',
        'geneSet_PSD': 'PSD',
        'geneSet_DDD': 'DDD',
        'geneSet_BE': 'BE',
        'geneSet_CHD8_Common': 'CHD8Common',
        'geneSet_FMRP_Darnell': 'FMRPDarnell',
        'geneSet_Protein_Coding': 'ProteinCoding',
        'geneSet_Pseudogene': 'Pseudogene',
        'geneSet_lincRNA': 'lincRNA',
        'geneSet_Antisense': 'Antisense',
        'geneSet_Processed_Transcript': 'ProcessedTranscript',
     }

    region_term_conv = {
        'Any': 'Any',
        'ChmmState15_E1_Brain': 'ChmE1',
        'ChmmState15_E2_Brain': 'ChmE2',
        'ChmmState15_E3_Brain': 'ChmE3',
        'ChmmState15_E4_Brain': 'ChmE4',
        'ChmmState15_E5_Brain': 'ChmE5',
        'ChmmState15_E6_Brain': 'ChmE6',
        'ChmmState15_E7_Brain': 'ChmE7',
        'ChmmState15_E8_Brain': 'ChmE8',
        'ChmmState15_E9_Brain': 'ChmE9',
        'ChmmState15_E10_Brain': 'ChmE10',
        'ChmmState15_E11_Brain': 'ChmE11',
        'ChmmState15_E12_Brain': 'ChmE12',
        'ChmmState15_E13_Brain': 'ChmE13',
        'ChmmState15_E14_Brain': 'ChmE14',
        'ChmmState15_E15_Brain': 'ChmE15',
        'EpigenomeByGroup4_DNaseFDR001_Brain': 'EpiDNase',
        'EpigenomeByGroup4_H3K27ac_Brain': 'EpiH3K27ac',
        'EpigenomeByGroup4_H3K27me3_Brain': 'EpiH3K27me3',
        'EpigenomeByGroup4_H3K36me3_Brain': 'EpiH3K36me3',
        'EpigenomeByGroup4_H3K4me1_Brain': 'EpiH3K4me1',
        'EpigenomeByGroup4_H3K4me3_Brain': 'EpiH3K4me3',
        'EpigenomeByGroup4_H3K9ac_Brain': 'EpiH3K9ac',
        'EpigenomeByGroup4_H3K9me3_Brain': 'EpiH3K9me3',
        'H3K27ac_160407_multiInt_filtBy2_merge_3col': 'MidFetalH3K27ac',
        'atac_norep_160407_multiInt_filtBy2_merge_3col': 'MidFetalATAC',
        'HARs_Doan2016': 'HARs',
        'fantom5_enhancer_robust': 'EnhancerFantom',
        'EncodeDNaseClustersUCSC': 'EncodeDNase',
        'EncodeTfbsClusterV2UCSC': 'EncodeTFBS',
        'vistaEnhancerUCSC': 'EnhancerVista',
        'Yale_H3K27ac_CBC': 'YaleH3K27acCBC',
        'Yale_H3K27ac_DFC': 'YaleH3K27acDFC',
    }

    cwas_annot_term_conv = {
        'cons': cons_term_conv,
        'gene_list': gene_list_term_conv,
        'region': region_term_conv,
    }

    return cwas_annot_term_conv


# Functions for annotation of the variants
# The functions below make an annotation integer for each variant.
# Step 1: Annotate by the types (e.g. SNV) of variants
cdef cnp.ndarray annot_var_type(variant_df: pd.DataFrame, var_type_annot_idx_dict: dict):
    cdef cnp.ndarray refs, alts, is_snv_arr, annot_ints

    refs = variant_df['REF'].values
    alts = variant_df['ALT'].values

    is_snv_arr = (np.vectorize(len)(refs) == 1) & (np.vectorize(len)(alts) == 1)
    annot_int_conv = \
        lambda is_snv: 2 ** var_type_annot_idx_dict['SNV'] if is_snv else 2 ** var_type_annot_idx_dict['Indel']
    annot_ints = np.vectorize(annot_int_conv)(is_snv_arr)

    return annot_ints


# Step 2: Annotate by conservation scores
cdef cnp.ndarray annot_cons(variant_df: pd.DataFrame, cons_annot_idx_dict: dict):
    cdef cnp.ndarray phylop_scores, phast_scores, is_phylop_cons_arr, is_phast_cons_arr, annot_ints

    phylop_conv_func = lambda x: -2.0 if x == '' else max(map(float, x.split('&')))
    phast_conv_func = lambda x: 0.0 if x == '' else max(map(float, x.split('&')))

    phylop_scores = variant_df['phyloP46wayVt'].values
    phylop_scores = np.vectorize(phylop_conv_func)(phylop_scores)
    phast_scores = variant_df['phastCons46wayVt'].values
    phast_scores = np.vectorize(phast_conv_func)(phast_scores)

    is_phylop_cons_arr = phylop_scores >= 2.0
    is_phast_cons_arr = phast_scores >= 0.2

    annot_ints = \
        np.vectorize(lambda is_phylop_cons: 2 ** cons_annot_idx_dict['phyloP46wayVt'] if is_phylop_cons else 0)\
        (is_phylop_cons_arr)
    annot_ints += \
        np.vectorize(lambda is_phast_cons: 2 ** cons_annot_idx_dict['phastCons46wayVt'] if is_phast_cons else 0)\
        (is_phast_cons_arr)

    return annot_ints


# Step 3: Annotate by the names of the gene lists where the variant-associated genes are involved
cdef cnp.ndarray annot_gene_list(variant_df: pd.DataFrame, gene_annot_idx_dict: dict, gene_list_set_dict: dict):
    cdef cnp.ndarray gene_symbols, gene_nearests, effects
    cdef list annot_ints
    cdef dict annot_int_dict
    cdef str symbol, nearest, effect, gene
    cdef int annot_int
    cdef set gene_list_set

    gene_symbols = variant_df['SYMBOL'].values
    gene_nearests = variant_df['NEAREST'].values
    effects = variant_df['Consequence'].values  # GENCODE annotations

    annot_ints = []
    annot_int_dict = {}  # Key: a gene symbol, Value: its category integer (For memorization)

    for symbol, nearest, effect in zip(gene_symbols, gene_nearests, effects):
        gene = nearest if 'downstream_gene_variant' in effect or 'intergenic_variant' in effect else symbol
        annot_int = annot_int_dict.get(gene, 0)

        if annot_int == 0:
            gene_list_set = gene_list_set_dict.get(gene, set())

            if gene_list_set:
                for gene_cat in gene_annot_idx_dict:
                    if gene_cat in gene_list_set:
                        annot_int += 2 ** gene_annot_idx_dict[gene_cat]

        annot_ints.append(annot_int)

    return np.array(annot_ints)


# Step 4: Annotate by the effects (GENCODE annotations)
cdef cnp.ndarray annot_effect(variant_df: pd.DataFrame, effect_annot_idx_dict: dict, gene_list_set_dict: dict):
    cdef cnp.ndarray gene_symbols, gene_nearests, effects, polyphens
    cdef list annot_ints
    cdef str symbol, nearest, effect, polyphen, gene
    cdef int annot_int
    cdef bint is_in_coding

    gene_symbols = variant_df['SYMBOL'].values
    gene_nearests = variant_df['NEAREST'].values
    effects = variant_df['Consequence'].values
    polyphens = variant_df['PolyPhen'].values

    annot_ints = []

    for symbol, nearest, effect, polyphen in zip(gene_symbols, gene_nearests, effects, polyphens):
        gene = nearest if 'downstream_gene_variant' in effect or 'intergenic_variant' in effect else symbol
        gene_cat_set = gene_list_set_dict.get(gene, set())
        annot_int = 0
        is_in_coding = False

        if 'geneSet_Protein_Coding' in gene_cat_set:
            is_in_coding = True
            annot_int += 2 ** effect_annot_idx_dict['CodingRegion']

            # Coding region
            if 'stop_gained' in effect \
                    or 'splice_donor' in effect \
                    or 'splice_acceptor' in effect:
                annot_int += 2 ** effect_annot_idx_dict['LoFRegion']
            elif 'frameshift_variant' in effect \
                    or 'transcript_amplification' in effect \
                    or 'transcript_ablation' in effect:
                annot_int += 2 ** effect_annot_idx_dict['LoFRegion']
                annot_int += 2 ** effect_annot_idx_dict['FrameshiftRegion']
            elif 'missense_variant' in effect \
                    or 'start_lost' in effect \
                    or 'stop_lost' in effect:
                annot_int += 2 ** effect_annot_idx_dict['MissenseRegion']

                if 'probably_damaging' in polyphen:
                    annot_int += 2 ** effect_annot_idx_dict['MissenseHVARDRegionSimple']

            elif 'inframe_deletion' in effect \
                    or 'inframe_insertion' in effect:
                annot_int += 2 ** effect_annot_idx_dict['InFrameRegion']
            elif 'synonymous_variant' in effect:
                annot_int += 2 ** effect_annot_idx_dict['SilentRegion']
            elif 'stop_retained_variant' not in effect \
                    and 'incomplete_terminal_codon_variant' not in effect \
                    and 'protein_altering_variant' not in effect \
                    and 'coding_sequence_variant' not in effect:
                # Noncoding
                annot_int = 0
                is_in_coding = False  # This variant is in a noncoding region of the protein coding gene

        if not is_in_coding:
            annot_int += 2 ** effect_annot_idx_dict['NoncodingRegion']

            if '_UTR_' in effect:
                annot_int += 2 ** effect_annot_idx_dict['UTRsRegion']
            elif 'upstream_gene_variant' in effect:
                annot_int += 2 ** effect_annot_idx_dict['PromoterRegion']
            elif 'intron_variant' in effect:
                annot_int += 2 ** effect_annot_idx_dict['IntronRegion']
            elif 'splice_region_variant' in effect:
                annot_int += 2 ** effect_annot_idx_dict['SpliceSiteNoncanonRegion']
            elif 'downstream_gene_variant' in effect \
                    or 'intergenic_variant' in effect:
                annot_int += 2 ** effect_annot_idx_dict['IntergenicRegion']
            elif 'geneSet_Protein_Coding' not in gene_cat_set:
                if 'geneSet_Antisense' in gene_cat_set:
                    annot_int += 2 ** effect_annot_idx_dict['AntisenseRegion']
                elif 'geneSet_lincRNA' in gene_cat_set:
                    annot_int += 2 ** effect_annot_idx_dict['lincRnaRegion']
                else:
                    annot_int += 2 ** effect_annot_idx_dict['OtherTranscriptRegion']

        annot_ints.append(annot_int)

    return np.array(annot_ints)


# Step 5: Annotate by the annotation terms of the regions where the variants are
cdef cnp.ndarray annot_region(variant_df: pd.DataFrame, effect_annot_idx_dict: dict):
    cdef cnp.ndarray annot_ints, region_vals
    cdef str region

    annot_ints = np.zeros(len(variant_df.index), dtype=int)

    for region in effect_annot_idx_dict.keys():
        region_vals = variant_df[region].values

        if region.startswith('Yale_H3K27ac'):
            region_val_conv_func = lambda x: 0 if x == '' else max([int(y.split('_')[0]) for y in x.split('&')])
            region_vals = np.vectorize(region_val_conv_func)(region_vals)
            annot_int_conv_func = lambda x: 2 ** effect_annot_idx_dict[region] if x > 1 else 0
        else:
            annot_int_conv_func = lambda x: 0 if x == '' else 2 ** effect_annot_idx_dict[region]

        annot_ints += np.vectorize(annot_int_conv_func)(region_vals)

    return annot_ints


cdef list parse_annot_int(annot_int: int, all_annot_terms: list):
    """ From the list of annotation terms, choose the appropriate subset by parsing the annotation integer.
    """
    cdef int annot_idx
    cdef list annot_terms

    annot_idx = 0
    annot_terms = []

    while annot_int != 0:
        if annot_int % 2 == 1:
            annot_terms.append(all_annot_terms[annot_idx])

        annot_int >>= 1
        annot_idx += 1

    return annot_terms
