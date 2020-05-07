"""
Module for the categorization of the variants into CWAS categories and counting the number of variants in each category.
CWAS category is a combination of annotation terms from each of annotation groups.

There are currently 5 groups of annotation terms.

--- The groups of the annotation terms for CWAS ---
    1. Variant types (var_type)
    2. Conservation (cons)
    3. Gene lists (gene_list)
    4. GENCODE annotation categories (effect)
    5. Functional annotation categories (region)

"""
cimport numpy as cnp
import numpy as np
import pandas as pd


cpdef dict categorize_variant(variant_df: pd.DataFrame, dict category_dict, dict gene_list_dict):
    """ Categorize the variants in the input data frame into the combinations of the annotation terms in the 5 groups
    and return the dictionary that contains the distribution of the variants for the combinations. These combinations
    are defined as CWAS categories.

    :param variant_df: The DataFrame object that contains the annotated variants from VEP 
    :param category_dict: The dictionary from parsing the category configuration file
    :param gene_list_dict: The dictionary which key and value are a gene symbol and the set of the names of 
                               the gene lists where this gene is involved, respectively
    :return: The dictionary which key and value are a CWAS category and the number of the variants in this category, 
             respectively
    """
    cdef dict annot_terms_dict, annot_term_idx_dict, cat_result_dict
    cdef cnp.ndarray[long] var_type_annot_ints, cons_annot_ints, gene_list_annot_ints
    cdef cnp.ndarray[long] effect_annot_ints, region_annot_ints
    cdef long var_type_int, cons_int, gene_list_int, effect_int, region_int
    cdef list var_type_annots, cons_annots, gene_list_annots, effect_annots, region_annots
    cdef str var_type_annot, cons_annot, gene_list_annot, effect_annot, region_annot

    annot_terms_dict = {}
    annot_term_idx_dict = {}

    for annot_group in category_dict:
        annot_terms = list(category_dict[annot_group].keys())
        annot_terms_dict[annot_group] = annot_terms
        annot_term_idx_dict[annot_group] = get_idx_dict(annot_terms_dict[annot_group])

    # For annotating each variant with multiple annotation terms from each group efficiently,
    # "Annotation integer" (annot_int) is used.
    # Annotation integer: A bitwise representation of the annotation of each variant where each bit means
    # each annotation term
    # e.g. If a list of annotation terms is ['A', 'B', 'C'] and the annotation integer is 0b101, it means that
    # the variant is annotated as 'A' and 'B'.
    var_type_annot_ints = annot_var_type(variant_df, annot_term_idx_dict['var_type'])
    cons_annot_ints = annot_cons(variant_df, annot_term_idx_dict['cons'])
    gene_list_annot_ints = annot_gene_list(variant_df, annot_term_idx_dict['gene_list'], gene_list_dict)
    effect_annot_ints = annot_effect(variant_df, annot_term_idx_dict['effect'], gene_list_dict)
    region_annot_ints = annot_region(variant_df, annot_term_idx_dict['region'])

    # Categorize by the annotation terms for each variant
    cat_result_dict = {}

    for var_type_annot_int, cons_annot_int, gene_list_annot_int, effect_annot_int, region_annot_int in \
            zip(var_type_annot_ints, cons_annot_ints, gene_list_annot_ints, effect_annot_ints, region_annot_ints):
        var_type_annots = parse_annot_int(var_type_annot_int, annot_terms_dict['var_type'])
        cons_annots = parse_annot_int(cons_annot_int, annot_terms_dict['cons'])
        gene_list_annots = parse_annot_int(gene_list_annot_int, annot_terms_dict['gene_list'])
        effect_annots = parse_annot_int(effect_annot_int, annot_terms_dict['effect'])
        region_annots = parse_annot_int(region_annot_int, annot_terms_dict['region'])

        # Make combinations using the annotation terms
        for var_type_annot in var_type_annots:
            for cons_annot in cons_annots:
                for gene_list_annot in gene_list_annots:
                    for effect_annot in effect_annots:
                        for region_annot in region_annots:
                            cwas_cat = f'{category_dict["var_type"][var_type_annot]}' \
                                       f'_{category_dict["gene_list"][gene_list_annot]}' \
                                       f'_{category_dict["cons"][cons_annot]}' \
                                       f'_{category_dict["effect"][effect_annot]}' \
                                       f'_{category_dict["region"][region_annot]}'

                            if cat_result_dict.get(cwas_cat) is None:
                                cat_result_dict[cwas_cat] = 1
                            else:
                                cat_result_dict[cwas_cat] += 1

    return cat_result_dict


cdef dict get_idx_dict(list list_):
    """ Return a dictionary which key and value are an item of the input list and its index, respectively. """
    cdef dict idx_dict
    cdef int idx
    cdef str item

    idx_dict = {}
    idx = 0

    for item in list_:
        idx_dict[item] = idx
        idx += 1

    return idx_dict


# Functions for annotation of the variants
# The functions below make an annotation integer for each variant.
# Step 1: Annotate by the types (e.g. SNV) of variants
cdef cnp.ndarray[long] annot_var_type(variant_df: pd.DataFrame, dict var_type_annot_idx_dict):
    cdef cnp.ndarray[str] refs, alts
    cdef cnp.ndarray[bint] is_snv_arr
    cdef cnp.ndarray[long] annot_ints

    refs = variant_df['REF'].values
    alts = variant_df['ALT'].values

    is_snv_arr = ((np.vectorize(len)(refs) == 1) & (np.vectorize(len)(alts) == 1)).astype(np.int32)
    annot_int_conv = \
        lambda is_snv: 2 ** var_type_annot_idx_dict['SNV'] if is_snv else 2 ** var_type_annot_idx_dict['Indel']
    annot_ints = np.vectorize(annot_int_conv)(is_snv_arr)
    annot_ints += 2 ** var_type_annot_idx_dict['All']

    return annot_ints


# Step 2: Annotate by conservation scores
cdef cnp.ndarray[long] annot_cons(variant_df: pd.DataFrame, dict cons_annot_idx_dict):
    cdef cnp.ndarray[double] phylop_scores, phast_scores
    cdef cnp.ndarray[bint] is_phylop_cons_arr, is_phast_cons_arr
    cdef cnp.ndarray[long] annot_ints

    phylop_conv_func = lambda x: -2.0 if x == '' else max(map(float, x.split('&')))
    phast_conv_func = lambda x: 0.0 if x == '' else max(map(float, x.split('&')))

    phylop_scores = np.vectorize(phylop_conv_func)(variant_df['phyloP46wayVt'].values)
    phast_scores = np.vectorize(phast_conv_func)(variant_df['phastCons46wayVt'].values)

    is_phylop_cons_arr = (phylop_scores >= 2.0).astype(np.int32)
    is_phast_cons_arr = (phast_scores >= 0.2).astype(np.int32)

    annot_ints = \
        np.vectorize(lambda is_phylop_cons: 2 ** cons_annot_idx_dict['phyloP46wayVt'] if is_phylop_cons else 0)\
        (is_phylop_cons_arr)
    annot_ints += \
        np.vectorize(lambda is_phast_cons: 2 ** cons_annot_idx_dict['phastCons46wayVt'] if is_phast_cons else 0)\
        (is_phast_cons_arr)
    annot_ints += 2 ** cons_annot_idx_dict['All']

    return annot_ints


# Step 3: Annotate by the names of the gene lists where the variant-associated genes are involved
cdef cnp.ndarray[long] annot_gene_list(variant_df: pd.DataFrame, dict gene_annot_idx_dict, dict gene_list_set_dict):
    cdef cnp.ndarray[str] gene_symbols, gene_nearests, effects
    cdef cnp.ndarray[long] annot_ints
    cdef list annot_int_list
    cdef dict annot_int_dict
    cdef str symbol, nearest, effect, gene
    cdef long annot_int
    cdef set gene_list_set

    gene_symbols = variant_df['SYMBOL'].values
    gene_nearests = variant_df['NEAREST'].values
    effects = variant_df['Consequence'].values  # GENCODE annotations

    annot_int_list = []
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

        annot_int_list.append(annot_int)

    annot_ints = np.asarray(annot_int_list)
    annot_ints += 2 ** gene_annot_idx_dict['Any']

    return annot_ints


# Step 4: Annotate by the effects (GENCODE annotations)
cdef cnp.ndarray[long] annot_effect(variant_df: pd.DataFrame, dict effect_annot_idx_dict, dict gene_list_set_dict):
    cdef cnp.ndarray[str] gene_symbols, gene_nearests, effects, polyphens
    cdef cnp.ndarray[long] annot_ints
    cdef list annot_int_list
    cdef str symbol, nearest, effect, polyphen, gene
    cdef long annot_int
    cdef bint is_in_coding

    gene_symbols = variant_df['SYMBOL'].values
    gene_nearests = variant_df['NEAREST'].values
    effects = variant_df['Consequence'].values
    polyphens = variant_df['PolyPhen'].values

    annot_int_list = []

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

        annot_int_list.append(annot_int)

    annot_ints = np.asarray(annot_int_list)
    annot_ints += 2 ** effect_annot_idx_dict['Any']

    return annot_ints


# Step 5: Annotate by the annotation terms of the regions where the variants are
cdef cnp.ndarray annot_region(variant_df: pd.DataFrame, dict region_annot_idx_dict):
    cdef cnp.ndarray[long] annot_ints
    cdef cnp.ndarray[str] region_vals
    cdef str region

    annot_ints = np.zeros(len(variant_df.index), dtype=int)

    for region in region_annot_idx_dict:
        if region == 'Any':
            continue

        region_vals = variant_df[region].values

        if region.startswith('Yale_H3K27ac'):
            region_val_conv_func = lambda x: 0 if x == '' else max([int(y.split('_')[0]) for y in x.split('&')])
            annot_int_conv_func = lambda x: 2 ** region_annot_idx_dict[region] if x > 1 else 0
            annot_ints += np.vectorize(annot_int_conv_func)(np.vectorize(region_val_conv_func)(region_vals))
        else:
            annot_int_conv_func = lambda x: 0 if x == '' else 2 ** region_annot_idx_dict[region]
            annot_ints += np.vectorize(annot_int_conv_func)(region_vals)

    annot_ints += 2 ** region_annot_idx_dict['Any']

    return annot_ints


cdef list parse_annot_int(long annot_int, list all_annot_terms):
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
