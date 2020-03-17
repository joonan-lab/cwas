#!/usr/bin/env python
from functools import partial

import pandas as pd


# Functions for the information of annotation terms for CWAS
cpdef get_cwas_annot() -> (dict, dict):
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
    cdef list var_type_annots
    cdef list cons_annots
    cdef list gene_list_annots
    cdef list effect_annots
    cdef list region_annots

    cdef dict cwas_annot_dict
    cdef dict cwas_annot_idx_dict

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

    return cwas_annot_dict, cwas_annot_idx_dict


cpdef get_annot_term_conv() -> dict:
    """ Return the dictionary to convert the annotation terms in 'cons', 'gene_list', and 'region' groups to the 
    previous annotation terms. This is for backward compatibility with the previous CWAS.
    """
    cdef dict cons_term_conv
    cdef dict gene_list_term_conv
    cdef dict region_term_conv
    cdef dict cwas_cat_term_conv

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


## Category 1: variant type
cpdef check_varType(info, header_index):
    cdef dict out

    out = {'All': 1, 'SNV': 1, 'Indel': 0} if len(info[header_index['varType']['Cols_idx']['REF']]) == 1 and len(
        info[header_index['varType']['Cols_idx']['ALT']]) == 1 else {'All': 1, 'SNV': 0, 'Indel': 1}

    del info, header_index
    return out

## Category 3: Conservation scores
cpdef check_cons(info, header_index):
    cdef dict out
    cdef str score
    cdef float score1

    out = {'All': 1}
    for n in header_index['Cons']['Names'].keys():
        out[header_index['Cons']['Names'][n]] = 0
        score = ''
        score1 = 0

        if n == 'phyloP46wayVt':
            ## PhyloP
            score = info[header_index['Cons']['Cols_idx'][n]]
            if score == '':
                score1 = -2.0
            else:
                score1 = max([float(a) for a in info[header_index['Cons']['Cols_idx'][n]].split('&')])
                if float(score1) >= 2.0:
                    out[header_index['Cons']['Names'][n]] = 1

        elif n == 'phastCons46wayVt':
            ## PhastCons
            score = info[header_index['Cons']['Cols_idx'][n]]
            if score == '':
                score1 = 0.0
            else:
                score1 = max([float(a) for a in info[header_index['Cons']['Cols_idx'][n]].split('&')])
                if float(score1) >= 0.20:
                    out[header_index['Cons']['Names'][n]] = 1

        else:
            print('Unknown Cons')

    del info, header_index
    return out

cpdef check_effect_genelist(info, header_index):
    ## https://www.ensembl.org/info/genome/variation/predicted_data.html
    ## http://www.gencodegenes.org/releases/27.html
    cdef dict out_genes
    cdef dict out_effects
    cdef str genelist, e, Gene, NEAREST, PolyPhen, SYMBOL

    Gene = info[header_index['Effect']['Cols_idx']['SYMBOL']] \
        if info[header_index['Effect']['Cols_idx']['SYMBOL']] != '' \
        else info[header_index['Effect']['Cols_idx']['NEAREST']]
    out_genes = {'Any': 1}
    for genelist in sorted(header_index['geneList']['Names'].keys()):
        out_genes[header_index['geneList']['Names'][genelist]] = 1 if Gene in header_index['geneList']['List'][
            genelist] else 0

    ## Fill cats with zero
    out_genes = {'Any': 1}
    for genelist in header_index['geneList']['Names'].keys():
        out_genes[header_index['geneList']['Names'][genelist]] = 0

    out_effects = {'Any': 1}
    for e in header_index['Effect']['Names'].keys():
        out_effects[e] = 0

    e = info[header_index['Effect']['Cols_idx']['Consequence']]
    SYMBOL = info[header_index['Effect']['Cols_idx']['SYMBOL']]
    NEAREST = info[header_index['Effect']['Cols_idx']['NEAREST']]
    PolyPhen = info[header_index['Effect']['Cols_idx']['PolyPhen']]
    Gene = ''

    ## Category 2: Gene list
    if 'downstream_gene_variant' not in e and 'intergenic_variant' not in e:
        Gene = SYMBOL
        for genelist in sorted(header_index['geneList']['Names'].keys()):
            if Gene in header_index['geneList']['List'][genelist]:
                out_genes[header_index['geneList']['Names'][genelist]] = 1

    else:
        ## Intergenic
        Gene = NEAREST
        for genelist in sorted(header_index['geneList']['Names'].keys()):
            if Gene in header_index['geneList']['List'][genelist]:
                out_genes[header_index['geneList']['Names'][genelist]] = 1

    ## Category 4: Effect
    if Gene in header_index['geneList']['List']['geneSet_Protein_Coding']:
        # Coding
        if 'stop_gained' in e or 'splice_donor' in e or 'splice_acceptor' in e:
            out_effects['LoFRegion'] = 1
            out_effects['CodingRegion'] = 1

        elif 'frameshift_variant' in e or 'transcript_amplification' in e or 'transcript_ablation' in e:
            out_effects['LoFRegion'] = 1
            out_effects['FrameshiftRegion'] = 1
            out_effects['CodingRegion'] = 1

        elif 'missense_variant' in e or 'start_lost' in e or 'stop_lost' in e:
            out_effects['MissenseRegion'] = 1
            out_effects['CodingRegion'] = 1
            if 'probably_damaging' in PolyPhen:
                out_effects['MissenseHVARDRegionSimple'] = 1
            else:
                pass

        elif ('inframe_deletion' in e or 'inframe_insertion' in e):  # CHECK this consequence only in the protein coding transcript
            out_effects['InFrameRegion'] = 1
            out_effects['CodingRegion'] = 1

        elif 'synonymous_variant' in e:
            out_effects['SilentRegion'] = 1
            out_effects['CodingRegion'] = 1

        elif 'stop_retained_variant' in e or 'incomplete_terminal_codon_variant' in e or 'protein_altering_variant' in e or 'coding_sequence_variant' in e:
            out_effects['CodingRegion'] = 1

        ## Noncoding
        elif '_UTR_' in e:
            out_effects['NoncodingRegion'] = 1
            out_effects['UTRsRegion'] = 1

        elif 'upstream_gene_variant' in e:
            out_effects['NoncodingRegion'] = 1
            out_effects['PromoterRegion'] = 1

        elif 'intron_variant' in e:
            out_effects['NoncodingRegion'] = 1
            out_effects['IntronRegion'] = 1

        elif 'splice_region_variant' in e:
            out_effects['NoncodingRegion'] = 1
            out_effects['SpliceSiteNoncanonRegion'] = 1

        elif 'downstream_gene_variant' in e:
            out_effects['NoncodingRegion'] = 1
            out_effects['IntergenicRegion'] = 1

        elif 'intergenic_variant' in e:
            out_effects['NoncodingRegion'] = 1
            out_effects['IntergenicRegion'] = 1

        ## Noncoding transcripts
        elif 'non_coding_transcript_exon_variant' in e or 'mature_miRNA_variant' in e:
            out_effects['NoncodingRegion'] = 1
        else:
            print(e)

    elif Gene in header_index['geneList']['List']['geneSet_Antisense']:
        out_effects['NoncodingRegion'] = 1

        ## Noncoding but Antisense
        if '_UTR_' in e:
            out_effects['UTRsRegion'] = 1

        elif 'upstream_gene_variant' in e:
            out_effects['PromoterRegion'] = 1

        elif 'intron_variant' in e:
            out_effects['IntronRegion'] = 1

        elif 'splice_region_variant' in e:
            out_effects['SpliceSiteNoncanonRegion'] = 1

        elif 'downstream_gene_variant' in e:
            out_effects['IntergenicRegion'] = 1

        elif 'intergenic_variant' in e:
            out_effects['IntergenicRegion'] = 1

        else:
            out_effects['AntisenseRegion'] = 1

    elif Gene in header_index['geneList']['List']['geneSet_lincRNA']:
        out_effects['NoncodingRegion'] = 1

        ## Noncoding but lincRNA
        if '_UTR_' in e:
            out_effects['UTRsRegion'] = 1

        elif 'upstream_gene_variant' in e:
            out_effects['PromoterRegion'] = 1

        elif 'intron_variant' in e:
            out_effects['IntronRegion'] = 1

        elif 'splice_region_variant' in e:
            out_effects['SpliceSiteNoncanonRegion'] = 1

        elif 'downstream_gene_variant' in e or 'intergenic_variant' in e:
            out_effects['IntergenicRegion'] = 1

        else:
            out_effects['lincRnaRegion'] = 1

    else:
        ## Noncoding but other transcripts
        out_effects['NoncodingRegion'] = 1
        if '_UTR_' in e:
            out_effects['UTRsRegion'] = 1

        elif 'upstream_gene_variant' in e:
            out_effects['PromoterRegion'] = 1

        elif 'intron_variant' in e:
            out_effects['IntronRegion'] = 1

        elif 'splice_region_variant' in e:
            out_effects['SpliceSiteNoncanonRegion'] = 1

        elif 'downstream_gene_variant' in e or 'intergenic_variant' in e:
            out_effects['IntergenicRegion'] = 1

        else:
            out_effects['OtherTranscriptRegion'] = 1

    del info, header_index, genelist, e, Gene, NEAREST, PolyPhen, SYMBOL

    return [out_genes, out_effects]

## Category 5: Regional annotations
cpdef check_region(info, header_index):
    cdef dict out_reg
    cdef str r
    cdef int n_ind

    out_reg = {'Any': 1}
    for r in header_index['Reg']['Cols_idx'].keys():

        if 'Yale_H3K27ac' in r:
            if info[header_index['Reg']['Cols_idx'][r]] != '':
                n_ind = 0
                n_ind = max([int(a.split('_')[0]) for a in info[header_index['Reg']['Cols_idx'][r]].split('&')])
                out_reg[header_index['Reg']['Names'][r]] = 1 if n_ind > 1 else 0
            else:
                out_reg[header_index['Reg']['Names'][r]] = 0
        else:
            out_reg[header_index['Reg']['Names'][r]] = 1 if info[header_index['Reg']['Cols_idx'][r]] != '' else 0

    del info, header_index

    return out_reg

cpdef doCats(l, header_index):
    info = l.tolist()

    cdef dict varType
    cdef list gene_effect
    cdef dict genelist
    cdef dict cons
    cdef dict effect
    cdef dict reg

    varType, gene_effect, cons, reg = check_varType(info, header_index), check_effect_genelist(info,
                                                                                               header_index), check_cons(
        info, header_index), check_region(info, header_index)
    genelist, effect = gene_effect[0], gene_effect[1]

    cdef int a1, b1, c1, d1, e1
    catRes = pd.Series({'_'.join([a, b, c, d, e]): a1 * b1 * c1 * d1 * e1 \
                        for a, a1 in varType.items() \
                        for b, b1 in genelist.items() \
                        for c, c1 in cons.items() \
                        for d, d1 in effect.items() \
                        for e, e1 in reg.items()
                        })

    del l, header_index, varType, gene_effect, genelist, cons, effect, reg

    return catRes

def parCat(df, header_index):
    filename = '.'.join(['tmp_catego', df.SampleID.unique().tolist()[0], 'txt'])
    df.apply(partial(doCats, header_index=header_index), axis=1).sum(axis=0).to_csv(filename, sep=';', header=False)

    del df
    del header_index
