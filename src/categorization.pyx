#!/usr/bin/env python
from functools import partial

import pandas as pd


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
