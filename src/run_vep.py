#!/usr/bin/env python
__version__ = 0.3
__author__ = 'Joon An'
__date__ = 'October 5th, 2018'

description = '''
			Script for genomic and functional annotations using VEP.
			'''

import os,sys,glob,argparse
from os.path import expanduser

def main(infile, number_threads):

    ## Set the run
    if '/home/ec2-user' in expanduser("~"):
        vep_path = '/home/ec2-user/ensembl-vep/vep'
        custom_path = '/home/ec2-user/custom/'
    else:
        print(expanduser("~"))
        sys.exit(0)

    ## Split input file for a single chromosome
    chroms = ['chr' + str(n) for n in range(1,23)]

    for chrom in chroms:
        tmp = '.'.join(['tmp', chrom, 'vcf'])
        o = open(tmp, 'w')
        header = '\t'.join(['#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO'])
        o.write(header + '\n')
        with open(infile) as fh:
            for l in fh:
                l_chrom = l.split('\t')[0]
                if 'chr' not in l_chrom:
                    l_chrom = 'chr' + l_chrom
                else:
                    pass

                ## Check if the chromosome matching write to tmp file
                if l_chrom == chrom:
                    o.write(l)
                else:
                    pass
        o.close()

    ## Get a command for run
    cmds = []
    for chrom in chroms:
        tmp = '.'.join(['tmp', chrom, 'vcf'])
        outfile = 'output.vep.' + chrom + '.vcf'

        ## Add basic information
        cmd = [vep_path,
               '--assembly GRCh38 --offline',
               '--fork 2',
               '--force_overwrite',
               '--buffer_size 5000000000',
               '-i', tmp,
               '-o', outfile,
               '--vcf',
               '--no_stats',
               '--polyphen p',
               '--ccds',
               # '--hgvs', # it adds 50% run time as it checks the fasta file
               '--numbers',
               # '--domains',
               # '--regulatory',
               '--canonical',
               '--protein',
               '--biotype',
               '--uniprot',
               '--tsl',
               '--appris'
               # '--gene_phenotype',
               # '--af',
               # '--pubmed',
               # '--variant_class'
               # '--everything'
               ]

        ## Output only the most severe consequence per gene.
        cmd = cmd + ['--per_gene',
                     '--pick --pick_order canonical,appris,tsl,biotype,ccds,rank,length']

        ## Add options for nearest and distance
        cmd = cmd + ['--distance 2000',
                     '--nearest symbol',
                     '--symbol']

        ## Add custom annotations
        cmd = cmd + [
            ','.join(['-custom ' + custom_path + 'gnomad.genomes.r2.0.1.sites.GRCh38.noVEP.vcf.gz','gnomADg','vcf,exact,0,AF']),

            ','.join(['-custom ' + custom_path + 'phastCons46way.vertebrate.hg19ToHg38.bw','phastCons46wayVt','bigwig','overlap','0']),
            ','.join(['-custom ' + custom_path + 'phyloP46way.vertebrate.hg19ToHg38.bw','phyloP46wayVt','bigwig','overlap','0']),

            ','.join(['-custom ' + custom_path + 'ChmmState15.E1.Brain.hg19to38.sorted.bed.gz','ChmmState15_E1_Brain','bed','overlap','0']),
            ','.join(['-custom ' + custom_path + 'ChmmState15.E2.Brain.hg19to38.sorted.bed.gz','ChmmState15_E2_Brain','bed','overlap','0']),
            ','.join(['-custom ' + custom_path + 'ChmmState15.E3.Brain.hg19to38.sorted.bed.gz','ChmmState15_E3_Brain','bed','overlap','0']),
            ','.join(['-custom ' + custom_path + 'ChmmState15.E4.Brain.hg19to38.sorted.bed.gz','ChmmState15_E4_Brain','bed','overlap','0']),
            ','.join(['-custom ' + custom_path + 'ChmmState15.E5.Brain.hg19to38.sorted.bed.gz','ChmmState15_E5_Brain','bed','overlap','0']),
            ','.join(['-custom ' + custom_path + 'ChmmState15.E6.Brain.hg19to38.sorted.bed.gz','ChmmState15_E6_Brain','bed','overlap','0']),
            ','.join(['-custom ' + custom_path + 'ChmmState15.E7.Brain.hg19to38.sorted.bed.gz','ChmmState15_E7_Brain','bed','overlap','0']),
            ','.join(['-custom ' + custom_path + 'ChmmState15.E8.Brain.hg19to38.sorted.bed.gz','ChmmState15_E8_Brain','bed','overlap','0']),
            ','.join(['-custom ' + custom_path + 'ChmmState15.E9.Brain.hg19to38.sorted.bed.gz','ChmmState15_E9_Brain','bed','overlap','0']),
            ','.join(['-custom ' + custom_path + 'ChmmState15.E10.Brain.hg19to38.sorted.bed.gz','ChmmState15_E10_Brain','bed','overlap','0']),
            ','.join(['-custom ' + custom_path + 'ChmmState15.E11.Brain.hg19to38.sorted.bed.gz','ChmmState15_E11_Brain','bed','overlap','0']),
            ','.join(['-custom ' + custom_path + 'ChmmState15.E12.Brain.hg19to38.sorted.bed.gz','ChmmState15_E12_Brain','bed','overlap','0']),
            ','.join(['-custom ' + custom_path + 'ChmmState15.E13.Brain.hg19to38.sorted.bed.gz','ChmmState15_E13_Brain','bed','overlap','0']),
            ','.join(['-custom ' + custom_path + 'ChmmState15.E14.Brain.hg19to38.sorted.bed.gz','ChmmState15_E14_Brain','bed','overlap','0']),
            ','.join(['-custom ' + custom_path + 'ChmmState15.E15.Brain.hg19to38.sorted.bed.gz','ChmmState15_E15_Brain','bed','overlap','0']),

            ','.join(['-custom ' + custom_path + 'EpigenomeByGroup4.DNaseFDR001.Brain.hg19to38.sorted.bed.gz','EpigenomeByGroup4_DNaseFDR001_Brain','bed','overlap','0']),
            ','.join(['-custom ' + custom_path + 'EpigenomeByGroup4.H3K27ac.Brain.hg19to38.sorted.bed.gz','EpigenomeByGroup4_H3K27ac_Brain','bed','overlap','0']),
            ','.join(['-custom ' + custom_path + 'EpigenomeByGroup4.H3K27me3.Brain.hg19to38.sorted.bed.gz','EpigenomeByGroup4_H3K27me3_Brain','bed','overlap','0']),
            ','.join(['-custom ' + custom_path + 'EpigenomeByGroup4.H3K36me3.Brain.hg19to38.sorted.bed.gz','EpigenomeByGroup4_H3K36me3_Brain','bed','overlap','0']),
            ','.join(['-custom ' + custom_path + 'EpigenomeByGroup4.H3K4me1.Brain.hg19to38.sorted.bed.gz','EpigenomeByGroup4_H3K4me1_Brain','bed','overlap','0']),
            ','.join(['-custom ' + custom_path + 'EpigenomeByGroup4.H3K4me3.Brain.hg19to38.sorted.bed.gz','EpigenomeByGroup4_H3K4me3_Brain','bed','overlap','0']),
            ','.join(['-custom ' + custom_path + 'EpigenomeByGroup4.H3K9ac.Brain.hg19to38.sorted.bed.gz','EpigenomeByGroup4_H3K9ac_Brain','bed','overlap','0']),
            ','.join(['-custom ' + custom_path + 'EpigenomeByGroup4.H3K9me3.Brain.hg19to38.sorted.bed.gz','EpigenomeByGroup4_H3K9me3_Brain','bed','overlap','0']),

            ','.join(['-custom ' + custom_path + 'H3K27ac.160407.multiInt.filtBy2.merge.3col.hg19to38.sorted.bed.gz','H3K27ac_160407_multiInt_filtBy2_merge_3col','bed','overlap','0']),
            ','.join(['-custom ' + custom_path + 'atac.norep.160407.multiInt.filtBy2.merge.3col.hg19to38.sorted.bed.gz','atac_norep_160407_multiInt_filtBy2_merge_3col','bed','overlap','0']),

            ','.join(['-custom ' + custom_path + 'bamo_EncodeDNaseClustersUCSC.hg19to38.sorted.bed.gz','EncodeDNaseClustersUCSC','bed','overlap','0']),
            ','.join(['-custom ' + custom_path + 'bamo_EncodeTfbsClusterV2UCSC.hg19to38.sorted.bed.gz','EncodeTfbsClusterV2UCSC','bed','overlap','0']),
            ','.join(['-custom ' + custom_path + 'bamo_vistaEnhancerUCSC.hg19to38.sorted.bed.gz','vistaEnhancerUCSC','bed','overlap','0']),

            ','.join(['-custom ' + custom_path + 'fantom5.enhancer.robust.hg19to38.sorted.bed.gz','fantom5_enhancer_robust','bed','overlap','0']),

            ','.join(['-custom ' + custom_path + 'hg19_HARs_Doan2016.hg19to38.sorted.bed.gz','HARs_Doan2016','bed','overlap','0'])
        ]
        cmd = ' '.join(cmd)
        cmds.append(cmd)


    ## Run VEP in parallel
    import multiprocessing as mp
    pool = mp.Pool(number_threads)
    pool.map(os.system, cmds)
    pool.close()
    pool.join()

    ## Collates outputs into one

    outfile = infile.replace('txt','vep_gene.txt')
    o = open(outfile, 'w')
    fs = sorted(glob.glob('output*vcf'))
    for f in fs:
        with open(f) as fh:
            if fs.index(f) == 0:
                for l in fh:
                    o.write(l)
            else:
                for l in fh:
                    if l[0] == '#':
                        pass
                    else:
                        o.write(l)
    o.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('-i','--infile', required=True, type=str, help='Input File')
    parser.add_argument('-t','--number_threads', required=False, type=int, help='Number of threads', default=1)
    args = parser.parse_args()
    main(args.infile, args.number_threads)


