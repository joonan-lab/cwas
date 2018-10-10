#!/usr/bin/env python

__version__ = 0.2
__author__ = 'Joon An'
__date__ = 'October 5th, 2018'

description = '''
				Script for categorization and summerization.
			'''

import os,sys,argparse,glob,gzip
import pandas as pd
import multiprocessing as mp
from multiprocessing import Pool, cpu_count
from functools import partial
import pyximport; pyximport.install()
import pycatego_vep_cwas as ctest 


def main(infile, gene_matrix, number_threads, output_tag, AF_known, lof):
	## Print the run setting
	print '[Setting] Input file: %s' % infile
	print '[Setting] Gene matrix file: %s' % gene_matrix
	print '[Setting] Number of threads: %s' % number_threads
	print '[Setting] Output tag: %s' % output_tag
	print '[Progress] Loading the input file into the data frame' + '\n'

	## Load data
	colnames = ['CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO']
	if '.gz' in infile:
		df = pd.io.parsers.read_csv(infile, sep='\t', index_col=False, 
									comment='#', compression='gzip',
									names=colnames, header=None)
		## Get the columns of annotations
		with gzip.open(infile) as fh:
			for l in fh:
				## Get annotations from the commented line
				if '##INFO=<ID=CSQ' in l:
					colnames_annots = l.split(' Format: ')[1].split('"')[0].split('|')
					break
	else:
		df = pd.io.parsers.read_csv(infile, sep='\t', comment='#', index_col=False,
									names=colnames, header=None)
		## Get the columns of annotations
		with open(infile) as fh:
			for l in fh:
				if '##INFO=<ID=CSQ' in l:
					colnames_annots = l.split(' Format: ')[1].split('"')[0].split('|')
					break
	print 'DataFrame of input is %s' % (df.shape,) 


	## Remove columns that are not needed in categorization
	cols_rm = ["CHROM", "POS", "QUAL", "FILTER", "INFO",
				"Allele", "Allele_Rm", "IMPACT","Gene","Feature_type","Feature",
				"EXON","INTRON","HGVSc","HGVSp","cDNA_position","CDS_position","Protein_position","Amino_acids",
				"Codons","Existing_variation","STRAND","FLAGS","SYMBOL_SOURCE","HGNC_ID",
				"CANONICAL", "TSL", "APPRIS", "CCDS", "SOURCE",
				"CADD_RAW","CADD_PHRED", "gnomADg"]

	df[colnames_annots] = df['INFO'].str.split('|', expand=True)
	df[['SampleID', 'Batch', 'Allele_Rm']] = df['Allele'].str.split(';', expand=True)
	df = df[df.columns[~df.columns.isin(cols_rm)]]

	## (Optional) Filter known gnomad variants
	if AF_known == 'No':
		df = df[df['gnomADg_AF']=='']
		print '[Progress] Filtering known variants'
		print 'After filtering DataFrame of input is %s' % (df.shape,) 
	elif AF_known == 'Only':
		df = df[df['gnomADg_AF']!='']
		print '[Progress] Filtering unknown variants'
		print 'After filtering DataFrame of input is %s' % (df.shape,) 
	else:
		print '[Progress] Keep known variants'

	## Get sample information
	samples = sorted(df.SampleID.unique())
	print '[Progress] Total %s samples are ready for analysis' % str(len(samples))
	print samples[:10]
	print samples[-10:]

	# ## Save to local 
	# outfile_temp = '.'.join(['result','temp', output_tag, 'txt'])
	# df.to_csv(outfile_temp, sep='\t', index=False)

	## Creating the header information
	header_index = ctest.get_col_index(list(df.columns), gene_matrix)
	no_cats = ctest.buildCats(header_index)
	print '[Progress] Combined the annotations. Total %s domains' % len(no_cats)
	print '[Progress] Start processing' + '\n'

	## Split dataframes by samples
	s = df.groupby('SampleID')
	inputs = [s.get_group(x) for x in s.groups]
	print '[Progress] Split the dataframe by samples. Total %s dataframes' % len(inputs)

	## Creating a pool for parallel processing
	pool = mp.Pool(number_threads)
	pool.imap_unordered(partial(ctest.parCat, header_index=header_index), inputs) 
	pool.close() 
	pool.join() 
	print '[Progress] Calculation for each samples are done' + '\n'

	## Merging files after run
	fs = sorted(glob.glob('tmp_catego*'))
	print '[Progress] Start merging %s files' % str(len(fs))
	f = fs[0]
	fh = open(f).read().splitlines()
	out = []
	sampleID = f.replace('tmp_catego.','').replace('.txt','').replace('.','_')
	match = ['SampleID'] + [a.split(';')[0] for a in fh]
	out.append(match)
	match = [sampleID] + [ int(a.split(';')[1]) for a in fh ]
	out.append(match)

	for f in fs[1:]:
		match = [f.replace('tmp_catego.','').replace('.txt','').replace('.','_')]
		print fs.index(f)
		with open(f) as fh:
			for l in fh:
				match.append( int(l.rstrip('\n').split(';')[1]) )
		out.append(match)

	df = pd.DataFrame(out[1:],columns=out[0])

	## Writing out cat var matrix
	outfile_sumvar = '.'.join(['result','sumVar', output_tag, 'txt'])
	file_cats_zero = '.'.join(['result','cats_zero', output_tag, 'txt'])
	print '[Progress] Remove columns with zero variant' 
	df_zero = df.loc[:, (df == 0).any(axis=0)]
	o = open(file_cats_zero, 'w')
	for l in list(df_zero.columns):
		o.write(l + '\n')
	o.close()
	df = df.loc[:, (df != 0).any(axis=0)]
	print '[Progress] Writing the sumvar matrix to csv files' 
	df.to_csv(outfile_sumvar, sep='\t', index=False)

	## Clean up
	os.system('rm tmp*')

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description=description)
	parser.add_argument('-i','--infile', required=True, type=str, help='Input File')
	parser.add_argument('-g','--gene_matrix', required=False, type=str, help='Gene matrix File', default='geneMatrixV38_v1.txt')
	parser.add_argument('-t','--number_threads', required=False, type=int, help='Number of threads', default=1)
	parser.add_argument('-o','--output_tag', required=False, type=str, help='Output tag', default='output')
	parser.add_argument('-a','--AF_known', required=False, type=str, help='Keep known variants', default='Yes')
	parser.add_argument('-lof','--lof', required=False, type=str, help='Keep lof variants', default='No')
	args = parser.parse_args()
	main(args.infile, args.gene_matrix, args.number_threads, args.output_tag, args.AF_known, args.lof)

