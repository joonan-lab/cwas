#!/usr/bin/env python
import pandas as pd
import os,sys,gzip
from functools import partial
cimport cython
cimport numpy as np
import numpy as np
import scipy.stats as ss

def binomTest(np.ndarray[long] vals1, np.ndarray[long] vals2):
	cdef int n1, n2
	out = []
	for n1, n2 in zip(vals1, vals2):
		out.append(ss.binom_test(x = n1, n = n1 + n2, p=0.5, alternative='two-sided'))
	return out

def binomTest_onesided(np.ndarray[long] vals1, np.ndarray[long] vals2):
	cdef int n1, n2
	out = []
	for n1, n2 in zip(vals1, vals2):
		if n1 > n2:
			out.append(ss.binom_test(x = n1, n = n1 + n2, p=0.5, alternative='greater'))
		else:
			out.append(ss.binom_test(x = n2, n = n1 + n2, p=0.5, alternative='greater'))
	return out

def rr(np.ndarray[double] vals1, np.ndarray[double] vals2):
	cdef double n1, n2
	out = [np.float64(n1)/n2 for n1, n2 in zip(vals1, vals2)]
	return out

cpdef perm(np.ndarray[long] vals1, np.ndarray[long] vals2, list swap_index): # val1 is from probands ; val2 is from siblings
	cdef int n, val_p, val_s, i, total_p, total_s, perm_number

	perm_pro = []
	perm_sib = []
	perm_number = len(swap_index)

	for n in range(0, perm_number):
		index_swap_fam = np.where(np.array(swap_index[n])==1)
		index_keep_fam = np.where(np.array(swap_index[n])==0)
		perm_pro.append( vals1[index_keep_fam].sum() + vals2[index_swap_fam].sum() )
		perm_sib.append( vals1[index_swap_fam].sum() + vals2[index_keep_fam].sum() )

	perm_results = [binomTest( np.array(perm_pro), np.array(perm_sib) ), rr( np.array(perm_pro).astype(float), np.array(perm_sib).astype(float) ),
					perm_pro, perm_sib ]

	return perm_results

cpdef doperm(df_sumvar, df_burden, swap_index):
	cat = df_sumvar.columns.tolist()[0]
	print(cat)
	cdef double rr_original, perm_p, perm_rr, pval, rr
	cdef int perm_over, p, s
	cdef list perm_results

	permfile_pvalue = '.'.join(['perm_p',cat,'txt.gz'])
	permfile_rr = '.'.join(['perm_rr',cat,'txt.gz'])
	permfile_count = '.'.join(['perm_count',cat,'txt.gz'])

	df_cat = pd.merge(df_sumvar.loc[df_sumvar['Role']=='p1'][['Fam',cat]].rename(columns = {cat : 'p1'}), 
					df_sumvar.loc[df_sumvar['Role']=='s1'][['Fam',cat]].rename(columns = {cat : 's1'}), how='inner', on='Fam')
	
	rr_original = df_burden.loc[df_burden['Annotation_combo']==cat,'Adjusted_relative_risk']
	perm_number = 10000
	perm_over = 0

	perm_results = perm(df_cat['p1'].round(decimals=0).astype('int64').values, 
						df_cat['s1'].round(decimals=0).astype('int64').values, 
						swap_index[0:perm_number])
	
	## Calculate perm pvalue
	for perm_rr in perm_results[1]:
		if (rr_original >= 1 and perm_rr >= rr_original) or (rr_original < 1 and rr_original >= perm_rr):
			perm_over += 1 

	perm_p = perm_over/float(perm_number)
	pvals = perm_results[0]
	rrs = perm_results[1]
	perm_pro = perm_results[2]
	perm_sib = perm_results[3]

	## To get accurate p-value for categories with lower perm p
	if perm_p < 0.1:
		perm_number = 100000

		perm_results2 = perm(df_cat['p1'].round(decimals=0).astype('int64').values, 
						df_cat['s1'].round(decimals=0).astype('int64').values, 
						swap_index[10000:perm_number])

		## Calculate perm pvalue
		for perm_rr in perm_results2[1]:
			if (rr_original >= 1 and perm_rr > rr_original) or (rr_original < 1 and rr_original > perm_rr):
				perm_over += 1
		perm_p = perm_over/float(perm_number)
		pvals = pvals + perm_results2[0]
		rrs = rrs + perm_results2[1]
		perm_pro = perm_pro + perm_results2[2]
		perm_sib = perm_sib + perm_results2[3]

	## Save pvalues to file
	o = gzip.open(permfile_pvalue, 'wt')
	o.write(str(perm_p) + '\n')
	# o.write('\n'.join([str(pval) for pval, rr in zip(pvals, rrs) if rr >= 1 else str(-1 * pval)]))
	o.write('\n'.join([str(pval) if rr >=1 else str(-pval) for pval, rr in zip(pvals, rrs) ]))
	o.write('\n')
	o.close()

	## Save rrs to file
	o = gzip.open(permfile_rr, 'wt')
	o.write('\n'.join([str(rr) for rr in rrs]))
	o.write('\n')
	o.close()

	## Save count (pro and sib) to file
	o = gzip.open(permfile_count, 'wt')
	for p,s in zip(perm_pro, perm_sib):
		o.write('\t'.join([str(p), str(s)]) + '\n')
	o.close()


cpdef create_index(long no_fams):
	cdef int n
	list_idx = [[n for n in np.random.binomial(1, 0.5, size=no_fams)] for i in range(0,100000)]
	return list_idx

cpdef check_pToZ(np.ndarray[double] p):
	cdef np.ndarray[double] p1
	cdef np.ndarray[double] z
	# Transform to 1 - p
	p1 = 1 - p
	print(p1[:10])
	# Convert to z score
	z = ss.norm.ppf(p1)
	print(z[:10])
	return(z)