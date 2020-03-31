# Category-based burden test (An et al. 2018)

*Note: This is a modified version of the original CWAS repository.*
*The original CWAS repository: [sanderslab/cwas](https://github.com/sanderslab/cwas)*

This README contains the run command for the category-based burden test in the An et al. (2018). The scripts and materials are distributed in Amazon Machine Images and you can deploy the workflow using Amazon Web Services (AWS). The AMI ID is `ami-0e1757919181cfb66`, name 'An 2018 to share', and has 180 Gib size. For the instance setting, we recommend large AWS instance series like `m4.x2large` or `m4.x4large`.

The AMI includes:

- List of de novo variants from An et al. (2018). 
- Script for annotation
- Annotation files
- List of redundant categories to be trimmed
- List of gene sets 
- List of large category for burdenshift analysis
- Adjustment factor for burden test 


## Step 1. Annotation

Using Variant Effect Predictor (VEP; https://www.ensembl.org/vep), de novo variants will be annotated for genomic regions, functional regions, and functional/conservation scores.

##### Script: 
`run_vep.py`

##### Required arguments:

- -i, --infile = Path to file listing variants (VCF format).
- -t, --number_threads = Number of threads to use

```bash
# Example
python run_vep.py \
-i table.hq_dnvs.txt \
-t 2
```

## Step 2. Variant categorization

The annotated variants will be grouped by annotation categories (or CWAS categories) and individuals (or samples). 

##### Script: 
`categorize.py`

##### Requirement: 
Cython (http://cython.org/). Please note that this script will use Cython. So you need to compile a Cython script (`categorization.pyx`) to run run.

##### Required arguments:

- -i, --infile = Path to file listing variants annotated by VEP. Input file format described below.
- -g, --gene_matrix = Path to file listing gene sets (gene matrix file). You can find the file from the AMI.

##### Optional arguments:

- -r, --rdd_cat_file = Path to file listing redundant CWAS categories. Default is *'' (empty string)*, which means there are no redundant categories so any of categories will not be removed.
- -o, --outfile = Path to the categorization result. Default path is *cwas_cat_result.txt*.
- -p, --num_proc = Number of processes for this script. Default is *1*.
- -a, --af_known = Keep the variants with known allele frequencies by gnomAD. Possible values are *{yes, no, only}*. Default is *yes*.

##### Argument formats
 - **Input file format**:
A list of variants annotated by VEP, output from `run_vep.py`, above.

```bash
# Help
./categorize -h

# Usage
./categorization.py \
-i IN_VCF_PATH \
-g GENE_MAT_PATH \
[-r RDD_CAT_PATH] \
[-o OUTFILE_PATH]  \
[-p NUM_PROC] \
[-a {yes, no, only}]

# Note: '[]' means they are optional arguments. '{}' contains possible values for the argument. 
```

## Step 3. Burden tests on the CWAS categories

Each category from the step 2 will be subject to burden tests using binomial tests or permutation tests on total variant counts per category between cases and unaffected controls. 

##### Script: 
`burden_test.py`

##### Required arguments:

- -i, --infile = Path to a result of the categorization. Input file format described below.

##### Optional arguments:

- -a, --adj_file = Path to a file specifying adjustment factors the the number of variants of each individual. Adjustment file format described below. Default is *'' (empty string)*, which will bypass this adjustment step.
- -o, --outfile = Path to a result of the burden tests. Default is *cwas_burden_binom_result.txt* if you have excuted binomial tests, or *cwas_burden_perm_result.txt* if you have excuted permutation tests.

##### Optional arguments *only for permutation tests*:

- -n, --num_perm = Number of label-swapping permutations. Default is *10,000*.
- -p, --num_proc = Number of processes for this script. Default is *1*.
- -po, --perm_outfile = Path to a file listing relative risks after permutations for each category. Default is *'' (empty string)*, which will not save the relative risks after permutations.

##### Argument formats

- **Input file format**: File listing the number of variants per sample that belong to each annotation category, optionally redundant or user-specified categories removed (trimming these additional categories will speed up run time for the remaining steps). Output file from `categorize.py`, above. Rows for each sample, one of columns should be *Sample ID* (format: *familyID_phenotype*) and the other columns for each annotation category.

- **Adjustment file format**: File with a row for each sample listing adjustment factors (e.g. from regression) to multiply by total variant counts. Requires columns named *SampleID* and *AdjustFactor*.

```bash
# For binomial tests
# Help
./burden_test.py binom -h

# Usage
./burden_test.py binom \
-i CAT_RESULT_PATH \
[-a ADJ_FILE_PATH] \
[-o OUTFILE_PATH]


# For permutation tests
# Help
./burden_test.py perm -h

# Usage
./burden_test.py perm \
-i CAT_RESULT_PATH \
[-a ADJ_FILE_PATH] \
[-o OUTFILE_PATH] \
[-n NUM_PERM] \
[-p NUM_PROC] \
[-po PERM_RR_PATH]

# Note: '[]' means they are optional arguments.
```

## Step 4. Run burdenshift

This script will generate global burdenshift cross a large category (e.g. promoters, missenses). This is based on permutation files generated from the step 3. 

Script: run_burdenShift.R

```
Rscript  \
run_burdenShift.R \
result.perm_p.cwas.txt.gz \ # Permutation p-value matrix from the Step 5.
result.perm_burdenshift.cwas.txt.gz \ # Burden matrix from the Step 4.
list_catSetMembership.txt \ # Matrix for a large category membership
10000 \ # Number of permutation used in the Step 5.
0.05 \ # p-value threshold
cwas  # Tag for output
```

## Step 5. Risk score analysis and Annotation clustering 

The scripts for risk score analysis and annotation clustering is not located in this repository. Please refer to https://github.com/lingxuez/WGS-Analysis

