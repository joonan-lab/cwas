# Category-based burden test (An et al. 2018)

This README contains the run command for the category-based burden test in the An et al. (2018). The scripts and materials are distributed in Amazon Machine Images and you can deploy the workflow using Amazon Web Services (AWS). The AMI ID is `ami-0e1757919181cfb66`, name 'An 2018 to share', and has 180 Gib size. For the instance setting, we recommend large AWS instance series like `m4.x2large` or `m4.x4large`.


## Step 1. Annotation

Using Variant Effect Predictor (VEP; https://www.ensembl.org/vep), de novo variants will be annotated for genomic regions, functional regions, and functional/conservation scores.

Script: `run_vep.py`

Required arguments:

- -i, --infile = path to file listing variants (VCF format).
- -t, --number_threads = number of threads to use

```bash
# Example
python run_vep.py \
-i table.hq_dnvs.txt \
-t 2
```

## Step 2. Variant categorization

The annotated variants will be grouped by annotation categories and individuals. 

Script: `vep_pycatego.py`
Requirement: Cython (http://cython.org/). Please note that this script will use Cython. So you need to compile a Cython script (`pycatego_vep_cwas.pyx`) for the run.


Required arguments:

- -i, --infile = path to file listing annotated variants. Input file format described below.
- -g, --gene_matrix = path to file listing genesets. You can find the file from the AMI.

Optional arguments:

- -o, --output_tag = text string that will be used for naming output files, e.g. analysis name or date. Default is 'output'.
- -a, --AF_known = Keep the variants with known allele frequency by gnomAD. Default is 'Yes'.
- -lof, --lof = Keep LoF variants. Default is 'No'.
- -t, --number_threads = number of threads to use

Input file format:
Annotated variant list, output from run_vep.py, above.

```bash
python vep_pycatego.py \
-i table.hq_dnvs.vep_gene.txt \
-g geneMatrix_hg38.txt \
-t 16 \
-o cwas \
-a Yes
```

## Step 3. Trim cats

To speed up burden testing in the following step, categories with 0 variants will be removed from the file.

Script: `doperm.py`
Requirement: Cython (http://cython.org/). Please note that this script will use Cython. So you need to compile a Cython script (`doperm.pyx`) for the run.

```

python doperm.py \
-m trim \
-i result.sumVar.cwas.txt \
-r list_redundant_categories.txt

```

## Step 4. Run burden (binomial) test on annotation categories

Each remaining category after the trimming step 3 will be subject to burden test using the binomial test on total variant counts per category between cases and unaffected controls. 

Script: `getBinomtest.py`
Requirement: Cython (http://cython.org/). Please note that this script will use Cython. So you need to compile a Cython script (`doperm.pyx`) for the run.

Required arguments:
-i, --infile = path to file with filtered categories (or the non-filtered version of this file). Input file format described below.

Optional arguments:

- -o, --output_tag = text string that will be used for naming output files, e.g. analysis name or date. Default is 'out'.
- -a, --adj = file specifying adjustment to variant rate, e.g. from regression vs. covariates. Adjustment file format described below. Default is text string 'no', which will bypass this adjustment step.

Input file format:
File listing the number of variants per sample that belong to each annotation category, optionally with 0-variant and redundant or user-specified categories removed (trimming these additional categories will speed up run time for the remaining steps). Output file from run_filterCategories.py (or from run_categorizeVars.py), above. One row per sampleID, the first column should be `Sample ID` (format: `familyID_phenotype`) and one additional column for each annotation category.

Adjustment file format
File with 1 row per sampleID listing adjustment factor (e.g. from regression) to multiply by total variant counts. Requires columns named `SampleID` and "AdjustFactor"

```bash
python getBinomtest.py \
-i result.sumVar.cwas.trimmed.txt \
-a list_paternal_age.txt \
-r list_redundant_categories.txt \
-o cwas
```



## Step 5. Run permutation

Each remaining category after the trimming step 3 will be subject to permutated burden test using the binomial test on total variant counts per category between cases and unaffected controls.

Script: `doperm.py`
Requirement: Cython (http://cython.org/). Please note that this script will use Cython. So you need to compile a Cython script (`doperm.pyx`) for the run.

Required arguments:

- -i, --infile = path to file with filtered categories (or the non-filtered version of this file). Input file format described below.
- -m, --mode = Please choose a mode to do permutations (perm) or create family swap index (index). 

Optional arguments:

- -o, --output_tag = text string that will be used for naming output files, e.g. analysis name or date. Default is 'out'.
- -b, --burden_file = Non-permutaiton burden matrix file. Default is text string 'No'. This will only be required for burden testing.
- -a, --adj = file specifying adjustment to variant rate, e.g. from regression vs. covariates. Adjustment file format described below. Default is text string 'No', which will bypass this adjustment step.
- -r, --trim_file = File to remove redundant categories. Default is text string 'No', which will bypass this trimming step.
- -t, --number_threads = number of threads to use. Default it '4'
- -cats_start, --cats_start = Start position of categories. Default it '0'
- -cats_start, --cats_start = End position of categories. Default it '1'
- -s3_path, --s3_path = Copy path for s3. Default it 'No'


```
# Step 5-1. Creating family swap index
python doperm.py -m index -s swapFams.cwas.p -n 1902

# Step 5-2. Permutation
python doperm.py \
-m perm \
-i result.sumVar.cwas.trimmed.txt.gz \
-a list_paternal_age.txt \
-r list_redundant_categories.txt \
-b result.burden.cwas.noPerm.txt \
-s swapFams.cwas.p \
-cats_start 1 \
-cats_end 100 \
-t 16
```

## Step 6. Run burdenshift

This script will generate global burdenshift cross a large category (e.g. promoters, missenses). This is based on permutation files generated from the step 5. 

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

