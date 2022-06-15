# Category-wide association study (CWAS)

![CWAS CI Workflow](https://github.com/mwjjeong/cwas/actions/workflows/ci.yml/badge.svg)

**CWAS (Category-Wide Association Study)** is a data analytic tool to perform stringent association tests to find non-coding loci associated with autism spectrum disorder (ASD). CWAS runs category-based burden tests using de novo variants from whole genome sequencing data and diverse annotation data sets.

CWAS was used in the following papers.

- [An analytical framework for whole genome sequence association studies and its implications for autism spectrum disorder](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5961723/) (Werling _et al._, 2018)
- [Genome-wide de novo risk score implicates promoter variation in autism spectrum disorder](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6432922/) (An _et al._, 2018)

Here is _the original CWAS repository: [sanderslab/cwas](https://github.com/sanderslab/cwas)_

## Quickstart

### Data requirements

Users must prepare following data for CWAS because it is very essential but cannot be generated automatically. Here are details.

#### 1. Input VCF data (De novo variant list)

```
#CHROM  POS ID  REF ALT QUAL    FILTER  INFO
chr1    3747728 .        T       C       .       .       SAMPLE=11000.p1;BATCH=P231
chr1    38338861        .       C       A       .       .       SAMPLE=11000.p1;BATCH=P231
chr1    117942118       .      T       G       .       .       SAMPLE=11000.p1;BATCH=P231
```

- The input VCF data must follow the [specification of VCF](https://samtools.github.io/hts-specs/VCFv4.2.pdf).
- The _INFO_ field must contain a sample ID of each variant with this format `SAMPLE={sample_id}`.

#### 2. List of samples

|  SAMPLE  | FAMILY | PHENOTYPE |
| :------: | :----: | :-------: |
| 11000.p1 | 11000  |   case    |
| 11000.s1 | 11000  |   ctrl    |
| 11002.p1 | 11002  |   case    |
| 11002.s1 | 11002  |   ctrl    |

- CWAS requires the file like above listing sample IDs with its family IDs and phenotypes (Case=_case_, Control=_ctrl_).
- Here are details of the required format.
  - Tab separated
  - 3 essential columns: _SAMPLE_, _FAMILY_, and _PHENOTYPE_
  - A value in the _PHENOTYPE_ must be _case_ or _ctrl_.
- The values in the _SAMPLE_ must be matched with the sample IDs of variants in the input VCF file.

#### 3. List of adjustment factors (Optional)

|  SAMPLE  | AdjustFactor |
| :------: | :----------: |
| 11000.p1 |    0.932     |
| 11000.s1 |    1.082     |
| 11002.p1 |    0.895     |
| 11002.s1 |    1.113     |

- The file like above is required if you want to adjust the number of variants for each sample in CWAS.
- Here are details of the required format.
  - Tab separated
  - 2 essential columns: _SAMPLE_ and _AdjustFactor_
  - A value in the _AdjustFactor_ must be a float.
- The values in the _SAMPLE_ must be matched with the sample IDs of variants in the input VCF file.

You can get the examples of the above data requirements from **[joonan-lab/cwas-input-example](https://github.com/joonan-lab/cwas-input-example)**

#### 4. CWAS annotation files

You can install those file from this repository: **[joonan-lab/cwas-dataset](https://github.com/joonan-lab/cwas-dataset)**

```bash
git clone https://github.com/joonan-lab/cwas-dataset.git
```

Due to the sizes of BigWig files for conservation scores, you must install them manually. [Please follow this instruction](https://github.com/joonan-lab/cwas-dataset/blob/main/bw_recipe.md).

### Installation

CWAS uses _[conda virtual environment](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html)_ to build environment for CWAS. Run the following statements in your shell.

```bash
# In your directory where CWAS is installed
git clone https://github.com/mwjjeong/cwas.git
cd cwas
conda env create -f environment.yml -n cwas
conda activate cwas
python setup.py install
```

In addition, you must install _[Variant Effect Predictor (VEP)](https://www.ensembl.org/vep)_.

### CWAS Execution

#### 1. Start

Run this command.

```bash
cwas start
```

This command creates _CWAS workspace_ in your home directory. The path is `$HOME/.cwas`. `$HOME/.cwas/configuration.txt` has also generated.

```bash
.cwas
└── configuration.txt
```

#### 2. Configuration

Write the following information in the `$HOME/.cwas/configuration.txt`.

```bash
ANNOTATION_DATA_DIR=/path/to/your/dir
GENE_MATRIX=/path/to/your/file
ANNOTATION_KEY_CONFIG=/path/to/your/file
BIGWIG_CUTOFF_CONFIG=/path/to/your/file
VEP=/path/to/your/vep
```

The `ANNOTATION_DATA` is a directory that contains all the BED files and BigWig files from **[joonan-lab/cwas-dataset](https://github.com/joonan-lab/cwas-dataset)**.

After writing the above file, run this command.

```bash
cwas configuration
```

Following files will be generated in your home directory.

```bash
.cwas
├── annotation-data
├── annotation_cutoff_bw.yaml
├── annotation_key_bed.yaml
├── annotation_key_bw.yaml
├── category_domain.yaml
├── configuration.txt
├── gene_matrix.txt
└── redundant_category.txt
.cwas_env
```

#### 3. Preparation

This step merges the BED files to annotate variants. Run the following command.

```bash
cwas preparation -p 4
```

`4` is the number of worker processes. You can adjust this.

After running this, Merged BED file and its index will be generated in your CWAS workspace.

```bash
.cwas
...
├── merged_annotation.bed.gz
├── merged_annotation.bed.gz.tbi
...
```

#### 4. Annotation

This step annotate your VCF file using _VEP_. Run this command.

```bash
cwas annotation -v /path/to/your/vcf
```

Here is the result file.

```bash
.cwas
...
├── {Your VCF filename}.annotated.vcf
...
```

#### 5. Categorization

This step categorize your variants using the annotation datasets. Run this command.

```bash
cwas categorization -p 4
```

`4` is the number of worker processes. You can adjust this.

After running this, you will get...

```bash
.cwas
...
├── {Your VCF filename}.categorization_result.txt
...
```

#### 6. Burden Test (Binomial Test)

This step runs category-based burden tests using the categorization result. A type of this burden test is binomial test. Run this command.

```bash
cwas binomial_test -s /path/to/your/samples [-a /path/to/your/adj_factors]
```

`[]` means that this is optional. If `-a` option does not specified, this step will bypass the adjustment step.

After running this, you will get...

```bash
.cwas
...
├── {Your VCF filename}.burden_test.txt
...
```
