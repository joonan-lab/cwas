# Category-wide association study (CWAS)

![CWAS CI Workflow](https://github.com/joonan-lab/cwas/actions/workflows/ci.yml/badge.svg)

**CWAS (Category-Wide Association Study)** is a data analytic tool to perform stringent association tests to find non-coding loci associated with autism spectrum disorder (ASD). CWAS runs category-based burden tests using de novo variants from whole genome sequencing data and diverse annotation data sets.

CWAS was used in the following papers.

- [An analytical framework for whole genome sequence association studies and its implications for autism spectrum disorder](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5961723/) (Werling _et al._, 2018)
- [Genome-wide de novo risk score implicates promoter variation in autism spectrum disorder](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6432922/) (An _et al._, 2018)
- CWAS-Plus: Estimating genome-wide evaluation of noncoding variation from whole genome sequencing data. (Kim et al., in preperation)

Here is _the original CWAS repository: [sanderslab/cwas](https://github.com/sanderslab/cwas)_

## Quickstart

### Data requirements

Users must prepare following data for CWAS because it is very essential but cannot be generated automatically. Here are details.

#### 1. Input VCF data (De novo variant list)

```
#CHROM  POS ID  REF ALT QUAL    FILTER  INFO
chr1    3747728 .        T       C       .       .       SAMPLE=11000.p1
chr1    38338861        .       C       A       .       .       SAMPLE=11000.p1
chr1    117942118       .      T       G       .       .       SAMPLE=11000.p1
```

- The input VCF data must follow the [specification of VCF](https://samtools.github.io/hts-specs/VCFv4.2.pdf).
- The _INFO_ field must contain a sample ID of each variant with this format `SAMPLE={sample_id}`.

#### 2. List of samples

|  SAMPLE  | PHENOTYPE |
| :------: | :-------: |
| 11000.p1 |   case    |
| 11000.s1 |   ctrl    |
| 11002.p1 |   case    |
| 11002.s1 |   ctrl    |

- CWAS requires the file like above listing sample IDs with its family IDs and phenotypes (Case=_case_, Control=_ctrl_).
- Here are details of the required format.
  - Tab separated
  - 2 essential columns: _SAMPLE_ and _PHENOTYPE_
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

We recomment using _[conda virtual environment](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html)_ to build environment for CWAS. Installing CWAS-Plus within a conda environment will prevent its installation in the global environment. When creating a conda environment, also install Python to enable local installations using pip. We recommend installing R too. Run the following statements in your shell.

##### pip

```bash
conda create -n cwas python=3.10 r-base=4.2.2
conda activate cwas
pip install cwas
```

##### github

```bash
conda create -n cwas python=3.10 r-base=4.2.2
conda activate cwas
git clone https://github.com/joonan-lab/cwas.git
cd cwas
pip install .
```


In addition, you must install _[Variant Effect Predictor (VEP)](https://www.ensembl.org/vep)_.

### CWAS Execution

#### 1. Start

Run this command.

```bash
cwas start
```

As default, this command creates _CWAS workspace_ in your home directory. The path is `$HOME/.cwas`. `$HOME/.cwas/configuration.txt` has also generated.

Alternatively, _CWAS workspace_ can be specified with `-w` option. The `configuration.txt` will also be generated in the specified _CWAS workspace_.

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
VEP=/path/to/your/vep
VEP_CACHE_DIR=/path/to/your/vep/cache/dir
VEP_CONSERVATION_FILE=/path/to/your/vep/resource
VEP_LOFTEE=/path/to/your/vep/resource
VEP_HUMAN_ANCESTOR_FA=/path/to/your/vep/resource
VEP_GERP_BIGWIG=/path/to/your/vep/resource
VEP_MIS_DB=/path/to/your/missense/database
VEP_MIS_INFO_KEY=
VEP_MIS_THRES=
```

The `ANNOTATION_DATA` is a directory that contains all the BED files from **[joonan-lab/cwas-dataset](https://github.com/joonan-lab/cwas-dataset)**.
The `VEP_MIS_DB` is a database that is used to define damaging missense variants. The `VEP_MIS_INFO_KEY` is an user-defined name of the database used to annotate variants. The `VEP_MIS_THRES` is a threshold for missense variants (missense variants with value>=threshold are defined as damaging missense variants).


After writing the above file, run this command.

```bash
cwas configuration
```

Following files will be generated in your home directory as default. If you specify _CWAS workspace_, the files will be located in the same directory as the `configuration.txt`.

```bash
.cwas
├── annotation-data
├── annotation_keys.yaml
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

You can adjust the number of worker processes with `-p`.

After running this, merged BED file and its index will be generated in your _CWAS workspace_.

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
cwas annotation -v /path/to/your/vcf -p 4
```

You can adjust the number of worker processes with `-p`.

Here is the result file.

```bash
.cwas
...
├── {Your VCF filename}.annotated.vcf.gz
...
```

#### 5. Categorization

This step categorize your variants using the annotation datasets. Run this command.

```bash
cwas categorization -i /path/to/your/annotated/vcf -p 4
```

You can adjust the number of worker processes with `-p`.

After running this, you will get...

```bash
.cwas
...
├── {Your VCF filename}.categorization_result.zarr
...
```

Categorized results are generated in [zarr format](https://zarr.readthedocs.io/en/stable/index.html). Outputs are easily stored and loaded with zarr.

#### 6. Burden Test (Binomial Test)

This step is for calculation of relative risks and p-values for each category. As a default, these tests are based on variant-level analysis. The `--use_n_carrier` option can be used for sample-level analysis.

##### Binomial test

This step runs category-based burden test using the categorization result. The type of the test is binomial test. Run this command.

```bash
cwas binomial_test -i /path/to/your/categorization/result -s /path/to/your/samples [-a /path/to/your/adj_factors]
```

`[]` means that this is optional. If `-a` option is not specified, this step will bypass the adjustment step.

After running this, you will get...

```bash
.cwas
...
├── {Your VCF filename}.burden_test.txt
├── {Your VCF filename}.burden_test.volcano_plot.pdf
├── {Your VCF filename}.category_counts.txt
├── {Your VCF filename}.category_info.txt
...
```

##### Permutation test

This step runs category-based permutations using the categorization result. Run this command.

```bash
cwas permutation_test -i /path/to/your/categorization/result -s /path/to/your/samples [-a /path/to/your/adj_factors] [-b]
```

If `-b` option is specified, this step will generate binomial p-values for each permutation. This p-values will be used for burden shift and DAWN analysis.

After running this, you will get...

```bash
.cwas
...
├── {Your VCF filename}.permutation_test.txt
├── {Your VCF filename}.binom_pvals.txt.gz
...
```
