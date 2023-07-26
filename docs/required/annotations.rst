.. _overview:

=================================
Annotations
=================================

0. Overview of annotation datasets
########################################


To start CWAS-Plus, users should prepare a subset of or total variants that are selected as input from qualifying variants in Variant Call Format (VCF) format (for example, noncoding de novo variants). These variants contain the position, reference allele, alternate allele, and the sample identifications. For CWAS-Plus to utilize the input, the users need to annotate the variants in advance to the variant categorization. The Ensembl Variant Effect Predictor (VEP) is recommended for general output. Besides the basic resources needed for VEP (such as, cache file with human genome database), user-provided annotation datasets are required.

The annotation datasets for CWAS-Plus are divided into five major groups, regarding their genetic features. The groups are (1) variant type, (2) gene biotype, (3) gene set, (4) functional annotation, and (5) functional score. The variant type indicates whether variants are SNVs or insertion-deletions (indels). The gene biotype is defined by the location of the variant relative to genes, such as coding regions or noncoding regions. Conversely, the remaining three groups are based on the user-provided datasets, which are selected based on their connection with the phenotype. For example, the users can add disease risk genes to the gene set group and genomic intervals of disease-specific regulatory elements to the functional annotation group. The users can also utilize datasets that annotate scores into genome sequences, such as conservation scores, to the functional score group. Through the customization of phenotype-specific annotation datasets, the users can test multiple hypothesis regarding the phenotype (details in the Methods section of the paper).

With VEP annotation, variants will be annotated with their consequences and genes or nearest genes, when the variant is at downstream or upstream. Now, we have the variant type and gene biotype based on the length of the alleles and the location of the variant. As these two groups of categories are pre-determined, the users mostly focus on the other three category groups, gene set, functional annotations, and functional scores. These datasets are required to be prepared in text, bed, and bigwig formats, respectively.

For bed custom annotation, CWAS-Plus gathers all functional annotations and scores and merge them into a single bed file. This file contains three columns of positional information and one column of binary code, which indicates which annotation this region falls into. CWAS-Plus encodes this information in binary format in order to save space and reduce execution time. For example, 8 in the fourth column of the merged file refers to 00100, which is a binary code in opposite direction to the normal binary scale. 00100 means that is interval only falls into the third functional annotation dataset. This order is fixed by CWAS-Plus in advance, with alphabetical order of the names of the datasets.


1. Gene sets
###################

The gene set is a matrix with genes in each row identifying whether a gene falls into each gene set in columns. This matrix must contain Ensembl gene ID as a column to match genes with the gene annotated VCF from annotation output. The user-provided gene sets are added as columns with each row filled with binary code (1, if the gene belongs to the gene set. 0, if not.). The gene set file should also contain protein-coding genes, pseudogenes, and long noncoding RNAs as columns so that CWAS-Plus can use these gene sets while categorizing. The updated columns should be added with these columns.

+--------------------+-----------+---------------+---------+--------------+
| gene_id            | gene_name | ProteinCoding | lincRNA | ASDTADAFDR03 |
+====================+===========+===============+=========+==============+
| ENSG00000000003.15 | TSPAN6    |1              | 0       | 0            |
+--------------------+-----------+---------------+---------+--------------+
| ENSG00000000005.6  | TNMD      |1              | 0       | 0            |
+--------------------+-----------+---------------+---------+--------------+
| ENSG00000000419.14 | DPM1      |1              | 0       | 0            |
+--------------------+-----------+---------------+---------+--------------+


2. Functional annotations
############################

The users should prepare functional annotations in bed format. In preparation, users should be cautious that bed format is different from VCF, as the position of VCF starts from 1 while position of bed starts from 0. Note that if users wish to use datasets in formats other than bed for functional annotations, the conversion process may require adjusting chromosomal positions.

The bed file contains intervals of genomic features, such as epigenomic marks. The first three columns contain positional information, and the fourth column contains features.

+------+--------+--------+
|chr1  | 10000  |  10600 |
+------+--------+--------+
|chr1  | 79200  |  80000 |
+------+--------+--------+
|chr1  | 610420 | 612020 |
+------+--------+--------+
|chr1  | 631820 | 632020 |
+------+--------+--------+


3. Functional scores
#########################

The users can utilize functional scores, such as conservation scores or pathogenicity scores. This file is in bed format. Functional scores contain genomic intervals, same as functional annotations. The dataset should be filtered by a reliable threeshold that is required for assuring that a specific region have this feature, such as conservation. For example, cutoff of 2 was applied to filter conserved regions in phyloP database.

+----+--------+-------+
|chr1|  12180 | 12181 |
+----+--------+-------+
|chr1|  12181 | 12182 |
+----+--------+-------+
|chr1|  12183 | 12184 |
+----+--------+-------+
|chr1|  12185 | 12186 |
+----+--------+-------+

