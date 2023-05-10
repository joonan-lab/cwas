============================
Gene list
============================


The gene list is a matrix with genes in each row identifying whether a gene falls into each gene set in columns. This matrix must contain Ensembl gene ID as a column to match genes with the gene annotated VCF from annotation output. The user-provided gene sets are added as columns with each row filled with binary code (1, if the gene belongs to the gene set. 0, if not.). The gene list file should also contain protein-coding genes, pseudogenes, and long noncoding RNAs as columns so that CWAS-Plus can use these gene sets while categorizing. The updated columns should be added with these columns.



+--------------------+-----------+---------------+---------+--------------+
| gene_id            | gene_name | ProteinCoding | lincRNA | ASDTADAFDR03 |
+====================+===========+===============+=========+==============+
| ENSG00000000003.15 | TSPAN6    |1              | 0       | 0            |
+--------------------+-----------+---------------+---------+--------------+
| ENSG00000000005.6  | TNMD      |1              | 0       | 0            |
+--------------------+-----------+---------------+---------+--------------+
| ENSG00000000419.14 | DPM1      |1              | 0       | 0            |
+--------------------+-----------+---------------+---------+--------------+

