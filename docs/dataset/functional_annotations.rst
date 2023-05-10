============================
Functional annotation
============================


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

