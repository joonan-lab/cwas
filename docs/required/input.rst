.. _input:

*********************************
Input
*********************************


1. Input vcf file (variant list)
###################################

Prepare sorted variants in `vcf <https://samtools.github.io/hts-specs/VCFv4.2.pdf>`_ format. The INFO field in each entry includes a sample ID in ``SAMPLE={sample_id}`` format.

.. code-block:: solidity
    
    #CHROM  POS ID  REF ALT QUAL    FILTER  INFO
    chr1    3747728 .        T       C       .       .       SAMPLE=11000.p1
    chr1    38338861        .       C       A       .       .       SAMPLE=11000.p1
    chr1    117942118       .      T       G       .       .       SAMPLE=11000.p1


2. Sample information
###################################

Prepare the sample information in a text file, such as in txt or tsv format. The file must be tab separated. It also must contain three columns, *SAMPLE* and *PHENOTYPE*. A value in the *PHENOTYPE* muse be *case* or *ctrl*.
The values in the SAMPLE column are matched to the sample IDs of variants in the input vcf file.

+----------+-----------+
|  SAMPLE  | PHENOTYPE |
+==========+===========+
| 11000.p1 |   case    |
+----------+-----------+
| 11000.s1 |   ctrl    |
+----------+-----------+
| 11002.p1 |   case    |
+----------+-----------+
| 11002.s1 |   ctrl    |
+----------+-----------+
  
1. Adjustment factors
###################################

Adjustment factors are required if the users want to adjust the number of variants for each sample in CWAS-Plus. The file must be tab separated and must contain two columns, *SAMPLE* and *AdjustFactor*. A value in the *AdjustFactor* must be a float.
The values in the SAMPLE column must be matched to the sample IDs of variants in the input vcf file.

+----------+--------------+
| SAMPLE   | AdjustFactor |
+==========+==============+
| 11000.p1 | 0.932        |
+----------+--------------+
| 11000.s1 | 1.082        |
+----------+--------------+
| 11002.p1 | 0.895        |
+----------+--------------+
| 11002.s1 | 1.113        |
+----------+--------------+
  
  
  For example run, the above data are available at `joonan-lab/cwas-input-example <https://github.com/joonan-lab/cwas-input-example>`_.

.. code-block:: solidity
    
    cd $HOME
    git clone https://github.com/joonan-lab/cwas-input-example.git

