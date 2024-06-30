.. _faq:

**********************
FAQ
**********************


Easy installation of VEP
######################################################

For easy installation of VEP, we provide the following commands for installation through conda.

For LOFTEE plugin, the installation location for cpanm, perl, and VEP should be the same. For example, for conda, all tools should be located under the same conda environment for VEP to locate other tools.

To install VEP version 110, users can use the following command. If users want to install other versions, users can modify the number `110` to their desired specific version.

(The provided commands are referenced from `VEP document <https://asia.ensembl.org/info/docs/tools/vep/script/vep_download.html>`.)

.. code-block:: solidity
  
  conda install -c bioconda perl-bioperl
  conda install perl-App-cpanminus
  conda install -c bioconda ensembl-vep=110
  cpanm --force Bio::Perl

  wget https://github.com/ucscGenomeBrowser/kent/archive/v335_base.tar.gz
  tar xzf v335_base.tar.gz

  export KENT_SRC=$PWD/kent-335_base/src
  export MACHTYPE=$(uname -m)
  export CFLAGS="-fPIC"
  export MYSQLINC=`mysql_config --include | sed -e 's/^-I//g'`
  export MYSQLLIBS=`mysql_config --libs`

  cd $KENT_SRC/lib
  echo 'CFLAGS="-fPIC"' > ../inc/localEnvironment.mk

  make clean && make
  cd ../jkOwnLib
  make clean && make

  ln -s $KENT_SRC/lib/x86_64/* $KENT_SRC/lib/

  cpanm Bio::DB::BigFile
  cpanm DBD::SQLite

When the following message appears when using VEP,

.. code-block:: solidity
  
  Compress::Raw::Zlib version 2.201 required--this is only version 2.105

Try running the following command.

.. code-block:: solidity
  
  conda update -c conda-forge perl-compress-raw-zlib

Users can check whether VEP is installed through `vep --help`. If VEP is installed, the message will appear as following.

.. code-block:: solidity

  #----------------------------------#
  # ENSEMBL VARIANT EFFECT PREDICTOR #
  #----------------------------------#

  Versions:
    ensembl              : 110.584a8f3
    ensembl-funcgen      : 110.24e6da6
    ensembl-io           : 110.b1a0d57
    ensembl-variation    : 110.d34d25e
    ensembl-vep          : 110.1

  Help: dev@ensembl.org , helpdesk@ensembl.org
  Twitter: @ensembl

  http://www.ensembl.org/info/docs/tools/vep/script/index.html

  Usage:
  ./vep [--cache|--offline|--database] [arguments]

  Basic options
  =============

  --help                 Display this message and quit

  -i | --input_file      Input file
  -o | --output_file     Output file
  --force_overwrite      Force overwriting of output file
  --species [species]    Species to use [default: "human"]

  --everything           Shortcut switch to turn on commonly used options. See web
                        documentation for details [default: off]
  --fork [num_forks]     Use forking to improve script runtime

  For full option documentation see:
  http://www.ensembl.org/info/docs/tools/vep/script/vep_options.html


How to configure ANNOTATION_KEY_CONFIG yaml file
######################################################

The **ANNOTATION_KEY_CONFIG** yaml file contains *functional_score* and *functional_annotation*. Each file is in bed format and should be written with their short name that represents them. This name will be further used in other analyses.

When setting the short name, avoid using underscores ('_'). Underscores are used to distinguish different domains within a single category. For example, a category 'A_B_C_D_E' will be recognized as five domains, but if the name of the category is 'A_B_C_D_E_F', it will cause error while association testing, as the category is divided into six domains.

An example for **ANNOTATION_KEY_CONFIG** yaml file looks like below.

.. code-block:: solidity

  functional_score:
    bed1.bed.gz: annot1
    bed2.bed.gz: annot2
  functional_annotation:
    bed3.bed.gz: annot_3 # Do not use underscores like this. Users can use 'annot3' instead.
    bed4.bed.gz: annot4

The files should be located inside **ANNOTATION_DATA_DIR**. For preparation step, these files should be indexed using tabix.

As an example, users can sort and index their bed file like below.

.. code-block:: solidity
  
  cat bed1.bed | sort -k1,1V -k2,2n -k3,3n -t$'\t' | bgzip -c > sorted.bed1.bed.gz
  tabix -p bed sorted.bed1.bed.gz




How to add or remove a gene set
######################################################

Modify the txt file used for gene set (**GENE_MATRIX**).

The gene set contains gene ID, gene name and columns that represents each set of genes. To add or remove gene sets, users can add or remove columns.

+--------------------+-----------+---------------+---------+--------------+
| gene_id            | gene_name | ProteinCoding | lincRNA | ASDTADAFDR03 |
+====================+===========+===============+=========+==============+
| ENSG00000000003.15 | TSPAN6    |1              | 0       | 0            |
+--------------------+-----------+---------------+---------+--------------+
| ENSG00000000005.6  | TNMD      |1              | 0       | 0            |
+--------------------+-----------+---------------+---------+--------------+
| ENSG00000000419.14 | DPM1      |1              | 0       | 0            |
+--------------------+-----------+---------------+---------+--------------+

Then, configure again.

.. code-block:: solidity

  cwas configuration -f

If users already have the annotated VCF (`*annotated.vcf.gz`) and gene set is the only thing they changed, then they can start from categorization step. Gene sets are first used in the categorization step.


How to add or remove a functional annotation or score
######################################################

Modify the **ANNOTATION_KEY_CONFIG** yaml file. Add a new line or remove previous lines from the file.

Then, configure again.

.. code-block:: solidity

  cwas configuration -f


The annotation of CWAS-Plus contains two steps: (1) VEP annotation, (2) BED custom annotation. If the output of each step already exists, then the step is skipped.

If users already have the VEP annotated file (`*.vep.vcf.gz`), they can start from annotation step.

CWAS-Plus skips VEP annotation if the VEP annotated file already exists.

.. code-block:: solidity

  cwas annotation -v INPUT.vcf -o_dir OUTPUT_DIR -p 8

However, before annotation, please **remove the annotated VCF** (`*annotated.vcf.gz`). If the annotated VCF exists, CWAS-Plus will also skip BED custom annotation step.

