.. CWAS-Plus documentation master file, created by
   sphinx-quickstart on Wed May 10 00:04:11 2023.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

=====================================
Welcome to CWAS-Plus's documentation!
=====================================

CWAS-Plus (Category-Wide Association Study) is a data analysis tool designed to conduct rigorous association tests for discovering noncoding associations in complex genomic disorders. It runs category-based burden tests using variants from whole-genome sequencing data and various annotation datasets. CWAS-Plus provides a user-friendly interface for efficient hypothesis testing and has promising implications for uncovering the pathophysiology of various genomic disorders.

Here are the reference papers:

* `An analytical framework for whole genome sequence association studies and its implications for autism spectrum disorder <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5961723/>`_ (Werling et al., 2018)
* `Genome-wide de novo risk score implicates promoter variation in autism spectrum disorder <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6432922/>`_ (An et al., 2018)
* CWAS-Plus: Estimating genome-wide evaluation of noncoding variation from whole-genome sequencing data. (Kim et al., in preperation)

Here is the original CWAS repository: `sanderslab/cwas <https://github.com/sanderslab/cwas>`_


--------------------
CWAS-Plus workflow
--------------------

.. figure:: ../images/CWAS_workflow.jpg
   :alt: CWAS-Plus workflow
   :width: 90%
   :align: center

   **A.** Overall workflow of CWAS-Plus. Variants extracted from whole-genome sequencing data of samples (Sample WGS) and random variants are used as inputs. The color of the circles next to the box indicates which input each process uses (light blue, Sample WGS; yellow, Random variants). Three steps inside blue box refers to the final outputs obtained from CWAS-Plus. **B-I**. Graphic descriptions of each process in CWAS-Plus. Red (case) and blue (control) indicate the phenotype. H. The network inside the purple circle shows an example of a subset of a network in a magnified scale. The color indicates the direction of burden in each category (red, case burden; blue, control burden). **I**. The color of the clusters indicates the scale of z-score, which represent the significance and the relative risk of the cluster. Darker red refers to higher case burden.


.. toctree::
   :maxdepth: 1
   :caption: CWAS-Plus tutorial

   quickstart/tutorial.rst

.. toctree::
   :maxdepth: 1
   :caption: CWAS-Plus requirements

   required/installation.rst

.. toctree::
   :maxdepth: 1
   :caption: CWAS-Plus datasets

   dataset/overview_dataset.rst
   dataset/genelist.rst
   dataset/functional_annotations.rst
   dataset/functional_scores.rst

.. toctree::
   :maxdepth: 1
   :caption: CWAS-Plus configuration

   config/configuration.rst

.. toctree::
   :maxdepth: 1
   :caption: CWAS-Plus steps

   steps/annotation.rst
   steps/categorization.rst
   steps/burden.rst
   steps/generation_randomset.rst

.. toctree::
   :maxdepth: 1
   :caption: Other useful commands

   utils/extract_variants.rst

.. toctree::
   :maxdepth: 1
   :caption: Connect

   GitHub <https://github.com/joonan-lab/cwas>

