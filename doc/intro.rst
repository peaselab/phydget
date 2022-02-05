.. _intro:

###############
Getting Started
###############

What is PhyDGET?
=========================
Phylogenetic Differential Gene Expression Tool (PhyDGET) is a method for analyzing the changes in transcriptome-wide expression levels gene by gene on a phylogeny. PhyDGET is a merger in method and thinking between Phylogenetic Comparative Methods and Differential Gene Expression.  PhyDGET uses the **R** Bioconductor packages *limma*, *voom*, and *edgeR* to normalize and filter heterogeneous RNA-seq data and prepare it for phylogenetic testing.  PhyDGET then parallelizes the passage of these data to BayesTrait-v3, which tests each gene's expression level as a quantitative trait evolving on the tree.  Using BayesTrait's likelihood ratio test framework, you can specify a range of branch rate-shifting models to test what genes' expression levels are changing on the targeted branches.  Note that this is NOT a traditional differential expression framework using linear regression models with control vs. treatment.  The goal of PhyDGET is to estimate species-specific expression levels for each gene and examine them phylogenetically as quantative traits. 

You can learn more in our paper: <URL>.

How do I cite this program?
===========================
Citation forthcoming. 

Please also include the URL <https://www.github.com/jbpease/phydget> in your methods section where the program is referenced.

**Please also be sure to cite R (CRAN), edgeR, limma, and BayesTrait**

Installation
============
No installation of the package itself is required, PhyDGET scripts should work as long as the Requirements below are satisfied.  The repository can be cloned or downloaded as a .zip file from GitHub.

Requirements
------------
* Python 3.x (2.7 will not work) https://www.python.org/downloads/
* Numpy for Python3 http://www.numpy.org
* Scipy for Python3 https://www.scipy.org 
* BayesTrait (V3.0.1) http://www.evolution.rdg.ac.uk/BayesTraitsV3.0.1/BayesTraitsV3.0.1.html
* R (latest version recommended) https://cran.r-project.org/
* Bioconductor for R (latest version recommended) https://www.bioconductor.org/
* *limma* via Bioconductor https://bioconductor.org/packages/release/bioc/html/limma.html
* *edgeR* via Bioconductor https://bioconductor.org/packages/release/bioc/html/edgeR.html

.. important:: Note that you can also prepare your data table with R prior to the main PhyDGET process for BayesTrait.  Thus it is possible to prepare your data on a local system with R installed, then migrate the voom-normalized data to a cluster with Python3 and BayesTrait for the main process.

Preparing your data
===================

Phylogeny
---------

An ultametric phylogenetic tree in Nexus format should be used.  We recommend preparation of trees using the *ape* package from R. 

RNA-Seq Count File
------------------

Data file should be a tab-separated file with a single header line and gene names in the first column.

Basic usage
===========

::
  python3 phydget.py --data DATAFLIE --tree TREEFILE --out OUTPUT --groups A+B

Required Paramters
------------------
``--tree``: File containing a single Nexus-formatted phylogeny.

``--data`` Tab-separated count file with column headers and first-column as gene names.

``--out`` Name of the output file

``--groups`` This specifies the branch or branches to be tested in the alternative model.

``--threads`` Number of threads for Python to use in parallel (does not pass through to BayesTrait, this is for Python multiprocessing of per-branch replicates in parallel)

Recommended Parameters
----------------------

``--bg``: All taxa not included in ``--groups`` should be listed here comma-separated.  

.. important:: Any taxa not in either ``--bg`` or ``--groups`` will be ignored for the analysis (pruned from tree and columns filtered from teh data).
