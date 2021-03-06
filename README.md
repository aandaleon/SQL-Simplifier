# SQL Simplifier 

## Overview
This repository contains a Python wrapper for SQLite3 to take in parameters from the user and automate the queries to PrediXcan database files, producing .csv files ready for further analysis.

## Software Requirements
* Linux
* [Python 3.6.7](https://www.python.org/downloads/) with the libraries:
  * [argparse](https://docs.python.org/3/library/argparse.html) 1.1
  * [csv](https://docs.python.org/3/library/csv.html) 1.0
  * [numpy](http://www.numpy.org/) 1.16.1
  * [os](https://docs.python.org/3/library/os.html) 3.6.7
  * [pandas](https://pandas.pydata.org/) 0.24.2
  * [sqlite3](https://docs.python.org/3/library/sqlite3.html) 2.6.0
  * [sys](https://docs.python.org/3/library/sys.html) 3.6.7

## Downloading the project
Open a terminal session and enter: `git clone https://github.com/aandaleon/SQL-Simplifier.git`

## Input files
### Required
* One database file or one folder containing database files
  * This program is calibrated for GTEx V7 and MESA .dbs
### Optional
* List of genes (Ensembl ids), one per row
* List of gene names, one per row
* Flags containing the information you want queried (see Program options)

When the program is run without any parameters or gene lists, it will query every gene in a model and output `db`, `gene`, `genename`, `cv_R2_avg`, `rsid`, and `weight`, the most common metrics utilized by the Wheeler lab. Default output will be in `SQL_Simplifier_output.csv`

## Example
* Query a list of genes in a folder of .db files without any query flags (will output into `SQL_Simplifier_output.csv`)
  * `python3 master.py --dbs example_data/ --genenames example_data/genenames.txt`

| db                                                     | gene               | genename | cv_R2_avg          | rsid        | weight               |
|--------------------------------------------------------|--------------------|----------|--------------------|-------------|----------------------|
| gtex_v7_Whole_Blood_imputed_europeans_tw_0.5_signif.db | ENSG00000130203.5  | APOE     | 0.0135600336620361 | rs2356537   | -0.151237337241466   |
| gtex_v7_Whole_Blood_imputed_europeans_tw_0.5_signif.db | ENSG00000130203.5  | APOE     | 0.0135600336620361 | rs11668687  | -0.00241744847729031 |
| gtex_v7_Whole_Blood_imputed_europeans_tw_0.5_signif.db | ENSG00000130203.5  | APOE     | 0.0135600336620361 | rs11673170  | -0.00232453795322572 |

* Query all genes in a single .db file and their `genename`, `cv_R2_avg`, `n.snps.in.model`, and `pred.perf.R2`, outputting into `gene_info.csv`
  * `python3 master.py --dbs example_data/AFA_imputed_10_peer_3_pcs_v2.db --genename_col --cv_R2_avg --n.snps.in.model --pred.perf.R2 --out_prefix gene_info`

| genename | n.snps.in.model | cv_R2_avg          | pred.perf.R2       |
|----------|-----------------|--------------------|--------------------|
| FUCA2    | 21              | 0.219505222129763  | 0.239011989288677  |
| ENPP4    | 82              | 0.396811312007829  | 0.411548308924569  |
| ANKIB1   | 26              | 0.0961397809054118 | 0.0890519595368423 |

* Query all genes in a single .db file with cv_R2_avg > 0.1 and their genenames, cv_R2_avg, rsids, weights, outputting into `cv_R2_avg_0.1.csv`
  * `python3 master.py --dbs example_data/gtex_v7_Whole_Blood_imputed_europeans_tw_0.5_signif.db --genename_col --cv_R2_avg --rsid --weight --cv_R2_avg_thres 0.1 --out_prefix cv_R2_avg_0.1`

| genename     | cv_R2_avg         | rsid       | weight                |
|--------------|-------------------|------------|-----------------------|
| ISG15        | 0.154111838799616 | rs1058161  | -0.11337283758081     |
| ISG15        | 0.154111838799616 | rs11804831 | -0.0126092887627783   |
| ISG15        | 0.154111838799616 | rs2477782  | -0.0644525079361206   |

## Program options
* **Input/output files**
  * `--db`: path to .db file or folder path you want to query
  * `--genes`: file containing gene (Ensembl IDs) separated by line
  * `--genenames`: file containing gene names separated by line
  * `--out_prefix`: output file prefix; will end in .csv

* **Inclusion parameters**
  * `--db_col`: Output the column of .db file of origin.
  * `--gene_col`: Output the column of genes (Ensembl IDs).
  * `--genename_col`: Output the column of gene names.
  * `--n.snps.in.model`: Output the number of SNPs within the cis window that have non-zero weights, as found by elastic net.
  * `--test_R2_avg`: Output the average coefficient of determination when predicting values of the hold out fold during nested cross validation.
  * `--cv_R2_avg`: Output the average coefficient of determination for each of the hold out folds when cross-validation was performed on the entire data set.
  * `--rho_avg`: Output the average correlation between predicted and observed on the hold out folds when doing nested cross-validation.
  * `--rho_zscore`: Output the transformation of rho_avg into a z-score using Stouffer's Method.
  * `--pred.perf.R2`: Output the rho_avg squared.
  * `--pred.perf.pval`: Output the p-value for rho_zscore.
  * `--rsid`: Output the rsids in the models of queried genes.
  * `--varID`: Output the variant IDs in the models of queried genes. These are string labels of the format chromosome_position_allele1_allele2_build. All varIDs are from build 37 of the HUman Reference Genome.
  * `--ref_allele`: Output the reference alleles of the SNPs in the models of the queried genes.
  * `--eff_allele`: Output the effect alleles of the SNPs in the models of the queried genes.
  * `--weight`: Output the effect alleles of the SNPs in the models of the queried genes.
  * `--n_samples`: Output the number of samples used the make the .db file.
  * `--population`: Output the population studied.
  * `--tissue`: Output the tissue or MESA population from which RNA was sequenced.

* **Filtering parameters**
  * `--test_R2_avg_thres` (default = 0): Restrict the test_R2_avg to values above this threshold.
  * `--cv_R2_avg` (default = 0): Restrict the cv_R2_avg to values above this threshold.
  * `--rho_avg_thres` (default = 0): Restrict the rho_avg to values above this threshold.
  * `--pred.perf.R2_thres` (default = 0): Restrict the test_R2_avg to values above this threshold.
  * `--pred.perf.pval_thres` (default = 1): Restrict the pred_perf_pval to values below this threshold.

## Project summaries
* [Project prompt](https://docs.google.com/presentation/d/1Xarn0oowpogUH9NmHpkTC-sKIEeIR__ac2_Azgp5Ilo/edit?usp=sharing)
* [Design document](https://github.com/aandaleon/SQLite3-wrapper/wiki/Design-Document)
* [Presentation 1](https://docs.google.com/presentation/d/1lDZIZd-aw6z8_7F-tAtBdKWFPR-5bLE_pI3pmGNPjFM/edit?usp=sharing)
* [Applications note](https://docs.google.com/document/d/1zZdlgaizWUCQ0v088a9LqwBZrsxtPLGiHGKr7nPlDOQ/edit?usp=sharing)
* [Final presentation](https://docs.google.com/presentation/d/19DFuks-hMrekXAK4OjyANLwLJQVzDdleuS0ILW1yA5k/edit#slide=id.g50c2a10715_0_3)

## Quick background and resources
Our project queries information from database files used by the program [PrediXcan](https://github.com/hakyim/PrediXcan). PrediXcan predicts gene expression by aggregate precalculated weights based on an individual's genotype that are stored in database files. These weights are calculated in various tissues and cohorts, such as the [Genotype-Tissue Expression Project (GTEx)](https://gtexportal.org/home/documentationPage) and the [Multi-Ethnic Study of Atherosclerosis](https://github.com/WheelerLab/DivPop), and all public database files are available at [predictdb.org](predictdb.org). A general layout of database files and descriptions for all information stored is available [here](https://s3.amazonaws.com/predictdb2/contributed/MESA-2018-05-v2/MESAdb_2018-05-28_updated_README.txt). We give the users the ability to query information from these models without prior knowledge of SQL and in a simple command line format. For more detail on the context, goals, and milestones of the project, please consult the [design document](https://github.com/aandaleon/SQLite3-wrapper/wiki/Design-Document).

## Authors
This program and documentation were created by BS and MS [Bioinformatics](https://www.luc.edu/bioinformatics/index.shtml) students Angela Andaleon, Carlee Bettler, and Sherya Wadhwa for Computational Biology (COMP 383/483) Spring 2019 with Dr. Catherine Putonti at Loyola University Chicago. The original project idea was proposed by Angela Andaleon, Peter Fiorica, Ryan Schubert, and Dr. Heather Wheeler for use by the [Wheeler Lab](https://hwheeler01.github.io/).

## References

* Gamazon ER‡, Wheeler HE‡, Shah KP‡, Mozaffari SV, Aquino-Michaels K, Carroll RJ, Eyler AE, Denny JC, GTEx Consortium, Nicolae DL, Cox NJ, Im HK. (2015) [A gene-based association method for mapping traits using reference transcriptome data.](https://www.nature.com/articles/ng.3367) Nature Genetics 47(9):1091-8. ‡Contributed equally.
* GTEx Consortium. (2013) [The Genotype-Tissue Expression (GTEx) project.](https://www.nature.com/articles/ng.2653) Nature Genetics 45, 580–585.
* Mogil LS, Andaleon A, Badalamenti A, Dickinson SP, Guo X, Rotter JI, Johnson WC, Im HK, Liu Y, Wheeler HE. (2018) [Genetic architecture of gene expression traits across diverse populations.](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1007586) PLOS Genetics 14(8):e1007586.
