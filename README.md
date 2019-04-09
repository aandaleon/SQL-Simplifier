# SQLite3-wrapper

## Overview
This repository contains a Python wrapper for SQLite3 to take in parameters from the user and automate the queries to PrediXcan database files, producing .csv files ready for further analysis.

## Software Requirements
* Linux
* [Python 3.x](https://www.python.org/downloads/) with the libraries: **find version numbers for all of these*
  * [argparse](https://docs.python.org/3/library/argparse.html)
  * [numpy](http://www.numpy.org/)
  * [os](https://docs.python.org/3/library/os.html)
  * [pandas](https://pandas.pydata.org/)
  * [sqlite3](https://docs.python.org/3/library/sqlite3.html)
  * [sys](https://docs.python.org/3/library/sys.html)

## Downloading the project
Open a terminal session and enter: `git clone https://github.com/aandaleon/SQLite3-wrapper.git; cd SQLite3-wrapper`

## Input files
### Required
* One database file or one folder containing database files
### Optional
* List of genes (Ensembl ids), one per row
* List of gene names, one per row

## Example

## Program options
**describe flags here**



## Project summaries
* [Project prompt](https://docs.google.com/presentation/d/1Xarn0oowpogUH9NmHpkTC-sKIEeIR__ac2_Azgp5Ilo/edit?usp=sharing)
* [Presentation 1](https://docs.google.com/presentation/d/1lDZIZd-aw6z8_7F-tAtBdKWFPR-5bLE_pI3pmGNPjFM/edit?usp=sharing)

## Quick background
Our project queries information from database files used by the program [PrediXcan](https://www.nature.com/articles/ng.3367). PrediXcan predicts gene expression by aggregate precalculated weights based on an individual's genotype that are stored in database files. These weights are calculated in various tissues and cohorts, such as the [Genotype-Tissue Expression Project (GTEx)](https://gtexportal.org/home/documentationPage), and all public database files are available at [predictdb.org](predictdb.org). A general layout of database files is available [here])(https://s3.amazonaws.com/predictdb2/contributed/MESA-2018-05-v2/MESAdb_2018-05-28_updated_README.txt). We give the users the ability to query information from these models without prior knowledge of SQL and in a simple command line format.

## Authors
This program and documentation were created by Angela Andaleon, Carlee Bettler, and Sherya Wadhwa for our Computational Biology class spring 2019 for Dr. Catherine Putonti. The original project idea was proposed by Angela Andaleon, Peter Fiorica, Ryan Schubert, and Dr. Heather Wheeler.
