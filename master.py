#this file can be run to take user input and pass it to the querying file
import argparse
import csv
import numpy as np
import os
import pandas as pd
import sys
import sqlite3
import warnings
warnings.simplefilter(action='ignore', category = FutureWarning) #pandas is picky about .loc

#ARGUMENTS
#if the python script is run without any flags, just output genenames, cv_R2_avg, rsid, and weights (most common things we use)
parser = argparse.ArgumentParser()

#string inputs
#if none, query all available
parser.add_argument("--dbs", type = str, action = "store", dest = "dbs", required = False, help = ".db you want to query.") #"db/gtex_v7_Whole_Blood_imputed_europeans_tw_0.5_signif.db"
parser.add_argument("--genes", type = str, action = "store", dest = "genes", required = False, help = "File containing genes (Ensembl IDs) separated by line.") #"practice_python_queries/genenames.txt"
parser.add_argument("--genenames", type = str, action = "store", dest = "genenames", required = False, help = "File containing gene names separated by line.") #"practice_python_queries/genenames.txt"
parser.add_argument("--out_prefix", type = str, action = "store", dest = "out_prefix", required = False, default = "SQL_Simplifier_output", help = "Prefix of the output .csv file. Default = 'SQL_Simplifier_output'") #"practice_python_queries/genenames.txt"

#boolean inputs
#include db, gene, genename columns
parser.add_argument("--db_col", action = "store_true", dest = "db_col", default = False, help = "Output the column of .db file of origin.")
parser.add_argument("--gene_col", action = "store_true", dest = "gene_col", default = False, help = "Output the column of genes (Ensembl IDs).")
parser.add_argument("--genename_col", action = "store_true", dest = "genename_col", default = False, help = "Output the column of gene names.")

#EXTRA
parser.add_argument("--n.snps.in.model", action = "store_true", dest = "n_snps_in_model", default = False, help = "Output the number of SNPs within the cis window that have non-zero weights, as found by elastic-net.")
parser.add_argument("--test_R2_avg", action = "store_true", dest = "test_R2_avg", default = False, help = "Output the average coefficient of determination when predicting values of the hold out fold during nested cross validation.")
parser.add_argument("--cv_R2_avg", action = "store_true", dest = "cv_R2_avg", default = False, help = "Output the average coefficient of determination for each of the hold out folds when cross-validation was performed on the entire data set.")
parser.add_argument("--rho_avg", action = "store_true", dest = "rho_avg", default = False, help = "Output the average correlation between predicted and observed on the hold out folds when doing nested cross-validation.")
parser.add_argument("--rho_zscore", action = "store_true", dest = "rho_zscore", default = False, help = "Output the transformation of rho_avg into a z-score using Stouffer's Method.")
parser.add_argument("--pred.perf.R2", action = "store_true", dest = "pred_perf_R2", default = False, help = "Output the rho_avg squared.")
parser.add_argument("--pred.perf.pval", action = "store_true", dest = "pred_perf_pval", default = False, help = "Output the p-value for rho_zscore.")

#WEIGHTS (this table depends on the EXTRA table)
parser.add_argument("--rsid", action = "store_true", dest = "rsid", default = False, help = "Output the rsids in the models of queried genes.")
parser.add_argument("--varID", action = "store_true", dest = "varID", default = False, help = "Output the variant IDs in the models of queried genes. These are string labels of the format chromosome_position_allele1_allele2_build. All varIDs are from build 37 of the HUman Reference Genome.")
parser.add_argument("--ref_allele", action = "store_true", dest = "ref_allele", default = False, help = "Output the reference alleles of the SNPs in the models of the queried genes.")
parser.add_argument("--eff_allele", action = "store_true", dest = "eff_allele", default = False, help = "Output the effect alleles of the SNPs in the models of the queried genes.")
parser.add_argument("--weight", action = "store_true", dest = "weight", default = False, help = "Output the weights for the SNPs that are used to calculate predicted expression for the gene. In predicting the expression for the gene, the weight is multiplied by the count, or estimated count, of the effect allele in individual. This value is added to all other weighted SNPs in the model.")

#SAMPLE INFO
parser.add_argument("--n_samples", action = "store_true", dest = "n_samples", default = False, help = "Output the number of samples used the make the .db file.")
parser.add_argument("--population", action = "store_true", dest = "population", default = False, help = "Output the population studied.")
parser.add_argument("--tissue", action = "store_true", dest = "tissue", default = False, help = "Output the tissue or MESA population from which RNA was sequenced.")

#THRESHOLDS
parser.add_argument("--test_R2_avg_thres", type = float, dest = "test_R2_avg_thres", default = 0, help = "Restrict the test_R2_avg to values above this threshold. Default = 0.")
parser.add_argument("--cv_R2_avg_thres", type = float, dest = "cv_R2_avg_thres", default = 0, help = "Restrict the cv_R2_avg to values above this threshold. Default = 0.")
parser.add_argument("--rho_avg_thres", type = float, dest = "rho_avg_thres", default = 0, help = "Restrict the rho_avg to values above this threshold. Default = 0.")
parser.add_argument("--pred.perf.R2_thres", type = float, dest = "pred_perf_R2_thres", default = 0, help = "Restrict the test_R2_avg to values above this threshold. Default = 0.")
parser.add_argument("--pred.perf.pval_thres", type = float, dest = "pred_perf_pval_thres", default = 1, help = "Restrict the pred_perf_pval to values below this threshold. Default = 1.")

args = parser.parse_args() #then pass these arguments to further things

###INPUT SANITATION
if args.dbs is None:
    print("No .db destination detected. Please input a .db destination using the --db flag.")
    sys.exit(1)

###GENES, GENENAMES
if args.genes is None and args.genenames is None:
    print("No list of genes has been supplied with --genes or --genenames. All genes in the model(s) will be queried.")
    query_genes = []
elif args.genes is not None and args.genenames is not None:
    print("Please select an input for only genes or only genenames and not both.")
    sys.exit(1)
elif args.genes is not None:
    query_genes = list(np.loadtxt(args.genes, dtype = "str", ndmin = 1))
elif args.genenames is not None:
    query_genes = list(np.loadtxt(args.genenames, dtype = "str", ndmin = 1))

###COLS
col_flags = []
if args.db_col:
    col_flags.append("db")
if args.gene_col:
    col_flags.append("gene")
if args.genename_col:
    col_flags.append("genename")

###EXTRA
extra_flags = [] #store the flags the user passes
if args.n_snps_in_model:
    extra_flags.append("n.snps.in.model")
if args.test_R2_avg:
    extra_flags.append("test_R2_avg")
if args.cv_R2_avg:
    extra_flags.append("cv_R2_avg")
if args.rho_avg:
    extra_flags.append("rho_avg")
if args.rho_zscore:
    extra_flags.append("rho_zscore")
if args.pred_perf_R2:
    extra_flags.append("pred.perf.R2")
if args.pred_perf_pval:
    extra_flags.append("pred.perf.pval")

###WEIGHTS
weights_flags = []
if args.rsid:
    weights_flags.append("rsid")
if args.varID:
    weights_flags.append("varID")
if args.ref_allele:
    weights_flags.append("ref_allele")
if args.eff_allele:
    weights_flags.append("eff_allele")
if args.weight:
    weights_flags.append("weight")

###SAMPLE INFO
sample_info_flags = []
if args.n_samples:
    sample_info_flags.append("n_samples")
if args.population:
    sample_info_flags.append("population")
if args.tissue:
    sample_info_flags.append("tissue")

#make sure the user gets what they want
query_flags = col_flags + extra_flags + weights_flags + sample_info_flags
if len(query_flags) == 0:
    print("No query flags have been passed. Program will output cv_R2_avg, rsid, and weights of all genes in the input models.")
else: 
    print("Flags queried: " + ", ".join(query_flags))

#making input .dbs into a list
args_dbs = args.dbs #.endswith doesn't like arguments
if args_dbs.endswith(".db"): #its a single .db file
    dbs = [(args_dbs)]
    print("Model queried: " + args_dbs.split("/")[-1]) #don't print the full path
else: #its (I assume) a folder
    if args_dbs.endswith("/"):
        folder_name = args_dbs
    else:
        folder_name = args_dbs + "/"
    dbs = []
    if not os.path.exists(folder_name):
        print("The .db path is invalid. Please input a valid .db path.")
        sys.exit(1)
    for file in os.listdir(folder_name): #find files in a folder - https://stackoverflow.com/questions/3964681/find-all-files-in-a-directory-with-extension-txt-in-python
        if file.endswith(".db"):
            dbs.append(folder_name + file)
    if len(dbs) == 0:
        print("No .db models were found in the input destination. Program exiting.")
        sys.exit(1)
    print("Models queried: " + ", ".join([_.replace(folder_name, "") for _ in dbs])) #no need to print the folder name multiple times

#THRES FLAGS (for filtering later)
test_R2_avg_thres = args.test_R2_avg_thres
cv_R2_avg_thres = args.cv_R2_avg_thres
rho_avg_thres = args.rho_avg_thres
pred_perf_R2_thres = args.pred_perf_R2_thres
pred_perf_pval_thres = args.pred_perf_pval_thres

#QUERYING
data = [] #List of lists .db files info to output for further pandas filtering and parsing
rsid = None
varID = None
ref_allele = None
eff_allele= None
weight = None
n_snps_in_model = None
test_R2_avg = None
cv_R2_avg = None
rho_avg = None
rho_zscore = None
pred_perf_R2 = None
pred_perf_pval = None
n_samples = None
population = None
tissue = None

print("Beginning querying.")
for db in dbs:
    conn = sqlite3.connect(db) 
    c = conn.cursor()
    query_db_genes = query_genes.copy() #specific to this iteration of running dbs (otherwise runs into issues)

    if args.genenames is not None: #translate gene names to ensembl ids
        query_genenames = []
        for genename in query_db_genes:
            c.execute("select gene from extra where genename = '" + genename + "';")
            for row in c:
                query_genenames.append(row[0])
        query_db_genes = query_genenames
    if len(query_db_genes) == 0: #if no query genes input, query all genes
        c.execute("select gene from extra;")
        for row in c:
            query_db_genes.append(row[0])
            
    c.execute("select n_samples, population, tissue from sample_info;") #model-level data
    for row in c:
        n_samples = row[0]
        population = row[1]
        tissue = row[2]        
        
    for gene in query_db_genes:
        c.execute("select [n.snps.in.model], test_R2_avg, cv_R2_avg, rho_avg, rho_zscore, [pred.perf.R2], [pred.perf.pval], genename from extra where gene = '" + gene + "';") #gene-level data
        for row in c:
            n_snps_in_model = row[0]
            test_R2_avg = row[1]
            cv_R2_avg = row[2]
            rho_avg = row[3]
            rho_zscore = row[4]
            pred_perf_R2 = row[5]
            pred_perf_pval = row[6]
            genename = row[7]

        c.execute("select rsid, varID, ref_allele, eff_allele, weight from weights where gene = '" + gene + "';") #SNP-level data
        for row in c:
            rsid = row[0]
            varID = row[1]
            ref_allele = row[2]
            eff_allele = row[3]
            weight = row[4]
            data.append([db.split("/")[-1], gene, genename, n_samples, population, tissue, n_snps_in_model, test_R2_avg, cv_R2_avg, rho_avg, rho_zscore, pred_perf_R2, pred_perf_pval, rsid, varID, ref_allele, eff_allele, weight])
    conn.close()
print("Completed querying. Parsing SQL output.")

#FILTERING
data_frame = pd.DataFrame(data) #make list of lists into dataframe
data_frame.columns = ["db", "gene", "genename", "n_samples", "population", "tissue", "n.snps.in.model", "test_R2_avg", "cv_R2_avg", "rho_avg", "rho_zscore", "pred.perf.R2", "pred.perf.pval", "rsid", "varID", "ref_allele", "eff_allele", "weight"] #give column names so user knows what they're looking at

#subset through thresholds
if test_R2_avg_thres > 0:
    data_frame = data_frame.loc[data_frame['test_R2_avg'] > test_R2_avg_thres]
if cv_R2_avg_thres > 0:
    data_frame = data_frame.loc[data_frame['cv_R2_avg'] > cv_R2_avg_thres]
if rho_avg_thres > 0:
    data_frame = data_frame.loc[data_frame['rho_avg'] > rho_avg_thres]
if pred_perf_R2_thres > 0: 
    data_frame = data_frame.loc[data_frame['pred.perf.R2'] > pred_perf_R2_thres]
if pred_perf_pval_thres < 1: 
    data_frame = data_frame.loc[data_frame['pred.perf.pval'] < pred_perf_pval_thres]

#picks out user specified flags from data frame
if len(query_flags) > 0:
    data_frame = data_frame.loc[:, query_flags]
else:
    data_frame = data_frame.loc[:, ["db", "gene", "genename", "cv_R2_avg", "rsid", "weight"]]
data_frame = data_frame.drop_duplicates() #remove duplicate rows

#print to csv
data_frame.to_csv(args.out_prefix + ".csv", na_rep = "NA", index = False, quoting = csv.QUOTE_NONE) 
print("Completed running SQL Simplifier. Have a nice day :).")  

'''
NOTE TO READER FROM ANGELA:
We decided to query all genes in the .db and then subsequently parse because:
1. The Wheeler Lab machine has so much memory size isn't an issue
2. Carlee seemed more comfortable with hardcoding than working with shifting variables
3. It's easier to have two people work on stuff if the stuff remains the same instead of variable
4. I'm more comfortable in Pandas than SQL
'''
