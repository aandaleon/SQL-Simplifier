#this file can be run to take user input and pass it to the querying file
import argparse
import numpy as np
import os
import pandas as pd
import sys
import sqlite3

########
#ANGELA#
########

#if the python script is run without any flags, just output genenames, cv_R2_avg, rsid, and weights (most common things we use)
parser = argparse.ArgumentParser()

#string inputs
#if none, query all available
parser.add_argument("--db", type = str, action = "store", dest = "db", required = False, help = ".db you want to query.") #"db/gtex_v7_Whole_Blood_imputed_europeans_tw_0.5_signif.db"
parser.add_argument("--genes", type = str, action = "store", dest = "genes", required = False, help = "File containing genes (Ensembl IDs) separated by line.") #"practice_python_queries/genenames.txt"
parser.add_argument("--genenames", type = str, action = "store", dest = "genenames", required = False, help = "File containing gene names separated by line.") #"practice_python_queries/genenames.txt"

#boolean inputs
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

#THRESHOLDS (may mess around with)
parser.add_argument("--test_R2_avg_thres", type = float, dest = "test_R2_avg_thres", default = 0, help = "Restrict the test_R2_avg to values above this threshold. Default = 0.")
parser.add_argument("--cv_R2_avg_thres", type = float, dest = "cv_R2_avg_thres", default = 0, help = "Restrict the cv_R2_avg to values above this threshold. Default = 0.")
parser.add_argument("--rho_avg_thres", type = float, dest = "rho_avg_thres", default = 0, help = "Restrict the rho_avg to values above this threshold. Default = 0.")
parser.add_argument("--pred.perf.R2_thres", type = float, dest = "pred_perf_R2_thres", default = 0, help = "Restrict the test_R2_avg to values above this threshold. Default = 0.")
parser.add_argument("--pred.perf.pval_thres", type = float, dest = "pred_perf_pval_thres", default = 1, help = "Restrict the pred_perf_pval to values below this threshold. Default = 1.")

args = parser.parse_args() #then pass these arguments to further things

###INPUT SANITATION
#List of dbs holds addresses to .db files user wants to query

if args.db is None:
    print("No .db destination detected. Please input a .db destination using the --db flag.")
    sys.exit(1)

###GENES, GENENAMES
if args.genes is None and args.genenames is None:
    print("No list of genes has been supplied with --genes or --genenames. All genes in the model(s) will be queried.")
elif args.genes is not None and args.genenames is not None:
    print("Please select an input for only genes or only genenames and not both.")
elif args.genes is not None:
    query_genes = list(np.loadtxt(args.genes, dtype = "str", ndmin = 1))
else:# args.genenames is not None:
    query_genes = list(np.loadtxt(args.genenames, dtype = "str", ndmin = 1))

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
#print(extra_flags)
  
###WEIGHTS
weights_query = False #default to not query weights
weights_flags = []
if args.rsid or args.varID or args.ref_allele or args.eff_allele or args.weight:
    weights_query = True #you will be querying the weights matrix
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
query_flags = extra_flags + weights_flags + sample_info_flags
if len(query_flags) == 0:
    print("No query flags have been passed. Program will output cv_R2_avg, rsid, and weights of all genes in the input models.")
else: 
    print("Flags queried: " + ", ".join(query_flags))

#making input .dbs into a list
args_db = args.db #.endswith doesn't like arguments
if args_db.endswith(".db"): #its a single .db file
    dbs = [(args_db)]
    print("Model queried: " + args_db.split("/")[-1]) #don't print the full path
else: #its (I assume) a folder
    if args_db.endswith("/"):
        folder_name = args_db
    else:
        folder_name = args_db + "/"
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
print(dbs)

#THRES FLAGS (for filtering later)
test_R2_avg_thres = args.test_R2_avg_thres
cv_R2_avg_thres = args.cv_R2_avg_thres
rho_avg_thres = args.rho_avg_thres
pred_perf_R2_thres = args.pred_perf_R2_thres
pred_perf_pval_thres = args.pred_perf_pval_thres

########
#CARLEE#
########

data = [] #List of lists .db files info to output for further pandas filtering and parsing
#Note: Are we supposed to add .db file name to data list for pandas? Talk to Shreya later
  #Angela: yes

genename = NA
rsid = NA
varID = NA
ref_allele = NA
eff_allele= NA
weight = NA
n_snps_in_model = NA
test_R2_avg = NA
cv_R2_avg = NA
rho_avg = NA
rho_zscore = NA
pred_perf_R2 = NA
pred_perf_pval = NA
n_samples = NA
population = NA
tissue = NA

#dbs is a list containing strings that are addresses of the .db files
for db in dbs: #dbs are in .dbs
    conn = sqlite3.connect(db) #We need to make a SQL connection to the database we are querying
    c = conn.cursor() #c connects to all the .db files
    c.execute("select genename from extra;") #Angela: why genename? We need gene
    for row in c:
        genename = row[0]
    
    c.execute("select [n.snps.in.model], test_R2_avg, cv_R2_avg, rho_avg, rho_zscore, [pred.perf.R2], [pred.perf.pval] from extra ;")
    for row in c:
        n_snps_in_model = row[0]
        test_R2_avg = row[1]
        cv_R2_avg = row[2]
        rho_avg = row[3]
        rho_zscore = row[4]
        pred_perf_R2 = row[5]
        pred_perf_pval = row[6]

    for gene in query_genes: #Angela: since we're doing all genes in Carlee's part this is unnecessary
        c.execute("select gene from extra;")
        for row in c:
            gene = row[0]
        #^Am I supposed to add this; Angela: we start with genes, not genenames
        
        c.execute("select rsid, varID, ref_allele, eff_allele, weight from weights ;")
        for row in c:
            rsid = row[0]
            varID = row[1]
            ref_allele = row[2]
            eff_allele = row[3]

    
    c.execute('select n_samples, population, tissue from sample_info ;')
    for row in c:
        n_samples = row[0]
        population = row[1]
        tissue = row[2]

#Angela: switch the order of sample info flags and weights flags and do this command during the iterations of weights flags
data.append([genename, gene, rsid, varID, ref_allele, eff_allele, weight, test_R2_avg, cv_R2_avg, rho_avg, rho_zscore, pred_perf_R2, pred_perf_pval, n_samples, population, tissue])
print(data)
conn.close()

########        
#SHREYA#
########

###PARSING QUERY OUTPUT
#alright so we're done querying shiz and we got a big boi list of lists
#convert list of lists into dataframe

#Don't we need to take in the genes that the user wants to query? Once we do that, I would add a line of code to the flag retrieving
#code to only keep the rows with those genes, correct?
  #Angela: yes, ask Carlee what she saved those as
data_frame = pd.DataFrame(data) #make list of lists into dataframe
data_frame.index = ["gene"]
data_frame.columns = ["rsid", "varID", "ref_allele", "eff_allele", "weight", "genename", "gene_type", "alpha",
               "n_snps_in_window", "n.snps.in.model", "lambda_min_mse", "test_R2_avg", "test_R2_sd", "cv_R2_avg",
               "cv_R2_sd", "in_sample_R2","nested_cv_fisher_pval", "rho_avg", "rho_se", "rho_zscore", "pred.perf.R2",
               "pred.perf.pval", "pred.perf.qval", "chromosome", "cv_seed", "n_samples", "population", "tissue"] #give column names so user knows what they're looking at

#Angela: this is already called query_flags so its redundant
user_specified_flags = []
num_of_flags = len(extra_flags) + len(weights_flags) + len(sample_info_flags)
for flag in extra_flags:
    user_specified_flags.append(flag)
for flag in weights_flags:
    user_specified_flags.append(flag)
for flag in sample_info_flags:
    user_specified_flags.append(flag)
    
#if user has flags 
  #get everything "true" the user wants
    #always include genename
    #ex. if they want cv_r2_avg, rho_avg, and rsid, make a list of ["genename", "cv_r2_avg", "rho_avg", "rsid"]
#if user has no flags
  #genename, cv_r2_avg, rs, and weights
#give column names
  #(this will all be in the same order so we just need to figure out what this order is)
#only pull columns of what the user wants

#restrict to thresholds the user wants (see threshold flags)
  #ex. if only want cv_r2_avg > 0.1, only keep those
   #(this will involve using a bunch of .loc)
#To do this I will use the query_weight_vals array created above (-Shreya)

if num_of_flags > 0:
    data_frame_mod = data_frame.loc[user_specified_flags]
else:
    data_frame_mod = data_frame.loc["genename", "cv_r2_avg", "rsid", "weights"]
  
#delete duplicate rows

#print to csv
  #remove indexes cause they're annoying
  
#additional ideas:
#I want to do this with something like Beautifulsoup4 (-Shreya)
  #pull what the genes have been implicated in the GWAS catalog - https://www.ebi.ac.uk/gwas/downloads
    
'''
NOTE TO READER FROM ANGELA:
We decided to query all genes in the .db and then subsequently parse because:
1. The Wheeler Lab machine has so much memory size isn't an issue
2. Carlee seemed more comfortable with hardcoding than working with shifting variables
3. It's easier to have two people work on stuff if the stuff remains the same instead of variable
4. I'm more comfortable in Pandas than SQL
'''
