import sqlite3
import pandas as pd

#which .db are we querying? SQLite3
#.db file is in Angela's folder
path = "/home/aandaleon/SQLite3-wrapper/db"
#open connection to db file
conn_to_db = sqlite3.connect(path + "gtex_v7_Whole_Blood_imputed_europeans_tw_0.5_signif.db") 
#start a list of lists to store data in (only have one dataframe to store this all in)
geneName_list = []
#load in list of gene names you want to query
  #APOE and PSRC1
for line in open(path + "gene"): #find gene from gene_name
    geneName_list.append(line.strip().split())
  #iterate through list
for gene in geneName_list:
    print gene
#what are the names of these genes (starting w/ ENSG)?
#what's its test R2 average?
#here we create a cursor to establish a connection to the SQLite db
curs = conn_to_db.cursor() 
curs.execute('SELECT (`gene` || " " || `genename` || " " || `pred.perf.r2`) FROM EXTRA;') #run search through python; "extra" database stores gene by gene information 
for row in curs: #iterate through the output of the query
  row = row[0].split(" ") #tuple is just one element, so split into two; happens b/c of the ` || " " || ` in the query
  genename_R2.append([row[0], row[1], row[2], pop]) #extract gene, gene name, and R2 from the row tuple and make into list to append to list of lists
conn_to_db.close() #close everything you open
#what SNPs are in the file
 #what are their rsids and weights?
 
#make list of lists into dataframe
#give column names so user knows what they're looking at
#write to .csv
