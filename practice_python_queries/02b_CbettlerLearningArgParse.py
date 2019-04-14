#load in list of gene names you want to query
#APOE and PSRC
#^Learning argparse to help accomplish this 

import sqlite3
import pandas
import os
from Bio import SeqIO
import argparse
pops = ["AFA", "AFHI", "ALL", "CAU", "HIS"]

currentDir = os.popen('pwd').read().rstrip()
#Gets path of current working directory to avoid hardcoding

db_path = (currentDir+ "/gtex_v7_Whole_Blood_imputed_europeans_tw_0.5_signif")
#Must be in current directory of db file. Will make this less hard-codey in the future.

#Running argparse tutorials to learn how to use it
'''
parser = argparse.ArgumentParser()
parser.add_argument("echo", help="echo the string you use here")
args = parser.parse_args()
print(args.echo)
#Note need to type echo after running
'''

#Running argparse tutorials to learn how to use it
parser = argparse.ArgumentParser(description='Short sample app')
parser.add_argument('-a', action="store_true", default=False)
parser.add_argument('-b', action="store", dest="b")
parser.add_argument('-c', action="store", dest="c", type=int)
print(parser.parse_args(['-a', '-bval', '-c', '3']))
#This runs
#Shreya try running this through the command line using a .py file


'''
currentDir = os.popen('pwd').read().rstrip()
geneFileName = input("What is your gene file exact name? Note that is must be a .txt file")
parser = argparse.ArgumentParser()
parser.add_argument(currentDir + "/"+geneFileName, type=argparse.FileType('r'))
#File is sample gene name file
args = parser.parse_args()
print(args.file.readlines())
'''
#^Trying to apply arg parse to parsing a file of gene names
#Still trying to figure out why this won't run
#Attempting to avoid hardcoding the file path, still not sure if this is the best way to do this

#Good video to learn argparse: https://www.youtube.com/watch?v=cdblJqEUDNo
#Note that this code does not run properly, this is rather evidence of me trying to learn a new skill and apply it to our data. 
