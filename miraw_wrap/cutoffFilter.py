
#!/usr/bin/env python
"""This function is used to 
filter miRAW results by MFE and/or Prediction 
with the following command:

    python cutoffFilter.py -f filepath/file -P value -E value

Currently, it only supports prediction probability and free energy filtering
-P the cutoff value for prediction probability filtering
-E the cutoff value for free energy filtering

*The values could be given either in negative or positive
*The function only compares the absolute value
*The kept absolute values for each prediction are larger than the given absolute cutoff value
*Users need to make sure the file they would like to apply the filter
"""

import argparse
import sys
import os
import logging
import csv
import datetime
import re

__author__ = "Yafei Xing"
__copyright__ = "Copyright 2018, AMG-OUS"
__version__ = "1.0.1"
__maintainer__ = "Yafei Xing"
__email__ = "yafei.xing@medisin.uio.no"
__status__ = "Production"

#global variables
targetFiles = ""
filename =""
file_tail = []
EnergyCutoff = 0
ProbabilityCutoff = 0
#data frame of the output files by miRAW
#GeneName   miRNA   SiteStart   SiteEnd Prediction  PairStartinSite SeedStart   SeedEnd Pairs   WC  Wob MFE Comment SiteTranscript  MatureMiRNATranscript   BracketNotation AdditionalProperties
#0              1       2           3       4               5           6          7       8    9   10  11  12          13              14                      15              16


logging.getLogger().setLevel(logging.INFO)

parser = argparse.ArgumentParser(description='show pairing according to predictions')

parser.add_argument("-H", "--HelpMe", action="store_true",
                    help="print detailed help")

parser.add_argument("-f", "--target_site_file", dest='fileToProcess',
                    help="miRAW detailed results folder")

parser.add_argument("-P", "--predProb", dest='probCutOff',
                    help="the cutoff value for prediction probability filtering")

parser.add_argument("-E", "--Energy", dest='energyCutOff',
                    help="the cutoff value for free energy filtering")
args = parser.parse_args()


def printLongHelpAndExit():
    logging.info("+" + "-" * 78 + "+")
    logging.info("+  To use cutoffFilter.py                                                      +")
    logging.info("+  you need to specify:                                                        +")
    logging.info("+                                                                              +")
    logging.info("+      a miRAW detailed results file as input           (-f/--target_site_file)+")
    logging.info("+                                                                              +")
    logging.info("+      the cuttoff value for prediction probability filtering   (-P/--predProb)+")
    logging.info("+                                                                              +")
    logging.info("+      the cuttoff value for free energy filtering             (-E/--Energy)   +")
    logging.info("+                                                                              +")
    logging.info("+" + "-" * 78 + "+")


def printHelpAndExit():
    parser.print_help()
    logging.info("stopping") 
    sys.exit()


def checkTargetsFile(): #check the validity of the commanded files   
    global foldername, targetFiles
    logging.info("check target file")
    if args.fileToProcess:
        foldername=args.fileToProcess
        ind = foldername.rfind("/")
        filename = foldername[ind+1:len(foldername)]
        ind_dot = foldername.rfind(".")
        file_tail = foldername[ind_dot+1:len(foldername)]
        targetFiles=foldername[0:ind_dot]+".cutoffFiltered."+file_tail
        if not os.path.isfile(foldername):
            logging.error("--can't find sites file at <" + file + ">")
            exit()
    logging.info("--OK")  


def checkBindingEnergyCutOff():
    logging.info("checkBindingEnergyCutOff") 
    
    global EnergyCutoff  
    if args.energyCutOff:
        logging.info(args.energyCutOff)
        if float(args.energyCutOff) <= 0.0:
            EnergyCutoff=abs(float(args.energyCutOff))
        else:
            logging.info(" binding energy cutoff (-e) needs to be less than 0" + str(args.energyCutOff) + ">")
            return()
    logging.info("--OK") 


def checkProbCutOff():
    logging.info("checkProbabilityCutOff") 
    
    global ProbabilityCutoff  
    if args.probCutOff:
        logging.info(args.probCutOff)
        if( abs(float(args.probCutOff)) >= 0.0 and abs(float(args.probCutOff)) <= 1.0):
            ProbabilityCutoff=abs(float(args.probCutOff))
        else:
            logging.info(" the probability cutoff (-P) needs to be between 0 and 1" + str(args.probCutOff) + ">")
            return()
    logging.info("--OK") 


def cutoffFiltering(col,cutoff):
    logging.info("Cutoff process")

    global foldername, targetFiles
    with open(foldername, 'r') as fin:
        reader = csv.reader(fin, delimiter='\t')
        with open(targetFiles, 'w') as fout:
            writer = csv.writer(fout, delimiter='\t')
            head = reader.next()
            writer.writerow(head)
            writer.writerows([row for row in reader if all(abs(float(row[col[index]]))>=cutoff[index] for index in range(0,len(col)))])
    logging.info("--done")  


def checkArgs():
    col2Process = []
    cutOff = []
    global EnergyCutoff, ProbabilityCutoff
    if args.HelpMe :
        printLongHelpAndExit() 
        exit()       
    checkTargetsFile()
    checkProbCutOff()
    if ProbabilityCutoff!= 0:
        col2Process.append(4)
        cutOff.append(ProbabilityCutoff)
    checkBindingEnergyCutOff()
    if EnergyCutoff!= 0:
        col2Process.append(11)
        cutOff.append(EnergyCutoff)
    cutoffFiltering(col2Process,cutOff)

checkArgs()