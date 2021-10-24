
#!/usr/bin/env python
"""This function is used to remove conflicts from miRAW
predictions and save all the conflicted observations in 
a separated file with the following command:
    python extractConflicts.py -f filepath 
-f to indicate the path of the files to be processed
The function will generate four additional files:
.allTargetSites.withoutConflicts.csv
.allTargetSites.onlyConflicts.csv
.positiveTargetSites.withoutConflicts.csv
.negativeTargetSites.withoutConflicts.csv

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
targetFiles = [] #a file list containing all the files needed to be processed
file_tail = [".targetPredictionOutput",".allTargetSites", ".positiveTargetSites", ".negativeTargetSites"]
genelist = []
mirnalist = []
conflictCount = 0

#data frame of the output files by miRAW
#GeneName   miRNA   SiteStart   SiteEnd Prediction  PairStartinSite SeedStart   SeedEnd Pairs   WC  Wob MFE Comment SiteTranscript  MatureMiRNATranscript   BracketNotation AdditionalProperties
#0              1       2           3       4               5           6          7       8    9   10  11  12          13              14                      15              16


logging.getLogger().setLevel(logging.INFO)

parser = argparse.ArgumentParser(description='remove conflicts in miRAW predictions')

parser.add_argument("-H", "--HelpMe", action="store_true",
                    help="print detailed help")

parser.add_argument("-f", "--folder_name", dest='folderName',
                    help="miRAW detailed results folder")


args = parser.parse_args()


def printLongHelpAndExit():
    logging.info("+" + "-" * 78 + "+")
    logging.info("+  addBracketNotation2PredParis                                                +")
    logging.info("+  you need to specify:                                                        +")
    logging.info("+                                                                              +")
    logging.info("+      a miRAW detailed results file path as input      (-f/--folder_name)     +")
    logging.info("+        (this is the file path, also the name of one case in the experiments) +")
    logging.info("+" + "-" * 78 + "+")


def printHelpAndExit():
    parser.print_help()
    logging.info("stopping") 
    sys.exit()


def checkTargetsFile(): #check the validity of the commanded files   
    global targetFiles, path_filename
    logging.info("check target file(s)")
    if args.folderName:
        foldername = args.folderName
        ind = foldername.rfind("/")
        filename = foldername[ind:len(foldername)]
        path_filename = foldername + filename
        targetFiles = [path_filename + tail + '.csv' for tail in file_tail]
        for file in targetFiles:
            if not os.path.isfile(file):
                logging.error("--can't find sites file at <" + file + ">")
                exit()
    logging.info("--OK")       


def extractConflicts(): #extract the conflicts
    global targetFiles
    conflictExist = checkSummary()
    if conflictCount>0:
        processFiles()
    else:
        logging.info("--there is no conflicts in this case")
        exit()


def checkSummary(): #check the prediction summary file to extract conflicted observations
    #GeneName   GeneId  miRNA   Prediction  HighestPredVal  LowestPredVal   PosSites    NegSites    RemovedSites    
    # 0            1        2       3           4               5               6           7           8
    logging.info("check the summary file")
    global conflictCount, genelist, mirnalist, targetFiles
    with open(targetFiles[0], 'r') as fin:
        reader = csv.reader(fin, delimiter='\t')
        reader.next()
        for row in reader:
            if int(row[6])>0 and int(row[7])>0:
                conflictCount = conflictCount+1
                genelist.append(row[0])
                mirnalist.append(row[2])
    logging.info("--done") 


def processFiles(): #remove conflicts in the prediction target files and save them separately
#data frame of the output files by miRAW
#GeneName   miRNA   SiteStart   SiteEnd Prediction  PairStartinSite SeedStart   SeedEnd Pairs   WC  Wob MFE Comment SiteTranscript  MatureMiRNATranscript   BracketNotation AdditionalProperties
#0              1       2           3       4               5           6          7       8    9   10  11  12          13              14                      15              16
    logging.info("process all target files")
    global targetFiles, genelist, mirnalist, path_filename, conflictCount
    for i in range(1,4):
        with open(targetFiles[i], 'r') as fin:
            reader = csv.reader(fin, delimiter='\t')
            head = reader.next()
            filename = path_filename+file_tail[i]+".withoutConflicts.csv"
            fout = open(filename, 'w')
            writer = csv.writer(fout, delimiter='\t')               
            writer.writerow(head)
            if i == 1: #for saving the conflicts
                filename_add = path_filename+file_tail[i]+".onlyConflicts.csv"
                f_add = open(filename_add, 'w')
                writer_add = csv.writer(f_add, delimiter='\t')
                writer_add.writerow(head)
            
            for row in reader:
                exist = 0
                for ind in range(0,conflictCount):
                    if row[0] == genelist[ind] and row[1] == mirnalist[ind]:   
                        exist = 1
                        break
                if exist == 0:
                    writer.writerow(row)
                elif exist == 1 and i == 1:
                    writer_add.writerow(row)

            f_add.close()
            fout.close()
    logging.info("--done")

               
def checkArgs():
    
    if args.HelpMe :
        printLongHelpAndExit() 
        exit()       
    checkTargetsFile()

    
checkArgs()
extractConflicts()



