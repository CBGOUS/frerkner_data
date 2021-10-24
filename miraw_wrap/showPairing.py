
#!/usr/bin/env python
"""This function is used to add pairing show 
in original miRAW output files with the following command:
    python showPairing.py -f filepath -p/-n/-a
the last three options could be a single one or multiple ones
-f to indicate the path of the files to be processed
-p positive targets file by miRAW (.positiveTargetSites.csv)
-n negative targets file by miRAW (.negativeTargetSites.csv)
-a all targets file by miRAW (.allTargetSites.csv)
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

#data frame of the output files by miRAW
#GeneName   miRNA   SiteStart   SiteEnd Prediction  PairStartinSite SeedStart   SeedEnd Pairs   WC  Wob MFE Comment SiteTranscript  MatureMiRNATranscript   BracketNotation AdditionalProperties
#0              1       2           3       4               5           6          7       8    9   10  11  12          13              14                      15              16


logging.getLogger().setLevel(logging.INFO)

parser = argparse.ArgumentParser(description='show pairing according to predictions')

parser.add_argument("-H", "--HelpMe", action="store_true",
                    help="print detailed help")

parser.add_argument("-f", "--target_site_file", dest='fileToProcess',
                    help="miRAW detailed results folder")

parser.add_argument("-p","--positiveTarget",action="store_true",
                    help="to show pairing in .positiveTargetSites.csv")

parser.add_argument("-n","--negativeTarget",action="store_true",
                    help="to show pairing in .negativeTargetSites.csv")

parser.add_argument("-a","--allTarget",action="store_true",
                    help="to show pairing in .allTargetSites.csv")

args = parser.parse_args()


def printLongHelpAndExit():
    logging.info("+" + "-" * 78 + "+")
    logging.info("+  addBracketNotation2PredParis                                                +")
    logging.info("+  you need to specify:                                                        +")
    logging.info("+                                                                              +")
    logging.info("+      a miRAW detailed results file path as input      (-f/--target_site_file)+")
    logging.info("+        (this is the file path, also the name of one case in the experiments) +")
    logging.info("+                                                                              +")
    logging.info("+      which file(s) will be processed    (-p/-n/-a)                           +")
    logging.info("+        -p, the file containing positive target sites,                        +")
    logging.info("+            with tail of .positiveTargetSites.csv                             +")
    logging.info("+                                                                              +")
    logging.info("+        -n, the file containing negative target sites,                        +")
    logging.info("+            with tail of .negativeTargetSites.csv                             +")
    logging.info("+                                                                              +")
    logging.info("+        -a, the file containing all predicted target sites,                   +")
    logging.info("+            with tail of .allTargetSites.csv                                  +")
    logging.info("+" + "-" * 78 + "+")


def printHelpAndExit():
    parser.print_help()
    logging.info("stopping") 
    sys.exit()


def checkTargetsFile(): #check the validity of the commanded files   
    global targetFiles, filename
    logging.info("check target file(s)")
    if args.fileToProcess:
        foldername=args.fileToProcess
        ind = foldername.rfind("/")
        filename = foldername + foldername[ind:len(foldername)]
        if args.positiveTarget:
            file_tail.append(".positiveTargetSites.csv")
        if args.negativeTarget:
            file_tail.append(".negativeTargetSites.csv")
        if args.allTarget:
            file_tail.append(".allTargetSites.csv")

        targetFiles=[filename + tail for tail in file_tail]
        for file in targetFiles:
            if not os.path.isfile(file):
                logging.error("--can't find sites file at <" + file + ">")
                exit()
    logging.info("--OK")       

      
def pairbyBracketNotation(): #read and process the commanded files 
    logging.info("process all target files")
    global filename, file_tail

    for filetail in file_tail:
        currentFile = filename + filetail
        with open(currentFile, 'r') as fin:
            reader = csv.reader(fin, delimiter='\t')
            with open(filename+".pairing"+filetail, 'w') as fout:
                writer = csv.writer(fout, delimiter='\t')
                # set headers here, grabbing headers from reader first
                head = reader.next()
                head.append('Pairing')
                writer.writerow(head)
                for row in reader:
                    #get the needed parameters
                    pair_start_in_site = 39-int(row[5])
                    utr = row[13]
                    utr_re = preprocessUtr(utr)
                    mirna = row[14]
                    bknotation = row[15]

                    #reform the bracket notation
                    bknotation_mirna = bknotation[len(utr)+1:len(bknotation)]

                    #remove self secondary structure in the 3'utr transcript and reverse the bracket notation
                    bknotation_utr = preprocessBN(bknotation[0:len(utr)])
                    bknotation_utr = alignAndExtendBN(bknotation_mirna,bknotation_utr,pair_start_in_site)
                                    
                    #show pairing using seed region bracket notation
                    new_value = pairing(utr_re,mirna,bknotation_utr,bknotation_mirna)
                    row.append(new_value)
                    writer.writerow(row)

    logging.info("--done")  

def preprocessUtr(utr):         #reverse the 3utr sequence and replace T by U, starting from the binding site
    logging.info(" preprocess 3utr transcript")

    re_utr=utr[::-1]
    processedUtr=re_utr.replace("T","U")

    logging.info("--done")
    return processedUtr

def preprocessBN(notation):     #remove secondary structure in the 3'utr transcript
    logging.info(" preprocess bracket notation")

    while notation.find(")")>0:
        listBN = list(notation)
        ind = notation.find(")")
        listBN[ind] = "."
        ind_re = notation[0:ind-1].rfind("(")
        listBN[ind_re] ="."
        notation = ''.join(listBN)

    logging.info("--done")
    return notation[::-1]    

def alignAndExtendBN(bnmirna,bnutr,pairstart): #reform the bracket notation according to the aligned mRNA transcript
    logging.info("reform bracket notation")

    ind_1st_pair = bnutr.find("(")
    if ind_1st_pair<pairstart:
        re_bnutr = "."*(pairstart-ind_1st_pair)+bnutr[0:(len(bnutr)-(pairstart-ind_1st_pair))]       
    else:
        re_bnutr = bnutr[(ind_1st_pair-pairstart) : len(bnutr)]+ "."*(ind_1st_pair-pairstart)
    logging.info("--done")
    return re_bnutr

def pairing(re_utr,mirnaseq,bnotation_utr,bnotation_mirna):
    logging.info("start the pairing process for each mirna:mrna")

    strline1 = mirnaseq
    strline2 = ""
    strline3 = re_utr
    listBNmirna = list(bnotation_mirna)
    listBNutr = list(bnotation_utr)
    for ind, item in enumerate(listBNutr):
        if ind < len(listBNmirna):
            if item == listBNmirna[ind]:
                strline2 += " "
            elif item == "." and listBNmirna[ind] == ")":
                listBNmirna.insert(ind,"-")
                if ind<bnotation_utr.find("("):
                    strline1 = " "+strline1
                else:
                    strline1 = strline1[0:ind]+"-"+strline1[ind:len(strline1)]
                strline2 += " "
            elif item == "(" and listBNmirna[ind] == ")":   
                strline2 += "|" 
            elif item == "(" and listBNmirna[ind] == ".":
                strline2 += " "
                strline3 = strline3[0:ind]+"-"+strline3[ind:len(strline3)]  
    str_pairing = "5' "+strline1+"  3' miRNA"+"\r\n"+"   "+strline2 +"\r\n"+"3' "+strline3+"  5' mRNA"

    logging.info("--done")
    return str_pairing
            
            
               
def checkArgs():
    
    if args.HelpMe :
        printLongHelpAndExit() 
        exit()       
    checkTargetsFile()

    
checkArgs()
pairbyBracketNotation()



