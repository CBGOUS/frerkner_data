
#!/usr/bin/env python
"""This function is used to
show if there is a binding at certain position(s) of 3'utr transcript
in original miRAW output files with the following command:

    python showPairing.py -f filepath/file -p1 position1

-f to indicate the abosulte path of the files to be processed
-p1 to indicate where a binding needs to be looked for

*Pay attention to parameter -p1, it needs to be the position counting from 0
*This version only supports for searching at one position
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

parser.add_argument("-p1","--pos1",dest="position_1",
                    help="to indicate the first position in 3'utr transcript")

#parser.add_argument("-p2","--pos2",dest="position_2",
#                    help="to indicate the second position in 3'utr transcript")

#parser.add_argument("-p3","--pos3",dest="position_3",
#                    help="to indicate the 3rd position in 3'utr transcript")

args = parser.parse_args()


def printLongHelpAndExit():
    logging.info("+" + "-" * 78 + "+")
    logging.info("+  showBindingAt.py                                                              +")
    logging.info("+  you need to specify:                                                        +")
    logging.info("+                                                                              +")
    logging.info("+      a miRAW detailed results file path as input      (-f/--target_site_file)+")
    logging.info("+        (this is the file path, also the name of one case in the experiments) +")
    logging.info("+                                                                              +")
    logging.info("+      which position will be processed    (-p1/--pos1)                        +")
    logging.info("+        -p1, the location to find a binding                                   +")
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
        
        if not os.path.isfile(foldername):
            logging.error("--can't find sites file at <" + file + ">")
            exit()

        targetFiles=foldername[0:ind_dot]+".BindingAt."+file_tail
    logging.info("--OK") 

def checkPosition():
    logging.info("check given position")
    if args.position_1:
        logging.info("--OK") 
    else:
        logging.error("please indicate the location")      

      
def procPair(): #process the positions

    logging.info("looking for bindings at given location")
    global foldername, targetFiles

    with open(foldername, 'r') as fin:
        reader = csv.reader(fin, delimiter='\t')
        with open(targetFiles, 'w') as fout:
            writer = csv.writer(fout, delimiter='\t')
            # set headers here, grabbing headers from reader first
            head = reader.next()
            head.append('BindingAtPos_'+args.position_1)
            writer.writerow(head)
            for row in reader:
                #get the needed parameters
                siteStart = int(row[2])
                siteEnd = int(row[3])-1
                if int(args.position_1)>=siteStart and int(args.position_1) <= siteEnd:                    
                    pairingStr = row[19]
                    ori_pos = int(args.position_1)
                    new_pairingStr = procString(pairingStr,siteStart,siteEnd,ori_pos)
                    row.append(new_pairingStr)
                else:
                    row.append('')
                writer.writerow(row)                                   

    logging.info("--done") 


def processPos(start, end, pos):

    if (end-start) < 39:
        if start == 0:
            new_pos =  end - pos 
        else:
            new_pos = 39 -(pos-start) 
    else:
        new_pos = end - pos 

    return new_pos


def procPosInStr(pos,strline):     #remove secondary structure in the 3'utr transcript
    logging.info("find the location in the string")

    pos_in_str = 3+pos
    count_bulg = strline[0:pos_in_str+1].count("-")
    if count_bulg!=0:
        pos_in_str=+count_bulg
        while strline[0:pos_in_str+1].count("-")> count_bulg:
            num = strline[0:pos_in_str+1].count("-")
            pos_in_str=pos_in_str+(num-count_bulg)
            count_bulg=num            

    logging.info("--done")
    return pos_in_str    


def procString(ori_str,start, end, pos):

    logging.info("process the string")
    ind = ori_str.rfind("\n")
    strline3 = ori_str[ind+1:len(ori_str)] #3utr
    reverse_pos = processPos(start,end,pos)
    pos_in_strline3 = procPosInStr(reverse_pos,strline3)
    ind2 = ori_str.find("\n")
    strline1 = ori_str[0:ind2+1]
    strline2 = ori_str[ind2+1:ind+1]
    if strline2[pos_in_strline3]=='|':
        strline2 = strline2[0:pos_in_strline3]+'{'+strline2[pos_in_strline3+1:len(strline2)]
        new_str = strline1+strline2+strline3
    else:
        new_str = ''
    logging.info("--done") 
    return new_str

               
def checkArgs():
    
    if args.HelpMe :
        printLongHelpAndExit() 
        exit()       
    checkTargetsFile()
    checkPosition()

    
checkArgs()
procPair()



