#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""This function is used to generate a batch shell script 
for adding pairing show in a series of original miRAW output files.
This code is used in the follwoing format

    python filepath1/extractConflictBatch.py -f filepath2/extractConflicts.py -o miRAWexperimentpath -e batchfilename 
    
the last three options could be a single one or multiple ones
-o to indicate the miraw experiment path which contains all mirna:mrna predictions
-f to indicate which function will be applied to miRAW predictions
-e to name the generated batch file
"""

import argparse
import sys
import os
import logging
import datetime
from Bio import SeqIO


__author__ = "Yafei Xing"
__copyright__ = "Copyright 2018, AMG-OUS"
__version__ = "1.0.1"
__maintainer__ = "Yafei Xing"
__email__ = "yafei.xing@medisin.uio.no"
__status__ = "Production"


MY_NEWLINE = "\n"
if os.name== "Windows":
    MY_NEWLINE ="\r\n"


logging.getLogger().setLevel(logging.INFO)


parser = argparse.ArgumentParser(description='set up batch commands to run seed pairing')

parser.add_argument("-e", "--exptName", dest='exptName',
                    help="specify name for this process")

parser.add_argument("-o", "--outFolder", dest='outFolder',
                    help="folder path for miRAW output and new output")

parser.add_argument("-f", "--funcLoc", dest='functionLocation',
                    help="location of the showPairing.py file")

parser.add_argument("-H", "--HelpMe", action="store_true",
                    help="print detailed help")

args = parser.parse_args()


def printLongHelpAndExit():
    logging.info("+" + "-" * 78 + "+")
    logging.info("+  PairBatch:                                                                  +")
    logging.info("+    wrapper code for generating files and scripts to run showPairing.py       +")
    logging.info("+                                                                              +")
    logging.info("+    you need to specify:                                                      +")
    logging.info("+                                                                              +")
    logging.info("+      a name for the batch file: (-e/--exptName)                              +")
    logging.info("+                                                                              +")
    logging.info("+      an OutputFolder to write results:           (-o/--outFolder)            +")
    logging.info("+                                                                              +")
    logging.info("+      the location of the extractConflicts.py file          (-f/--functionLoc)+")
    logging.info("+" + "-" * 78 + "+")


def printHelpAndExit():
    parser.print_help()
    logging.info("stopping")
    sys.exit()


def checkProgramLocation():
    #python functionpath/showPairing.py folderpath/folderpath.positiveTargetSites.csv
    global funcLocation
    if args.functionLocation:
        funcLocation =  args.functionLocation
    logging.info("program location is at <" + funcLocation + ">")


def checkExptName():
    
    logging.info("checking Name for the batch file:" )
    
    if args.exptName:
        logging.info(args.exptName)
        logging.info("--OK")
    else:
        logging.error("----you need to specify a name for the script using the -e/--exptName parameter")
        printHelpAndExit()


def checkOutFolder():
    
    logging.info("checking outFolder:")
    scriptFolder = os.path.dirname(os.path.abspath(__file__))#current absolute path
    
    if args.outFolder:
        args.outFolder=args.outFolder.strip()
        # is the file path a relative path (does it start with a '.' ?)
        if not os.path.isdir(args.outFolder):
            args.outFolder = os.path.join(scriptFolder, args.outFolder)
            logging.info("--OK")
        
        if not os.path.isdir(args.outFolder):
            logging.error("--the folder doesn't exist, can't continue. Try checking the specified path")
            logging.info("--stopping")
            sys.exit()           
    else:
        logging.error("----you need to specify an output folder using the -o/--outFolder parameter")
        printHelpAndExit()

        # check folder exists or is creatable


def checkArgs():
    
    if args.HelpMe :
        printLongHelpAndExit()
        exit()
    checkProgramLocation()
    checkExptName()
    checkOutFolder()


def writeScript():

    #python functionpath/extractConflicts.py -f folderpath/folderpath
    with open(os.path.join(args.outFolder, args.exptName  + '.sh'), "w") as f:
        for foldername in next(os.walk(args.outFolder))[1]:
            folderpath = os.path.join(args.outFolder, foldername)
            content = "python " + funcLocation + " -f " + folderpath
            f.write(content + MY_NEWLINE)



checkArgs()
writeScript()



