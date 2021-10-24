#!/usr/local/bin/python3.7
# encoding: utf-8
'''


@author:     Simon Rayner

@copyright:  2020 Oslo University Hospital. All rights reserved.

@license:    license

@contact:    simon.rayner@medisin.uio.no
@deffield    updated: Updated
'''

import sys
import csv
import os
from pathlib import Path

from datetime import datetime
import hashlib
import logging

import pandas as pd

from argparse import ArgumentParser
from argparse import RawDescriptionHelpFormatter
from Bio.Data.CodonTable import list_possible_proteins


__all__ = []
__version__ = 0.1
__date__ = '2020-12-19'
__updated__ = '2020-12-19'

DEBUG = 1
TESTRUN = 0
PROFILE = 0

class CLIError(Exception):
    '''Generic exception to raise and log different fatal errors.'''
    def __init__(self, msg):
        super(CLIError).__init__(type(self))
        self.msg = "E: %s" % msg
    def __str__(self):
        return self.msg
    def __unicode__(self):
        return self.msg
    
    
def initLogger(md5string):
    
    ''' setup log file based on project name'''
    projectBaseName = os.path.splitext(os.path.basename(resultfiles))[0]
    now = datetime.now()
    dt_string = now.strftime("%Y%m%d_%H%M%S")
    logFolder = os.path.join(os.getcwd(), "logfiles")
    if not os.path.exists(logFolder):
        print("--log folder <" + logFolder + "> doesn't exist, creating")
        os.makedirs(logFolder)   
    logfileName = os.path.join(logFolder, projectBaseName + "__" + dt_string + "__" + md5string +".log")
    handler = logging.StreamHandler(sys.stdout)
    logging.basicConfig(level=logging.DEBUG)
    
    fileh = logging.FileHandler(logfileName, 'a')
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    fileh.setFormatter(formatter)
    
    log = logging.getLogger()  # root logger
    for hdlr in log.handlers[:]:  # remove all old handlers
        log.removeHandler(hdlr)
    log.addHandler(fileh)      # set the new handler 
    log.addHandler(handler) 
    logging.info("+" + "*"*78 + "+")   
    logging.info("project log file is <" + logfileName + ">")         
    logging.info("+" + "*"*78 + "+")   
    logging.debug("debug mode is on")
    

def parseArgs(argv):
    
    '''parse out Command line options.'''

    program_name = os.path.basename(sys.argv[0])
    program_version = "v%s" % __version__
    program_build_date = str(__updated__)
    program_version_message = '%%(prog)s %s (%s)' % (program_version, program_build_date)
    program_shortdesc = __import__('__main__').__doc__.split("\n")[1]
    program_license = '''%s
    i
      Created by Simon Rayner on %s.
      Copyright 2020 Oslo University Hospital. All rights reserved.
    
      Licensed under the Apache License 2.0
      http://www.apache.org/licenses/LICENSE-2.0
    
      Distributed on an "AS IS" basis without warranties
      or conditions of any kind, either express or implied.
    
    USAGE
    ''' % (program_shortdesc, str(__date__))

    try:
        # Setup argument parser
        parser = ArgumentParser(description=program_license, formatter_class=RawDescriptionHelpFormatter)
        parser.add_argument("-r", "--resultfiles", dest="resultfiles", action="store", help="list of result files to process [default: %(default)s]")
        parser.add_argument("-u", "--upregulated", dest="upregulated", action="store", help="list of upregulated proteins [default: %(default)s]")
        parser.add_argument("-d", "--downregulated", dest="downregulated", action="store", help="list of downregulated proteins [default: %(default)s]")
        parser.add_argument("-p", "--probability", dest="probability", action="store", help="cut-off for min probability [default: %(default)s]")
        parser.add_argument("-e", "--energy", dest="energy", action="store", help="cut-off for min binding energy [default: %(default)s]")

        # Process arguments
        args = parser.parse_args()

        global resultfiles 
        global upregulatedProtFile
        global downregulatedProtFile
        global probability
        global energy 
        
        resultfiles = args.resultfiles
        upregulatedProtFile = args.upregulated
        downregulatedProtFile = args.downregulated
        probability = args.probability
        energy = args.energy

        print("resultfilelist <" + resultfiles + ">")
        print("upregulated <" + upregulatedProtFile + ">")
        print("downregulated <" + downregulatedProtFile + ">")
        print("probability <" + probability + ">")
        print("free energy <" + energy + ">")

    except KeyboardInterrupt:
        ### handle keyboard interrupt ###
        return 0
    except Exception as e:
        print(e)
        if DEBUG or TESTRUN:
            raise(e)
        indent = len(program_name) * " "
        sys.stderr.write(program_name + ": " + repr(e) + "\n")
        sys.stderr.write(indent + "  for help use --help")
        return 2
    
def loadProteinData():

    global dfMiRsVsProteins
    global proteins
    logging.info("loading UP protein file <" + upregulatedProtFile + ">")
    with open(upregulatedProtFile) as fUP:
        upProteins = fUP.read().splitlines()
    logging.info("done")
    logging.info("loading DOWN protein file <" + upregulatedProtFile + ">")
    with open(downregulatedProtFile) as fDOWN:
        downProteins = fDOWN.read().splitlines()
        
    # build data frame, one row / protein
    proteins = upProteins + downProteins
    dfMiRsVsProteins = pd.DataFrame(proteins)
    dfMiRsVsProteins = dfMiRsVsProteins.set_index(dfMiRsVsProteins[0])
    miRNames=[]
    miRNames = miRNAs['miRName'].tolist()
    dfMiRsVsProteins=pd.concat([dfMiRsVsProteins,pd.DataFrame(columns=miRNames)])
    dfMiRsVsProteins = dfMiRsVsProteins.drop([0], axis=1)
    logging.info("done")
    
def loadSampleList():
    global targetFiles
    global miRNAs
    # isos/karlsen_isos.MIMAT0000087_hsa-miR-30a-5p_GTAAACA_iso/allTargetSites.csv
    logging.info("processing target files")
    targetFiles=[]
    miRNames=[]
    
    miRNAs=pd.read_csv(resultfiles, sep='\t')
    miRNAs['miRName']='unknown'
    for index, row in miRNAs.iterrows():
        miRNAs.at[index, 'miRName'] = row['CHANGE'] + "__" + str(row['INDEX']) + "__" + os.path.basename(Path(row['FILE']).parents[0]).split('.')[1]

        
    
    #with open(resultfiles, 'r') as fd:
    #    csvrTargetFiles = csv.reader(fd)
    #    for targetFile in csvrTargetFiles:
    #        targetFiles.append(targetFile[0])
    #        miRNames.append(os.path.basename(Path(targetFile[0]).parents[0]).split('.')[1])
    logging.info("found <" + str(len(miRNAs))  + "> samples")
        

def processSamples():
    # read target list
    # loop through target list
    global targetFiles
    logging.info("processing target files")
    for index, row in miRNAs.iterrows():
        logging.info("-- processing file <" + miRNAs.at[index, 'FILE'] + ">")
        filterPredictionFile(miRNAs.at[index, 'FILE'], miRNAs.at[index, 'miRName'])
             
    
    outputFile = os.path.splitext(resultfiles)[0] + ".tsv"       
    dfMiRsVsProteins.to_csv(outputFile, sep='\t')
             
def filterPredictionFile(predFile, miRName):
    
    # grab miRNA name from parent folder name
    
    # read target predictions and filter by probability and energy
    logging.info("miR <" + miRName + "> " )
    preds = pd.read_csv(predFile, sep='\t')
    logging.info("--read <" + str(len(preds)) + "> lines" )
    predsFilter = preds[(preds['FreeEnergy']<float(energy)) & (preds['Prediction']>float(probability))]
    logging.info("--after processing,  <" + str(len(predsFilter)) + "> lines remain" )
    
    # add up- and down-regulated proteins as columns to the dataframe

    for protein in proteins:
        #logging.info("----protein:" + protein)
        if len(predsFilter[predsFilter['GeneName'].str.contains(protein)==True].index) > 0:
            dfMiRsVsProteins[miRName][protein]=1
        else:
            dfMiRsVsProteins[miRName][protein]=0
    
    logging.info("done")
          
    
    
    
    

def main(argv=None): # IGNORE:C0111
    
    
    #if argv is None:
    #    argv = sys.argv

    md5String = hashlib.md5(b"CBGAMGOUS").hexdigest()
    parseArgs(argv)
    initLogger(md5String)
    loadSampleList()
    loadProteinData()
    processSamples()
    

    logging.info("finished")

if __name__ == "__main__":
    if DEBUG:
        pass
        #sys.argv.append("-h")
        #sys.argv.append("-v")

    if TESTRUN:
        import doctest
        doctest.testmod()
    if PROFILE:
        import cProfile
        import pstats
        profile_filename = 'fairpype.virusPipe_profile.txt'
        cProfile.run('main()', profile_filename)
        statsfile = open("profile_stats.txt", "wb")
        p = pstats.Stats(profile_filename, stream=statsfile)
        stats = p.strip_dirs().sort_stats('cumulative')
        stats.print_stats()
        statsfile.close()
        sys.exit(0)
    sys.exit(main())