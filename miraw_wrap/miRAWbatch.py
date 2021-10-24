# -*- coding: utf-8 -*-
import argparse
import sys
import os
import logging
import datetime
from Bio import SeqIO

buildUnifiedFile = False

UNIFIEDFILE_EXISTS = 0
WRITE_BY_MIRNA = 1
WRITE_BY_UTR = 2
SPLIT_BY_MIRNA = 3
SPLIT_BY_UTR = 4
splitType = UNIFIEDFILE_EXISTS

CSSM_TS = "targetScan"
CSSM_PITA = "pita"
CSSM_REGULAR = "Regular"
CSSM_PERSONALIZED = "Personalized"

maximumSiteLength = 40
seedAlignmentOffset = 5
filterByAccessibilityEnergy = False
#I'm not clear how to specify this string, so for now we hard code as  
#it is specified in existing .properties files 
FILTER_ACCESSIBILITY_ENERGY_STRING="Filter0.AccessibilityEnergy;11" 

jarLocation = "miRAW.jar"  
unifiedFilePath = ""

miRheaderList = []
miRseqList = []
UTRheaderList = []
UTRseqList = []


splitData = False
MY_NEWLINE = "\n"
if os.name== "Windows":
    MY_NEWLINE ="\r\n"


logging.getLogger().setLevel(logging.INFO)


parser = argparse.ArgumentParser(description='set up batch commands to run miRAW')

parser.add_argument("-e", "--exptName", dest='exptName',
                    help="specify name for the analysis")

parser.add_argument("-o", "--outFolder", dest='outFolder',
                    help="output folder for results")

parser.add_argument("-u", "--unifiedFile", dest='unifiedFile',
                    help="absolute path for unifiedFile")

parser.add_argument("-d", "--dlModel", dest='dlModel',
                    help="absolute path to trained model file")

parser.add_argument("-c", "--cssm", dest='cssm',
                    help="Candidate Site Selection Model")

parser.add_argument("-m", "--mirnaFile", dest="miRFile",
                    help="absolute path to miRBase fastA file")

parser.add_argument("-3", "--3UTR", dest="utrFile",
                    help="absolute path to 3'UTR fastA file")

parser.add_argument("-t", "--split", dest="splitType",
                    help="how to split the unifiedFile")

parser.add_argument("-j", "--jarLoc", dest="jarFileLocation",
                    help="absolute path to the miRAW.jar file")

parser.add_argument("-H", "--HelpMe", action="store_true",
                    help="print detailed help")

parser.add_argument("-W", "--maximum_site_length", dest="maximumSiteLength",
                    help="maximum nucleotide length of a 3'UTR target")

parser.add_argument("-S", "--seed_alignment_offset", dest="seedAlignmentOffset",
                    help="stepwise distance for sliding target window along 3'UTR")

parser.add_argument("-E", "--energy_filtering", action="store_true",
                    help="perform energy filtering on results")

args = parser.parse_args()


def printLongHelpAndExit():
    logging.info("+" + "-" * 78 + "+")
    logging.info("+  RAWWrap:                                                                    +")
    logging.info("+    wrapper code for generating files and scripts to run miRAW                +")
    logging.info("+                                                                              +")
    logging.info("+    you need to specify:                                                      +")
    logging.info("+                                                                              +")
    logging.info("+      an ExperimentName for the .properties file: (-e/--exptName)             +")
    logging.info("+        (result files will also have this name)                               +")
    logging.info("+      an OutputFolder to write results:           (-o/--outFolder)            +")
    logging.info("+        (this can be a new or existing folder)                                +")
    logging.info("+      a DeepLearning model:                       (-d/--dlModel)              +")
    logging.info("+        (unless you want to build your own, use one of ours)                  +")
    logging.info("+        (see the README.md on how to build your own model)                    +")
    logging.info("+      a CSSM:                                     (-c/--cssm)                 +")
    logging.info("+        this can be Regular, Pita, targetScan or Personalized                 +")
    logging.info("+      a path to a unifiedFile:                    (-u/--unifiedFile)          +")
    logging.info("+        this can be an existing file or you can generate one by               +")
    logging.info("+        specifying both a miRNA                   (-m/--miRFile)              +")
    logging.info("+        and a 3'UTR file                          (-3/--3UTR)                 +")
    logging.info("+        (both these files should be in fasta format)                          +")
    logging.info("+      the location of the miRAW.jar file          (-j/--jarFileLocation)      +")
    logging.info("+        (if you want to run the predictions from a shell script)              +")
    logging.info("+                                                                              +")
    logging.info("+      NOTE: the file paths must be absolute, NOT relative                     +")
    logging.info("+                                                                              +")
    logging.info("+      You can also specify:                                                   +")
    logging.info("+      - how to handle the unifiedFile.                                        +")
    logging.info("+        If the file is large (i.e. requires many predictions) you may         +")
    logging.info("+        want to split into subfiles for running in parallel                   +")
    logging.info("+        you can do this using the -t/--split option                           +")
    logging.info("+                                                                              +")
    logging.info("+           Note: using -t 0   means you have an existing UnifiedFile you      +")
    logging.info("+                              want to use                                     +")
    logging.info("+                       -t 1-4 will generate one or more new UnifiedFiles      +")
    logging.info("+                              from the specified miRNA and 3'UTR input files  +")
    logging.info("+                                                                              +")
    logging.info("+        Options are:                                                          +")
    logging.info("+           0 : no split       (use (existing) specified unified File)         +")
    logging.info("+           1 : write by miRNA (single file for all predictions)               +")
    logging.info("+           2 : write by 3'UTR (single file for all predictions)               +")    
    logging.info("+           3 : split by miRNA (separate file for each miRNA)                  +")
    logging.info("+           4 : split by 3'UTR (separate file for each 3'UTR)                  +")    
    logging.info("+                                                                              +")
    logging.info("+        the program will output one or more shell scripts (to the specified   +")
    logging.info("+        output folder) containing the required commands to run miRAW          +")
    logging.info("+                                                                              +")
    logging.info("+      - energy filtering                          (-E/--energy_filtering)     +")
    logging.info("+        for conservative models (i.e. pita and targetScan) performance        +")
    logging.info("+        improves in the absence of energy filtering                           +")
    logging.info("+        for miRAW models using extended seed regions, performance improves    +")
    logging.info("+        when energy filtering is used                                         +")
    logging.info("+                                                                              +")
    logging.info("+          allowed values are (True/False)                                     +")
    logging.info("+                                                                              +")
    logging.info("+      - target site size (in nucleotides)         (--max_site_length)         +")
    logging.info("+          default 40nt                                                        +")
    logging.info("+      - step_size (in nucleotides)                (--seed_alignment_offset)   +")
    logging.info("+          default 5nt                                                         +")
    logging.info("+                                                                              +")
    logging.info("+        miRAW uses a sliding window to analyze a 3'UTR, so the                +")
    logging.info("+        default values will use a sliding window of 40nt with a step          +")
    logging.info("+        size of 5nt. Generally, there isn't much advantage in changing        +")
    logging.info("+        these values                                                          +")
    logging.info("+                                                                              +")
    logging.info("+                                                                              +")
    logging.info("+        simon.rayner@medisin.uio.no                                           +")
    logging.info("+                                                                              +")
    logging.info("+        you can also contact Albert Pla Planus (plaalbert@gmail.com)          +")
    logging.info("+        for questions about miRAW                                             +")
    logging.info("+                                                                              +")
    logging.info("+                                                                              +")
    logging.info("+" + "-" * 78 + "+")

    

def printHelpAndExit():
    parser.print_help()
    logging.info("stopping") 
    sys.exit()
    
    
    
def checkProgramLocation():
    #java -jar miRAW.jar GenePrediction predict ./Results/firstTry/firstTry.EF.properties
    global jarLocation        
    if args.jarFileLocation:
        jarLocation =  args.jarFileLocation
    logging.info("program location is at <" + jarLocation + ">")


def checkExptName():
    
    logging.info("checking Experimental Name:" )
    # get working folder
    
    if args.exptName:
        logging.info(args.exptName)
        logging.info("--OK")
    else:
        logging.error("----you need to specify an experimental Name using the -e/--exptName parameter") 
        printHelpAndExit()        
        
        
        
def checkOutFolder():  
          
    logging.info("checking outFolder:")
    scriptFolder = os.path.dirname(os.path.abspath(__file__))
        
    if args.outFolder:
        args.outFolder=args.outFolder.strip()
        # is the file path a relative path (does it start with a '.' ?)
        if not os.path.isdir(args.outFolder):
            os.makedirs(args.outFolder)
        
        if os.path.isabs(args.outFolder):
            logging.info("--created relative folder <" + args.outFolder + ">")
        else:
            logging.info("--created folder <" + os.path.join(scriptFolder, args.outFolder) + ">")    

        if not os.path.isdir(args.outFolder):
            logging.error("--create folder failed, can't continue. Try checking the specified path")
            logging.info("--stopping") 
            sys.exit()            
            
        args.outFolder = os.path.join(scriptFolder, args.outFolder)
        logging.info("--OK")
            
    else:
        logging.error("----you need to specify an output folder using the -o/--outFolder parameter") 
        printHelpAndExit()        
        
    # check folder exists or is creatable
    


    
def checkUnifiedFile():    
    # either the file exists, or we have to build it
    # if both the -m and -3 options are specified, we build a new one
    global buildUnifiedFile
    logging.info("checking UnifiedFile:")
    if args.unifiedFile:
        if args.miRFile and args.utrFile:
            buildUnifiedFile = True
            logging.info("--build unifiedFile")
            # check miR and 3'UTR files exist
            logging.info("--do all files exist?")
            if os.path.isfile(args.miRFile):
                if not os.path.isfile(args.utrFile):
                    logging.error("----can't find 3'UTR file at <" + args.utrFile + ">")
                    exit()
                else:
                    logging.info("----found all files")  
                    splitData = True              
            else:               
                logging.error("----can't find miRNA file at <" + args.miRFile + ">")
                exit()
            
        else:
            if os.path.isfile(args.miRFile) or os.path.isfile(args.utrFile):
                logging.info("--missing miR or 3'UTR file")
                logging.error("----you must specify both a miRFile and a 3'UTR file")
                printHelpAndExit()    
            else:
                logging.info("--using existing unifiedFile")
                        
    else:
        logging.info("----you need to specify a unified file using the -u/--unifiedFile parameter") 
        printHelpAndExit()        
    logging.info("--OK")
                
   
   
def checkCSSM():
    logging.info("checking CSSM:")
    if args.cssm:
        logging.info(args.cssm)
        if(args.cssm == CSSM_PITA):
            logging.info("--found pita CSSM")
            return()
        
        if(args.cssm == CSSM_TS):
            logging.info("--found targetScan CSSM")
            return()
        
        if(args.cssm == CSSM_REGULAR):
            logging.info("--found regular CSSM")
            return()
        
        if CSSM_PERSONALIZED in args.cssm:
            if args.cssm.split(";")[0]==CSSM_PERSONALIZED:
                if int(args.cssm.split(";")[1])>0:
                    if(int(args.cssm.split(";")[2])>0):
                        if(int(args.cssm.split(";")[3])>0):
                            logging.info("CSSM has the correct format")            
        else:
            logging.info("--unknown CSSM")
            logging.info("--Options are:")
            logging.info("----" + CSSM_PITA)
            logging.info("----" + CSSM_TS)
            logging.info("----" + CSSM_REGULAR)
            logging.info("----" + CSSM_PERSONALIZED)
            exit()
    else:
            logging.info("--missing CSSM")
            logging.error("----you must specify a CSSM using the -c/--cssm parameter")
            printHelpAndExit()    
    logging.info("--OK")
        
    
    
def checkDLModel():
    logging.info("checking DLModel:")
    if args.dlModel:
        logging.info(args.dlModel)
        if os.path.isfile(args.dlModel):
            logging.info(args.dlModel)
            logging.info("--DLModel file exists")
        else:
            logging.error("--can't find DLModel file")
            printHelpAndExit()
    else:
            logging.info("--missing DLModel")
            logging.error("----you must specify a DLModel using the -d/--dlModel parameter")
            printHelpAndExit()    
    logging.info("--OK")

    
    
def checkSplit():
    logging.info("checking split:")
    if args.splitType:
        logging.info(args.splitType)
        splitType=args.splitType
        
        if(int(args.splitType) == 0):
            logging.info("no split, keep as unified file as single file")
            return()
        
        if(int(args.splitType) == 1):
            logging.info("write unified file by miRNA")
            return()
        
        if(int(args.splitType) == 2):
            logging.info("write unified file by 3'UTR")
            return()
        
        if(int(args.splitType) == 3):
            logging.info("split unified file by miRNA")
            return()
        
        if(int(args.splitType) == 4):
            logging.info("split unified file by 3'UTR")
            return()
        
            
        logging.info("--unknown option")
        logging.info("--split can currently only take values 0 to 4")
        printHelpAndExit()    
    
                
    logging.info("--OK")


def checkEnergyFiltering():
    logging.info("EnergyFiltering")
    global filterByAccessibilityEnergy
    if args.energy_filtering:
        if args.energy_filtering==True:
            filterByAccessibilityEnergy=True
            logging.info("--energyFiltering is ON")
        else:
            logging.info("--energyFiltering is OFF")
    logging.info("--OK")
            

def checkWindowSize():
    logging.info("checking MaxSiteLength:")
    global maxSiteLength
    if args.maximumSiteLength:
        logging.info(args.maximumSiteLength)
        if(int(args.maximumSiteLength) > 0):
            maxSiteLength=int(args.maximumSiteLength)
        else:
            logging.info("--MaxSiteLength needs to be a positive number <" + args.maximumSiteLength + ">")
            return()
        
    logging.info("--OK")



def checkStepSize():
    logging.info("checking seedAlignmentOffset:")
    global seedAlignmentOffset
    if args.seedAlignmentOffset:
        logging.info(args.seedAlignmentOffset)
        if(int(args.seedAlignmentOffset) > 0):
            seedAlignmentOffset=int(args.seedAlignmentOffset)
        else:
            logging.info("--seedAlignmentOffset needs to be a positive number <" + args.seedAlignmentOffset + ">")
            return()
        
    logging.info("--OK")



def checkArgs():
    
    if args.HelpMe :
        printLongHelpAndExit() 
        exit()
    checkProgramLocation()
    checkExptName()
    checkOutFolder()
    checkUnifiedFile()
    checkCSSM()
    checkDLModel()
    checkSplit()
    checkEnergyFiltering()
    checkWindowSize()
    checkStepSize()
    
    
    
def readMiRNAFile():

    inFile = open(args.miRFile,'rU')
    for record in SeqIO.parse(inFile,'fasta'):
        miRheaderList.append(record.id)
        miRseqList.append(str(record.seq))
       
       
def readUTRFile():
    inFile = open(args.utrFile,'rU')
    for record in SeqIO.parse(inFile,'fasta'):
        UTRheaderList.append(record.id)
        UTRseqList.append(str(record.seq))
    
     
    
# Single Unified file: miRNAs are in the outer loop    
def writeUnifiedFileBymiRNA():
    # 
    logging.info("writeUnifiedFileBymiRNA")
    headerLine = "miRNA\tgene_name\tEnsemblId\tPositive_Negative\tMature_mirna_transcript\t3UTR_transcript" + MY_NEWLINE
    with open(unifiedFilePath, "w") as f:
        f.write(headerLine)
        m = 0
        for mHeader in miRheaderList:            
            u = 0
            for uHeader in UTRheaderList:               
                f.write('\t'.join([mHeader, uHeader, uHeader, "?", miRseqList[m], UTRseqList[u]]) + MY_NEWLINE)
                u+=1
            m+=1
 
            
    logging.info("--write Script File")
    #java -jar miRAW.jar GenePrediction predict ./Results/firstTry/firstTry.EF.properties        
    with open(os.path.join(args.outFolder, args.exptName  + '.sh'), "w") as f:
        f.write("java -jar " + jarLocation + " GenePrediction predict " 
                + os.path.join(args.outFolder, args.exptName  + ".properties") + MY_NEWLINE)    
  
   
    
    logging.info("done")
        
        
# Single Unified File: UTRs are in the outer loop    
def writeUnifiedFileByUTR():
    # 
    logging.info("writeUnifiedFileByUTR")
    headerLine = "miRNA\tgene_name\tEnsemblId\tPositive_Negative\tMature_mirna_transcript\t3UTR_transcript" + MY_NEWLINE
    with open(unifiedFilePath, "w") as f:
        f.write(headerLine)
        u = 0
        for uHeader in UTRheaderList:               
            m = 0
            for mHeader in miRheaderList:            
                f.write('\t'.join([mHeader, uHeader, uHeader, "?", miRseqList[m], UTRseqList[u]]) + MY_NEWLINE)
                m+=1
            u+=1


    logging.info("--write Script File")
    #java -jar miRAW.jar GenePrediction predict ./Results/firstTry/firstTry.EF.properties        
    with open(os.path.join(args.outFolder, args.exptName  + '.sh'), "w") as f:
        f.write("java -jar " + jarLocation + " GenePrediction predict " 
                + os.path.join(args.outFolder, args.exptName  + ".properties") + MY_NEWLINE)    
  
    logging.info("done")
    
    
# miRNAs are in the outer loop    
# we write one file per miRNA
# need to create a new .properties file for each miRNA as the unifiedFile name is stored in here
#
def splitUnifiedFileBymiRNA():
    # 
    logging.info("splitUnifiedFileBymiRNA")
    headerLine = "miRNA\tgene_name\tEnsemblId\tPositive_Negative\tMature_mirna_transcript\t3UTR_transcript" + MY_NEWLINE
    m = 0
    with open(os.path.join(args.outFolder, args.exptName + '.sh'), "w") as fSh:
        for mHeader in miRheaderList:       
            logging.info("--<" + mHeader + ">")
            #os.path.join(args.outFolder, args.exptName  + '.sh'), "w"
            featuresFile = os.path.join(args.outFolder, args.exptName  + "." + mHeader + '.unifiedFile.csv' )  
            propertiesFile = os.path.join(args.outFolder, args.exptName  + "." + mHeader + '.properties' )  
            with open(os.path.join(args.outFolder, featuresFile), 'w') as f:
                f.write(headerLine)
                u = 0
                for uHeader in UTRheaderList:               
                    f.write('\t'.join([mHeader, uHeader, uHeader, "?", miRseqList[m], UTRseqList[u]]) + MY_NEWLINE)
                    u+=1
    
            writePropertiesFileForFeature(mHeader)
            fSh.write("java -jar " + jarLocation + " GenePrediction predict " + propertiesFile + MY_NEWLINE)    
            m+=1
    #java -jar miRAW.jar GenePrediction predict ./Results/firstTry/firstTry.EF.properties        

  
    logging.info("done")
    
    
# UTRs are in the outer loop    
def splitUnifiedFileByUTR():
    # 
    logging.info("splitUnifiedFileByUTR")
    headerLine = "miRNA\tgene_name\tEnsemblId\tPositive_Negative\tMature_mirna_transcript\t3UTR_transcript" + MY_NEWLINE
    u = 0
    with open(os.path.join(args.outFolder, args.exptName + '.sh'), "w") as fSh:
        logging.info("")
        for uHeader in UTRheaderList:               
            logging.info("--<" + uHeader + ">")
            featuresFile = os.path.join(args.outFolder, args.exptName  + "." + uHeader.replace("|", "_") + '.unifiedFile.csv' )  
            propertiesFile = os.path.join(args.outFolder, args.exptName  + "." + uHeader.replace("|", "_") + '.properties' )  
            with open(os.path.join(args.outFolder, featuresFile), 'w') as f:
                m = 0
                f.write(headerLine)
                for mHeader in miRheaderList:            
                    f.write('\t'.join([mHeader, uHeader, uHeader, "?", miRseqList[m], UTRseqList[u]]) + MY_NEWLINE)
                    m+=1

            writePropertiesFileForFeature(uHeader.replace("|", "_"))
    
            fSh.write("java -jar " + jarLocation + " GenePrediction predict " + propertiesFile + MY_NEWLINE)    
            u+=1

    
    logging.info("done")
    
    
# format of .properties file used during miRAW testing and evaluation
    
# ########################################
# ###########GenePrediction###############
# ########################################
# ExperimentName=CV-MIRAW-Per
# ExperimentFolder=./Results
# UnifiedFile=./miRAWFile/unified.csv
# DLModel=./miRAWFile/PredModels/newModel.bin
# CandidateSiteFinder.Type=Personalized;10;1;7 <=This can be personalized to anything we want.  
def writePropertiesFile():
    
    # if the UnifiedFile contains a complete file path, we leave it alone
    # if it is just a filename, we append a folder name
    global unifiedFilePath
    global maximumSiteLength
    global seedAlignmentOffset
    unifiedFilePath = ""
    if os.path.dirname(args.unifiedFile) != "":
        unifiedFilePath = args.unifiedFile
    else:
        unifiedFilePath = os.path.join(args.outFolder, args.unifiedFile)
        
    with open(os.path.join(args.outFolder, args.exptName  + '.properties'), "w") as f:
        f.write("########################################" + MY_NEWLINE)
        f.write("# miRAW target prediction              #" + MY_NEWLINE)
        f.write("#                                      #" + MY_NEWLINE)
        f.write("#   generated from rawWrap             #" + MY_NEWLINE)
        f.write("#  " + str(datetime.datetime.now()) + "          #" + MY_NEWLINE)
        f.write("#                                      #" + MY_NEWLINE)
        f.write("########################################" + MY_NEWLINE)
        f.write("ExperimentName="+ args.exptName + MY_NEWLINE)
        f.write("ExperimentFolder=" + args.outFolder + MY_NEWLINE)
        f.write("UnifiedFile=" + unifiedFilePath + MY_NEWLINE)
        f.write("DLModel=" + args.dlModel + MY_NEWLINE)
        f.write("CandidateSiteFinder.Type=" + args.cssm + MY_NEWLINE)
        f.write("MaxSiteLength=" + str(maximumSiteLength) + MY_NEWLINE)
        f.write("SeedAlignmentOffset=" + str(seedAlignmentOffset) + MY_NEWLINE)
        if(filterByAccessibilityEnergy):
            f.write(FILTER_ACCESSIBILITY_ENERGY_STRING + MY_NEWLINE)
   

def writePropertiesFileForFeature(featureName):
    
    unifiedFilePath = os.path.join(args.outFolder, args.exptName  + "." + featureName + ".unifiedFile.csv")

        
        
    with open(os.path.join(args.outFolder, args.exptName  + "." + featureName + '.properties'), "w") as f:
        f.write("########################################" + MY_NEWLINE)
        f.write("# miRAW target prediction              #" + MY_NEWLINE)
        f.write("#                                      #" + MY_NEWLINE)
        f.write("#   generated from rawWrap             #" + MY_NEWLINE)
        f.write("#  " + str(datetime.datetime.now()) + "          #" + MY_NEWLINE)
        f.write("#                                      #" + MY_NEWLINE)
        f.write("########################################" + MY_NEWLINE)
        f.write("ExperimentName="+ args.exptName + "." + featureName + MY_NEWLINE)
        f.write("ExperimentFolder=" + args.outFolder + MY_NEWLINE)
        f.write("UnifiedFile=" + unifiedFilePath + MY_NEWLINE)
        f.write("DLModel=" + args.dlModel + MY_NEWLINE)
        f.write("CandidateSiteFinder.Type=" + args.cssm + MY_NEWLINE)
        f.write("MaxSiteLength=" + str(maximumSiteLength) + MY_NEWLINE)
        f.write("SeedAlignmentOffset=" + str(seedAlignmentOffset) + MY_NEWLINE)
        if(filterByAccessibilityEnergy):
            f.write(FILTER_ACCESSIBILITY_ENERGY_STRING + MY_NEWLINE)
   
        
    

def processAndBuild():
    if int(args.splitType)==UNIFIEDFILE_EXISTS:
        useExistingUnifiedFile()     
    else:
        buildNewUnifiedFile()



def useExistingUnifiedFile():

    # build properties file
    # write script file
    writePropertiesFile()
    #java -jar miRAW.jar GenePrediction predict ./Results/firstTry/firstTry.EF.properties        
    with open(os.path.join(args.outFolder, args.exptName  + '.sh'), "w") as f:
        f.write("java -jar " + jarLocation + " GenePrediction predict " 
                + os.path.join(args.outFolder, args.exptName  + '.properties') + MY_NEWLINE)       
    
    
    
    
def buildNewUnifiedFile():    
    readMiRNAFile()
    readUTRFile()
    
    if int(args.splitType) == UNIFIEDFILE_EXISTS:
        writePropertiesFile()
        writeUnifiedFileBymiRNA()
        return()

    if int(args.splitType) == WRITE_BY_MIRNA:
        writePropertiesFile()
        writeUnifiedFileBymiRNA()
        return()
    
    if int(args.splitType) == WRITE_BY_UTR:
        writePropertiesFile()
        writeUnifiedFileByUTR()
        return()
    
    if int(args.splitType) == SPLIT_BY_MIRNA:
        splitUnifiedFileBymiRNA()
        return()
    
    if int(args.splitType) == SPLIT_BY_UTR:
        splitUnifiedFileByUTR()
        return()
    
    
checkArgs()    
processAndBuild()


