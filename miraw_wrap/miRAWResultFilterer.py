import argparse
import sys
import os
import logging
import csv
import datetime

MY_NEWLINE = "\n"
if os.name== "Windows":
    MY_NEWLINE ="\r\n"
    
EXTRACT_POSITIVE_PREDICTIONS    = "EXPOS"
EXTRACT_NEGATIVE_PREDICTIONS    = "EXNEG"
FILTER_POSITIVE_PREDICTIONS     = "FTPOS"
FILTER_NEGATIVE_PREDICTIONS     = "FTNEG"
FILTER_BY_MIRNA                 = "FTMIR"

TARGET_FILE_EXTENSION           = "csv"
FILTERED_TARGET_FILE_EXTENSION  = "filtered.csv"
SITES_FILE_EXTENSION            = "csv"
FILTERED_SITES_FILE_EXTENSION   = "filtered.csv"
UNIFIED_FILE_EXTENSION          = "csv"
FILTERED_UNIFIED_FILE_EXTENSION = "filtered.csv"   
FILTERED_DETAILS_FILE_EXTENSION = "detailed.csv" 
    
UTR_QUERY_DELIMITER             = ":::"    
bindingEnergyCutoff = -1000000.0
positiveProbabilityCutoff = 0.0
negativeProbabilityCutoff = 0.0
removePositivePredictions = False
extractNegativePredictions = False
extractConflictPredictions = False
filterPredictionsByList = False
removePositivePredictions = False
keepPositivePredictions = False


miRAWtargetPredictionFile = ""
miRAWpositiveSitesFile = ""
miRAWunifiedFile = ""
mirnaFilterFile = ""
filteredTargetPredictionFile = ""
filteredPositiveSitesFile = ""
filteredUnifiedFile = ""
filteredPositiveSitesDetailedFile = ""


filterDropCount = 0
conflictDropCount = 0
keepDropCount = 0
posDropCount = 0
negDropCount = 0
posProbDropCount = 0
negProbDropCount = 0
emptyCount = 0
totalDrops = 0

headerString = "#" + MY_NEWLINE


#GeneName    GeneId    miRNA    Prediction    HighestPredVal    LowestPredVal    PosSites    NegSites    RemovedSites        miRNA    miRNA
GENE_NAME_COL          = 0
GENE_ID_COL            = 1
MIRNA_NAME_COL         = 2
PREDICTION_COL         = 3
HIGHEST_PREDICTION_COL = 4
LOWEST_PREDICTION_COL  = 5
POSITIVE_SITES_COL     = 6
NEGATIVE_SITES_COL     = 7
REMOVED_SITES_COL      = 8

GENE_NAME_ID           = "GeneName"
GENE_ID_ID             = "GeneID"
MIRNA_NAME_ID          = "miRNA"
PREDICTION_ID          = "Prediction"
HIGHEST_PREDICTION_ID  = "HighestPredVal"
LOWEST_PREDICTION_ID   = "LowestPredVal"
POSITIVE_SITES_ID      = "PosSites"
NEGATIVE_SITES_ID      = "NegSites"
REMOVED_SITES_ID       = "RemovedSites"

UNIFIED_UTRID_COL      = 2
UNIFIED_UTRSEQ_COL     = 5



HEADER_LINE = GENE_NAME_ID + "\t"\
                + GENE_ID_ID + "\t"\
                + MIRNA_NAME_ID + "\t"\
                + PREDICTION_ID + "\t"\
                + HIGHEST_PREDICTION_ID + "\t"\
                + LOWEST_PREDICTION_ID + "\t"\
                + POSITIVE_SITES_ID + "\t"\
                + NEGATIVE_SITES_ID + "\t"\
                + REMOVED_SITES_ID
                


geneNames = []
geneIDs = []
miRNANames = []
predictions = []
highestPredVal = []
lowestPredVal = []
positiveSites = []
negativeSites = []
removedSites = []
dropPredictions = []          # use this to keep track of which predictions to drop
positiveSiteLines = []
positiveSiteLocations = {}

filteredMiRNAList = []
utrSequences = {}
miRNASequences = {}

logging.getLogger().setLevel(logging.INFO)

parser = argparse.ArgumentParser(description='set up batch commands to run miRAW')

parser.add_argument("-H", "--HelpMe", action="store_true",
                    help="print detailed help")

parser.add_argument("-t", "--target_pred_file", dest='targetPredFile',
                    help="miRAW summary of results file (.targetPredictionOutput.csv)")

parser.add_argument("-s", "--target_site_file", dest='targetSiteFile',
                    help="miRAW detailed results file (.positiveTargetSites.csv)")

parser.add_argument("-u", "--unified_input_file", dest='unifiedInputFile',
                    help="miRAW input file used for target prediction (.unifiedFile.csv)")

parser.add_argument("-e", "--energy-cutoff", dest='energyCutOff',
                    help="remove predictions ABOVE this energy (i.e. weaker bindings)")

parser.add_argument("-p", "--pos-prob-cutoff", dest='posProbCutOff',
                    help="remove positive predictions BELOW this prediction probability")

parser.add_argument("-n ", "--neg-prob-cutoff", dest='negProbCutOff',
                    help="remove negative predictions ABOVE this prediction probability")

parser.add_argument("-L ", "--extract_from_list", dest='extractListFile',
                    help="remove negative predictions ABOVE this prediction probability")

parser.add_argument("-R", "--remove_conflicts", action="store_true",
                    help="remove predictions with both positive and negative predictions")

parser.add_argument("-K", "--keep_conflicts", action="store_true",
                    help="keep predictions with both positive and negative predictions")

parser.add_argument("-P", "--extract_pos_preds", action="store_true",
                    help="keep positive predictions, remove negative prediction")

parser.add_argument("-N", "--extract_neg_preds", action="store_true",
                    help="keep negative predictions, remove positive predictions")

args = parser.parse_args()


def printLongHelpAndExit():
    logging.info("+" + "-" * 78 + "+")
    logging.info("+  miRAWResultFilterer                                                         +")
    logging.info("+    wrapper code for generating files and scripts to run miRAW                +")
    logging.info("+                                                                              +")
    logging.info("+    you need to specify:                                                      +")
    logging.info("+                                                                              +")
    logging.info("+      a miRAW summary result file as input        (-t/--target_pred_file)     +")
    logging.info("+        (this is the file which ends in 'targetPredictionOutput.csv')         +")
    logging.info("+                                                                              +")
    logging.info("+      a miRAW detailed results file as input      (-s/--target_site_file)     +")
    logging.info("+        (this is the file which ends in 'positiveTargetSites.csv')            +")
    logging.info("+                                                                              +")
    logging.info("+      the input file used for prediction          (-u/--unified_input_file)   +")
    logging.info("+        (this is the file which ends in 'unifiedFile.csv')                    +")
    logging.info("+                                                                              +")
    logging.info("+      remove conflicts                            (-R/--remove_conflicts)     +")
    logging.info("+         (remove predictions that contain both                                +")
    logging.info("+          positive and negative predictions                                   +")
    logging.info("+          This isn't necessarily a problem, it's                              +")
    logging.info("+          just another way to filter the data)                                +")
    logging.info("+                                                                              +")
    logging.info("+      keep conflicts                              (-K/--keep_conflicts)       +")
    logging.info("+         (keep predictions that contain both                                  +")
    logging.info("+          positive and negative predictions)                                  +")
    logging.info("+                                                                              +")
    logging.info("+      binding energy cutoff:                      (-e/--energy-cutoff)        +")
    logging.info("+        (positive predictions with a probability below                        +")
    logging.info("+         this value will be removed                                           +")
    logging.info("+                                                                              +")
    logging.info("+      positive probability cutoff:                (-p/--pos_prob_cutoff)      +")
    logging.info("+        (positive predictions with a probability below                        +")
    logging.info("+         this value will be removed                                           +")
    logging.info("+         values can be between 0 --> 1                                        +")
    logging.info("+         higher values are better)                                            +")
    logging.info("+                                                                              +")
    logging.info("+      negative probability cutoff:                (-n/--neg_prob_cutoff)      +")
    logging.info("+        (negative predictions with a probability above                        +")
    logging.info("+         this value will be removed                                           +")
    logging.info("+         values can be between 0 --> -1                                       +")
    logging.info("+         lower values are better)                                             +")
    logging.info("+                                                                              +")
    logging.info("+      keep all positive predictions               (-P/--extract_positives)    +")
    logging.info("+      remove all negative predictions                                         +")
    logging.info("+                                                                              +")
    logging.info("+      keep all negative predictions               (-N/--extract_negatives)    +")
    logging.info("+      remove all positive predictions                                         +")
    logging.info("+                                                                              +")
    logging.info("+      extract miRNAs in the supplied list.        (-L/--extract_mirnas)       +")
    logging.info("+          remove miRNAs that are listed in the                                +")
    logging.info("+          supplied file                                                       +")
    logging.info("+          (use this to extract miRNAs you don't want                          +")
    logging.info("+          for example, high numbered/low confidence                           +")
    logging.info("+          or miRNAs in which you are particularly interested)                 +")
    logging.info("+                                                                              +")
    logging.info("+   filtering is performed in the order                                        +")
    logging.info("+      listFilter->Conflicts->PosCutOff->NegCutOff->ExtractPos/Neg Preds       +")
    logging.info("+                ->BindingEnergy                                               +")
    logging.info("+                                                                              +")
    logging.info("+" + "-" * 78 + "+")


def printHelpAndExit():
    parser.print_help()
    logging.info("stopping") 
    sys.exit()

def checkExtractConflicts(): 
    logging.info("extractConflicts")       
    global extractConflictPredictions, keepConflictPredictions
    if args.remove_conflicts & args.keep_conflicts: 
        logging.error("you can't set both '--keep_conflicts' and '--remove_conflicts'")
        logging.info("stopping") 
        sys.exit()
    if args.remove_conflicts:
        if args.remove_conflicts==True:
            extractConflictPredictions=True
            keepConflictPredictions = False
            logging.info("--remove_conflicts is True")
        else:
            logging.info("--remove_conflicts is False")
    if args.keep_conflicts:
        if args.keep_conflicts==True:
            extractConflictPredictions=False
            keepConflictPredictions=True
            logging.info("--keep_conflicts is True")
        else:
            logging.info("--keep_conflicts is False")
    logging.info("--OK")
    
    
def checkBindingEnergyCutOff():
    logging.info("checkBindingEnergyCutOff") 
    
    global bindingEnergyCutoff  
    if args.energyCutOff:
        logging.info(args.energyCutOff)
        if float(args.energyCutOff) >= 0.0:
            bindingEnergyCutoff=float(args.energyCutOff)
        else:
            logging.info(" binding energy cutoff (-e) needs to be less than 0" + str(args.energyCutOff) + ">")
            return()
     

def checkPositiveCutOff():
    logging.info("checkPositiveCutOff") 
    
    global positiveProbabilityCutoff  
    if args.posProbCutOff:
        logging.info(args.posProbCutOff)
        if( float(args.posProbCutOff) >= 0.0) & (float(args.posProbCutOff) <= 1.0):
            positiveProbabilityCutoff=float(args.posProbCutOff)
        else:
            logging.info(" positive probability cutoff (-p) needs to be between 0 and 1" + str(args.posProbCutOff) + ">")
            return()
     

def checkNegativeCutOff():
    logging.info("checknegativeCutOff")   
    
    global negativeProbabilityCutoff
    if args.negProbCutOff:
        logging.info(args.negProbCutOff)
        if (float(args.negProbCutOff) <= 0.0) & (float(args.negProbCutOff) >= -1.0):
            negativeProbabilityCutoff=float(args.negProbCutOff)
        else:
            logging.info(" negative probability cutoff (-n) needs to be between -1 and 0" + str(args.negProbCutOff) + ">")
            return()
    


def checkExtractPositives():
    logging.info("checkExtractPositives")
    global removePositivePredictions
    if args.extract_pos_preds:
        if args.extract_pos_preds==True:
            removePositivePredictions=True
            logging.info("--extract_pos_preds is True")
        else:
            logging.info("--extract_pos_preds is False")
    
    logging.info("--OK")


def checkExtractNegatives():
    logging.info("checkExtractNegatives")
    global extractNegativePredictions
    if args.extract_neg_preds:
        if args.extract_neg_preds==True:
            extractNegativePredictions=True
            logging.info("--extract_neg_preds is True")
        else:
            logging.info("--extract_neg_preds is False")
    logging.info("--OK")
    
    
def checkRemovemiRNAs():
    logging.info("check removeMiRNAsInList File")
    global mirnaFilterFile, filterPredictionsByList
    if args.extractListFile:
        filterPredictionsByList = True
        mirnaFilterFile=args.extractListFile
        logging.info("--extractListFile is " + mirnaFilterFile)
        if not os.path.isfile(mirnaFilterFile):
            logging.error("--can't find miRNA Extract List File at <" + mirnaFilterFile + ">")
            exit()
    logging.info("--OK")
        

def checkTargetPredictionFile():
    global miRAWtargetPredictionFile, filteredTargetPredictionFile            
    logging.info("check TargetPredictionFile")
    if args.targetPredFile:
        miRAWtargetPredictionFile=args.targetPredFile
        filteredTargetPredictionFile=miRAWtargetPredictionFile.replace(TARGET_FILE_EXTENSION, \
                                                                       FILTERED_TARGET_FILE_EXTENSION)
        logging.info("--Target Prediction File is " + miRAWtargetPredictionFile)
        logging.info("--Filtered Target Prediction File is " + filteredTargetPredictionFile)
        if not os.path.isfile(miRAWtargetPredictionFile):
            logging.error("--can't find Target Prediction File at <" + miRAWtargetPredictionFile + ">")
            exit()
    logging.info("--OK")


def checkPositiveTargetsFile():
    global miRAWpositiveSitesFile, filteredPositiveSitesFile
    global filteredPositiveSitesDetailedFile
    logging.info("check PositiveTargetsFile")
    if args.targetSiteFile:
        miRAWpositiveSitesFile=args.targetSiteFile
        filteredPositiveSitesFile=miRAWpositiveSitesFile.replace(SITES_FILE_EXTENSION, \
                                                                       FILTERED_SITES_FILE_EXTENSION)
        filteredPositiveSitesDetailedFile=miRAWpositiveSitesFile.replace(SITES_FILE_EXTENSION, \
                                                                       FILTERED_DETAILS_FILE_EXTENSION)
        
        logging.info("--Positive Sites File is " + miRAWpositiveSitesFile)
        logging.info("--Filtered Positive Sites File is " + filteredPositiveSitesFile)
        logging.info("--Detailed Positive Sites File is " + filteredPositiveSitesDetailedFile)
        if not os.path.isfile(miRAWpositiveSitesFile):
            logging.error("--can't find Positive Sites File at <" + miRAWpositiveSitesFile + ">")
            exit()
    logging.info("--OK")


def checkUnifiedInputFile():
    global miRAWunifiedFile, filteredUnifiedFile
    logging.info("check unifiedFile")
    if args.unifiedInputFile:
        miRAWunifiedFile=args.unifiedInputFile
        filteredUnifiedFile=miRAWunifiedFile.replace(UNIFIED_FILE_EXTENSION, \
                                                                       FILTERED_UNIFIED_FILE_EXTENSION)
        logging.info("--Unified Input File is " + miRAWunifiedFile)
        logging.info("--Unified Input File is " + filteredUnifiedFile)
        if not os.path.isfile(miRAWunifiedFile):
            logging.error("--can't find Unified Input File at <" + miRAWunifiedFile + ">")
            exit()
    logging.info("--OK")


def checkForBadFlagCombinations():
    if(removePositivePredictions & extractNegativePredictions):
        logging.warn("you have both --extract_pos_preds and --extract_neg_preds flags set")
        logging.warn("nothing will be filtered with this combination")
        


def readPositiveSitesFile():
    logging.info("read PositiveSiteFile")
    global miRAWpositiveSitesFile, positiveSiteLines
    utrNames = []

    with open(miRAWpositiveSitesFile, 'r') as fP:
        positiveSiteLines = fP.readlines()
        p=0
        for line in positiveSiteLines:
            if p==0:
                p+=1
                continue
            if not line.split("\t")[0].isdigit():
                positiveSiteLocations[line.split("\t")[0].strip() + UTR_QUERY_DELIMITER + line.split("\t")[1].strip()]=p
                miRNASequences[line.split("\t")[1].strip()] = line.split("\t")[2].strip()
                utrNames.append(line.split("\t")[0].strip() + UTR_QUERY_DELIMITER + line.split("\t")[1].strip())
            p+=1    
                
    logging.info("--read <" + str(len(positiveSiteLines)) + "> lines")            
    logging.info("--read <" + str(len(positiveSiteLocations)) + "> target pairs")            
    logging.info("--done")
        
 
    
    
    
def readUnifiedFile():
    logging.info("read UnifiedFile")
    global utrSequences
    
    with open(miRAWunifiedFile, 'r') as fF:
        utrLines = fF.readlines()
        u=0
        for line in utrLines:
            if u!=0:
                utrID = line.split("\t")[UNIFIED_UTRID_COL]
                if not utrID in utrSequences:
                    utrSequences[utrID] = line.split("\t")[UNIFIED_UTRSEQ_COL]
            u+=1       
    logging.info("--read " + str(len(utrSequences)) + " 3'UTR sequence(s)")   
    logging.info("--done")     

    


    
    
def readTargetPredictionFile():
    logging.info("read TargetPredictionFile")
    global miRAWtargetPredictionFile      

    with open(miRAWtargetPredictionFile,'r') as fT:
 
        targetLines = fT.readlines()
        l=0
        for line in targetLines:
            #logging.info(line)
            if l!=0:
                geneNames.append(line.split("\t")[GENE_NAME_COL])
                #reader=csv.reader(fT, delimiter='\t')
                #for geneName, geneID, mirName, pred, highPred, lowPred, posSite, negSite, remSite in reader:
                geneIDs.append(line.split("\t")[GENE_ID_COL])
                miRNANames.append(line.split("\t")[MIRNA_NAME_COL])
                predictions.append(line.split("\t")[PREDICTION_COL])
                highestPredVal.append(line.split("\t")[HIGHEST_PREDICTION_COL])
                lowestPredVal.append(line.split("\t")[LOWEST_PREDICTION_COL])    
                positiveSites.append(line.split("\t")[POSITIVE_SITES_COL])
                negativeSites.append(line.split("\t")[NEGATIVE_SITES_COL])
                removedSites.append(line.split("\t")[REMOVED_SITES_COL])     
                dropPredictions.append(0)
            l+=1
    logging.info("--read " + str(len(geneNames)) + " predictions")   
    logging.info("--done")     
            
            
            
def readFilterFile():
    logging.info("readFilterFile")
    global mirnaFilterFile, filteredMiRNAList
    with open(mirnaFilterFile, 'r') as fM:        
        filteredMiRNAList = fM.readlines()
    logging.info("--read " + str(len(filteredMiRNAList)) + " entries")
    
    
                
def processTargetResults():
    logging.info("process target results")



def filterByList():
    logging.info("filter by list")
    global filterDropCount, filteredMiRNAList, totalDrops
    filterDropCount = 0
    for filterMiRName in filteredMiRNAList:
        filterMiRName = filterMiRName.strip()
        r = 0
        for miRNAName in miRNANames:
            logging.debug(filterMiRName + "\t" + miRNAName)
            if(filterMiRName == miRNAName):
                dropPredictions[r] = 1
                filterDropCount += 1
            r += 1

    logging.info("--Dropped " + str(filterDropCount) + " entries")
    totalDrops += filterDropCount
    
    

def filterByProbability():
    logging.info("filter by prediction probability")
    global posProbDropCount, negProbDropCount, totalDrops
    global positiveProbabilityCutoff, negativeProbabilityCutoff
    
    r = 0
    posProbDropCount = 0
    negProbDropCount = 0
    while r<len(dropPredictions):        
        #logging.info(str(positiveSites[r]) + "|"  \
        #             + str(negativeSites[r]) + "\t" + str(positiveProbabilityCutoff) \
        #             + "|"+ str(negativeProbabilityCutoff) + "\t" + str(predictions[r]) )
        if int(positiveSites[r]) > 0:
            if positiveProbabilityCutoff > float(predictions[r]): 
                dropPredictions[r]=1
                posProbDropCount += 1
        else:
            if int(negativeSites[r]) > 0:
                if negativeProbabilityCutoff < float(predictions[r]):
                    dropPredictions[r]=1
                    negProbDropCount += 1
        r += 1
        
    logging.info("--Dropped " + str(posProbDropCount) + " positive entries")
    logging.info("--Dropped " + str(negProbDropCount) + " negative entries")
    totalDrops += posProbDropCount + negProbDropCount

    

                                
def filterByPositive():
    logging.info("filter by positives")
    global posDropCount, totalDrops
    
    r = 0
    posDropCount = 0
    while r<len(dropPredictions):
        if int(positiveSites[r]) == 0 :
            dropPredictions[r]=1
            posDropCount += 1
        r += 1
        
    logging.info("--Dropped " + str(posDropCount) + " entries")
    totalDrops += posDropCount


                    
def filterByNegative():
    logging.info("filter by negatives")
    global negDropCount, totalDrops
    
    r = 0
    negDropCount = 0
    while r<len(dropPredictions):
        if int(negativeSites[r]) == 0 :
            dropPredictions[r]=1
            negDropCount += 1
        r += 1
        
    logging.info("--Dropped " + str(negDropCount) + " entries")
    totalDrops += negDropCount


def removeEmptyPredictions():
    global emptyCount, totalDrops
    
    r=1
    while r<len(dropPredictions):        

        if int(positiveSites[r]) == 0 & int(positiveSites[r]) == 0:
            dropPredictions[r]=1
            emptyCount += 1
        r += 1
        
    logging.info("--Dropped " + str(emptyCount) + " empty entries")
    totalDrops += emptyCount
    

    
def keepByConflicts():
    global keepDropCount, totalDrops
    
    r = 0
    keepDropCount = 0
    while r<len(dropPredictions):
        if not (int(positiveSites[r]) != 0) & (int(negativeSites[r]) != 0):
            dropPredictions[r]=1
            keepDropCount += 1
        r += 1
        
    logging.info("--Kept " + str(keepDropCount) + " entries")
    totalDrops += keepDropCount


def filterByConflicts():
    global loggingconflictDropCount, totalDrops
    
    r = 0
    conflictDropCount = 0
    while r<len(dropPredictions):
        if (int(positiveSites[r]) != 0) & (int(negativeSites[r]) != 0):
            dropPredictions[r]=1
            conflictDropCount += 1
        r += 1
        
    logging.info("--Dropped " + str(conflictDropCount) + " entries")
    totalDrops += conflictDropCount


def summarizeFiltering():
    global totalDrops
    totalDrops = 0
    for drops in dropPredictions:
        totalDrops += int(drops)
        
    logging.info("--Dropped a total of " + str(totalDrops) + " entries from " \
                 + str(len(dropPredictions)) + " predictions")





def grabTargetPairLines(queryLine): 
    global positiveSiteLocations, positiveSiteLines
    logging.info(queryLine)
    logging.info(positiveSiteLocations[queryLine])
    
    lineNo = positiveSiteLocations[queryLine]
    queryLines = []
    l = lineNo + 1
    try:
        while(positiveSiteLines[l].split("\t")[1].isdigit()):
            queryLines.append(positiveSiteLines[l])
            l+=1
        return queryLines
    
    except IndexError: # in case of empty line at end of file
        return queryLines
    
    
     



def writeFilteredDetailedData():
    global filteredPositiveSitesDetailedFile
    logging.info(" write filtered detailed Data")
    
    r = 0
    with open(filteredPositiveSitesDetailedFile,'w') as fD:
        while r<len(dropPredictions):
            if dropPredictions[r]==0:
                targetQuery=geneIDs[r]+ UTR_QUERY_DELIMITER + miRNANames[r]
                
                predictionOutputString = geneNames[r] + "\t" + \
                          geneIDs[r] + "\t" + \
                          miRNANames[r] + "\t" + \
                          predictions[r] + "\t" + \
                          highestPredVal[r] + "\t" + \
                          lowestPredVal[r] + "\t" + \
                          positiveSites[r] + "\t" + \
                          negativeSites[r] + "\t" + \
                          removedSites[r] + "\t" 
                # Need to collect all the other data related to this entry
                # from the positive sites file
                # and the sequence from the unifiedFile
                
                targetLines = grabTargetPairLines(targetQuery)
                for targetLine in targetLines:
                    utrStart = int(targetLine.split("\t")[0])-1
                    utrStop = int(targetLine.split("\t")[1])-1
                    logging.info(str(utrStart) + ":" + str(utrStop))
                    fD.write(predictionOutputString + "\t" + \
                             targetLine.strip() + "\t" + \
                             miRNASequences[miRNANames[r]] + "\t" + \
                             utrSequences[geneIDs[r]][utrStart:utrStop] + \
                             MY_NEWLINE)

            r += 1
    logging.info("--finished")
    
    

def writeFilteredSummaryData():
    logging.info(" write filtered Data")
    
    global filteredTargetPredictionFile, headerString, totalDrops
    global conflictDropCount, keepDropCount, emptyCount, filterCount
    global negDropCount, posDropCount, negProbDropCount, posProbDropCount
    global negativeProbabilityCutoff, positiveProbabilityCutoff    
    headerString = headerString + "#             " + MY_NEWLINE
    headerString = headerString + "#   filtered target prediction data" + MY_NEWLINE
    headerString = headerString + "#   generated by miRAWResultFilterer" + MY_NEWLINE
    headerString = headerString + "#  " + str(datetime.datetime.now()) + MY_NEWLINE
    headerString = headerString + "#             " + MY_NEWLINE
    if keepDropCount>0:
        headerString = headerString + "#   kept " + str(keepDropCount) + \
                    " entries with both positive + negative predictions#" + MY_NEWLINE
    else:
        headerString = headerString + "#   removed " + str(conflictDropCount) + \
                    " entries with both positive + negative predictions#" + MY_NEWLINE
    headerString = headerString + "#   removed " + str(emptyCount) + \
                " entries with no predictions#" + MY_NEWLINE
    headerString = headerString + "#   removed " + str(negDropCount) + \
                " entries by dropping negative entries#" + MY_NEWLINE
    headerString = headerString + "#   dropped " + str(posDropCount) + \
                " entries by dropping positive entries#" + MY_NEWLINE
    headerString = headerString + "#   removed " + str(posProbDropCount) + \
                " positive entries with prediction probability < "  + str(negativeProbabilityCutoff) + MY_NEWLINE
    headerString = headerString + "#   removed " + str(negProbDropCount) + \
                " negative entries with prediction probability > " + str(positiveProbabilityCutoff) + MY_NEWLINE
    headerString = headerString + "#   removed " + str(filterDropCount) + \
                " entries with matches in miRNA name filter file <" + \
                mirnaFilterFile + ">" + MY_NEWLINE
    headerString = headerString + "#   total removed entries: " + str(totalDrops) +  \
                " out of " + str(len(geneIDs)) + " predictions"
    headerString = headerString + "#             " + MY_NEWLINE
    
    r = 0
    with open(filteredTargetPredictionFile,'w') as fT:
        fT.write(headerString + MY_NEWLINE)
        fT.write(HEADER_LINE+ MY_NEWLINE)
        while r<len(dropPredictions):
            if dropPredictions[r]==0:
                fT.write( geneNames[r] + "\t" + \
                          geneIDs[r] + "\t" + \
                          miRNANames[r] + "\t" + \
                          predictions[r] + "\t" + \
                          highestPredVal[r] + "\t" + \
                          lowestPredVal[r] + "\t" + \
                          positiveSites[r] + "\t" + \
                          negativeSites[r] + "\t" + \
                          removedSites[r] + MY_NEWLINE)

            r += 1
    logging.info("--finished")
            
            
            
    
def checkArgs():
    
    if args.HelpMe :
        printLongHelpAndExit() 
        exit()
        
    checkTargetPredictionFile()
    checkPositiveTargetsFile()
    checkUnifiedInputFile()
    checkRemovemiRNAs()
    checkExtractConflicts()
    checkBindingEnergyCutOff()
    checkPositiveCutOff()
    checkNegativeCutOff()
    checkExtractPositives()
    checkExtractNegatives()
    checkForBadFlagCombinations()


def filterData():
    global extractConflictPredictions, keepConflictPredictions
    logging.info("filterData")
    removeEmptyPredictions()
    if filterPredictionsByList:
        readFilterFile()
        filterByList()
    if extractConflictPredictions:
        filterByConflicts()
    if keepConflictPredictions:
        keepByConflicts()
    if (positiveProbabilityCutoff > 0) or (negativeProbabilityCutoff < 0):
        filterByProbability()
    if removePositivePredictions:
        filterByPositive()
    if extractNegativePredictions:
        filterByNegative()
    summarizeFiltering()

    


checkArgs()
readUnifiedFile()
readPositiveSitesFile()
readTargetPredictionFile()

filterData()
writeFilteredSummaryData()
writeFilteredDetailedData()


