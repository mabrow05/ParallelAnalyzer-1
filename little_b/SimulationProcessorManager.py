#!/usr/bin/python

#######################################################################
# This script, as is, is meant to manage the process for producing sets of 
# of linearity twiddles, running Ce, Sn, and Bi simulations through 
# the processor with these twiddles, testing whether the set of parameters
# was good or bad (peaks fit in calibration envelope), saving the good sets
# for each source, compiling a list of the parameters which work for all sources
# into one file, and also running the beta events through the processor using
# this final set of "good" parameters.
####################################################################### 

import os
import sys
from optparse import OptionParser

from math import *
import numpy as np


### This is the default error envelope as given by Michael M. (2008,2010) and Michael B. (2011,2012)
envelope = {2008:[(0,5.0),(250,5.0),(500,500*0.013),(900,900*0.025),(1000,1000*0.025),(1200,1200*0.025)],
            2010:[(0,2.5),(200,200*0.0125),(500,500*0.0125),(1000,500*0.0125)],
            2011:[(0,0.025*130.3),(130.3,130.3*0.025),(368.4938,0.018*368.4938),(993.789,993.789*0.013),(1000,1000.*0.013)],
            2012:[(0,2.5),(200,200*0.0125),(500,500*0.0125),(1000,500*0.0125)] }

### make sure proper directories exist to store trees and stuff
os.system("mkdir -p linCurves/")
os.system("mkdir -p analyzed_files/")


def binary_search_bool(arr,val,low=0,high=None): #Put in to check if parameter sets are good for the different sources
    high = high if high is not None and high<len(arr) else len(arr)-1

    if high==-1:
        return False

    comp = (high)/2
    if arr[comp]==val:
        return True 
    elif high>0:
        if val<arr[comp]:
            return binary_search_bool(arr[:comp],val)
        else:
            return binary_search_bool(arr[comp+1:],val)
    else: 
        return False



##### Make linearity param file 
def makeLinearityParamFile():
    paramDeltaRanges = {"p0":(-0.,0.), "p1":(-0.0,0.0), "p2":(-5.e-5,5.e-5), "p3":(-3.e-7,3.e-7)} #Range of parameter space one wants to cover for cubic linearity twiddle
    paramNSteps = {"p0":1, "p1":1, "p2":11, "p3":7} # number of steps over the parameter range

    paramFile = open("linCurves/parameters.dat",'w');
    paramSetIndex = 0

    #Gererates all combinations of parameters in the ranges above for each detector
    for Ep0 in np.linspace(paramDeltaRanges["p0"][0],paramDeltaRanges["p0"][1], paramNSteps["p0"]):
        for Ep1 in np.linspace(paramDeltaRanges["p1"][0],paramDeltaRanges["p1"][1], paramNSteps["p1"]):
            for Ep2 in np.linspace(paramDeltaRanges["p2"][0],paramDeltaRanges["p2"][1], paramNSteps["p2"]):
                for Ep3 in np.linspace(paramDeltaRanges["p3"][0],paramDeltaRanges["p3"][1], paramNSteps["p3"]):
                    for Wp0 in np.linspace(paramDeltaRanges["p0"][0],paramDeltaRanges["p0"][1], paramNSteps["p0"]):
                        for Wp1 in np.linspace(paramDeltaRanges["p1"][0],paramDeltaRanges["p1"][1], paramNSteps["p1"]):
                            for Wp2 in np.linspace(paramDeltaRanges["p2"][0],paramDeltaRanges["p2"][1], paramNSteps["p2"]):
                                for Wp3 in np.linspace(paramDeltaRanges["p3"][0],paramDeltaRanges["p3"][1], paramNSteps["p3"]):
                   
                                    paramFile.write("%i\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n"%(paramSetIndex,Ep0,Ep1,Ep2,Ep3,Wp0,Wp1,Wp2,Wp3))
                                    paramSetIndex+=1


# This will run all the source simulations through the Simulation Processor. If linCorr is true, it applies the twiddles in parameters.dat 
# and creates a file with the passing parameters for that source called passingParams_{geometry}_{src}.dat. Else, it runs the 
# simulation processor for that source with no corrections as a baseline
def runAllSourceSims(geometry="2010", numEvents=5000,linCorr=True):
    srcs = ["Bi207","Sn113", "Ce139"] 
    for src in srcs:
        if linCorr:
            os.system("rm linCurves/passingParams_%s_%s.dat"%(geometry[0:4],src))
            infile = open("linCurves/parameters.dat",'r')

            for line in infile:
                params = line.split()
                os.system("./SimulationProcessor %s %s %i %i %s %s %s %s %s %s %s %s %s"%(src,geometry,numEvents,linCorr,params[1],params[2],params[3],params[4],params[5],params[6],params[7],params[8],params[0]))
    
            infile.close()
        else:
            os.system("./SimulationProcessor %s %s %i %i"%(src,geometry,numEvents,linCorr))


##### NEEDS TO BE TESTED ############
# This will check the existing passingParams files for each source for common parameter sets and create a new file
# with these sets listed. It will not overwrite previous "good" sets, but it's the users job to keep them straight
# as they are simply numbered sequentially
def makeFinalParamsFile(geometry):
    CeFile = open("linCurves/passingParams_%s_Ce139.dat"%geometry[:4], 'r')
    SnFile = open("linCurves/passingParams_%s_Sn113.dat"%geometry[:4], 'r')
    BiFile = open("linCurves/passingParams_%s_Bi207.dat"%geometry[:4], 'r')

    Ce = []
    Sn = []
    Bi = []
    BiParamsList = []
    goodParams = []

    for line in CeFile:
        CeParams = line.split() 
        Ce.append(int(CeParams[0]))
    for line in SnFile:
        SnParams = line.split() 
        Sn.append(int(SnParams[0]))
    for line in BiFile:
        BiParams = line.split()
        BiParamsList.append(line)
        Bi.append(int(BiParams[0]))
        
    CeFile.close()
    SnFile.close()
    BiFile.close()

    num=0
    fileName = "linCurves/matchingParams_%s_%i.dat"%(geometry[:4],num)
    if (os.path.isfile(fileName)):
        num+=1
        fileName = "linCurves/matchingParams_%s_%i.dat"%(geometry[:4],num)

    finalFile = open(fileName,'w')

    for i,el in enumerate(Bi):
        print i,el 
        if binary_search_bool(Sn,el):
            if binary_search_bool(Ce,el):
                print BiParamsList[i]
                finalFile.write("%s\n"%BiParamsList[i])

    print "DONE writing final parameter file for %s"%geometry

# Function to run the beta decay simulations (both normal betas and fierz betas) through the processor.
# if linCorr is true, they will be run for every set of "good" parameters in the paramSet given. Else
# they will just be run with no twiddles as a baseline.
def runBetaSims(geometry="2010", numEvents=5000,linCorr=True,paramSet=0):
    srcs = ["Beta","Beta_fierz"]
    for src in srcs:
        if linCorr:
            infile = open("linCurves/matchingParams_%s_%i.dat"%(geometry[:4],paramSet),'r')
        
            for line in infile:
                params = line.split()
                os.system("./SimulationProcessor %s %s %i %i %s %s %s %s %s %s %s %s"%(src,geometry,numEvents,linCorr,params[1],params[2],params[3],params[4],params[5],params[6],params[7],params[8]))
    
            infile.close()
        else:
            os.system("./SimulationProcessor %s %s %i %i"%(src,geometry,numEvents,linCorr))






if __name__ == "__main__":

    parser = OptionParser()
    parser.add_option("--makeLinCorrFile",dest="makeLinCorrFile",action="store_true",default=False,
                      help="Makes file containing set of linearity corrections given the settings passed in main")
    parser.add_option("--runSourcesBaseline",dest="runSourcesBaseline",action="store_true",default=False,
                      help="Run each source with zero correction to get a baseline for the peaks")  
    parser.add_option("--runSourcesLinCorr",dest="runSourcesLinCorr",action="store_true",default=False,
                      help="Run sources through analyzer applying corrections from parameters.dat")
    parser.add_option("--makeGoodParamFile",dest="makeGoodParamFile",action="store_true",default=False,
                      help="Crosschecks linearity twiddles which are good for all source peaks and creates file")
    parser.add_option("--runBetaBaseline",dest="runBetaBaseline",action="store_true",default=False,
                      help="Runs the betas (both little b on and off) without Linearity twiddles")
    parser.add_option("--runBetaLinCorr",dest="runBetaLinCorr",action="store_true",default=False,
                      help="Runs the betas for all sets of parameters deemed \"good\"")
    parser.add_option("--nEvents","-n",type="int", action="store",dest="nEvents",default=5000,
                      help="Set the number of events to be run through the processor");
    parser.add_option("--geometry","-g", action="store",dest="geometry",default="2010",
                      help="Set the geometry of the simulation to be used in the processor");
    parser.add_option("--betaParamFile","-p",type="int", action="store",dest="paramSet",default=0,
                      help="Set the parameter set to be used for processing the betas");



    options, args = parser.parse_args()


    if options.makeLinCorrFile:
        makeLinearityParamFile()

    if options.runSourcesBaseline:
        runAllSourceSims(geometry=options.geometry,numEvents=options.nEvents,linCorr=False)

    if options.runSourcesLinCorr:
        runAllSourceSims(geometry=options.geometry,numEvents=options.nEvents,linCorr=True)

    if options.makeGoodParamFile:
        makeFinalParamsFile(geometry=options.geometry)

    if options.runBetaBaseline:
        runBetaSims(geometry=options.geometry,numEvents=options.nEvents,linCorr=False)

    if options.runBetaLinCorr:
        runBetaSims(geometry=options.geometry,numEvents=options.nEvents,linCorr=True,paramSet=options.paramSet)


    
 
    
