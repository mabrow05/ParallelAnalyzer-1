#!/usr/bin/python

import os
import sys
from optparse import OptionParser
from math import *
import numpy as np
#import ../scripts/MButils

### This is the default error envelope as given by Michael M. (2008,2010) and Michael B. (2011,2012)
envelope = {2008:[(0,5.0),(250,5.0),(500,500*0.013),(900,900*0.025),(1000,1000*0.025),(1200,1200*0.025)],
            2010:[(0,2.5),(200,200*0.0125),(500,500*0.0125),(1000,500*0.0125)],
            2011:[(0,0.025*130.3),(130.3,130.3*0.025),(368.4938,0.018*368.4938),(993.789,993.789*0.013),(1000,1000.*0.013)],
            2012:[(0,2.5),(200,200*0.0125),(500,500*0.0125),(1000,500*0.0125)] }

### make sure proper directories exits directory exists to store trees and stuff
os.system("mkdir -p linCurves/")
os.system("mkdir -p analyzed_files/")


def binary_search_bool(arr,val,low=0,high=None): #Put in to check if parameter sets are good for the different sources
    high = high if high is not None else len(arr)
    comp = high/2
    
    if arr[comp]==val
        return True 
    else:
        if val<arr[comp]:
            binary_search_bool(arr[:(comp/2),val,
        comp = comp/2 if val<arr[comp] else (len(arr)+comp)/2
        binary_search_bool(arr[

            

##### Make linearity param file 
def makeLinearityParamFile():
    paramDeltaRanges = {"p0":(-0.,0.), "p1":(-0.0,0.0), "p2":(-5.e-5,5.e-5), "p3":(-3.e-7,3.e-7)} 
    paramNSteps = {"p0":1, "p1":1, "p2":11, "p3":7}
    
    paramFile = open("linCurves/parameters.dat",'w');
    paramSetIndex = 0
    for p0 in np.linspace(paramDeltaRanges["p0"][0],paramDeltaRanges["p0"][1], paramNSteps["p0"]):
        for p1 in np.linspace(paramDeltaRanges["p1"][0],paramDeltaRanges["p1"][1], paramNSteps["p1"]):
            for p2 in np.linspace(paramDeltaRanges["p2"][0],paramDeltaRanges["p2"][1], paramNSteps["p2"]):
                for p3 in np.linspace(paramDeltaRanges["p3"][0],paramDeltaRanges["p3"][1], paramNSteps["p3"]):
                   
                    paramFile.write("%i\t%e\t%e\t%e\t%e\n"%(paramSetIndex,p0,p1,p2,p3))
                    paramSetIndex+=1


def runAllSourceSims(geometry="2010", numEvents=5000,linCorr=True):
    srcs = ["Sn113", "Ce139", "Bi207"]
    for src in srcs:
        if linCorr:
            os.system("rm linCurves/passingParams_%s_%s.dat"%(geometry,src))
            infile = open("linCurves/parameters.dat",'r')

            for line in infile:
                params = line.split()
                os.system("./SimulationAnalyzer %s %s %i %b %s %s %s %s %s"%(src,geometry,numEvents,linCorr,params[1],params[2],params[3],params[4],params[0],))
    
            infile.close()
        else:
            os.system("./SimulationAnalyzer %s %s %i %b"%(src,geometry,numEvents,linCorr))


def makeFinalParamsFile(geometry):
    CeFile = open("linCurves/passingParams_%s_Ce139.dat"%geometry, 'r')
    SnFile = open("linCurves/passingParams_%s_Sn113.dat"%geometry, 'r')
    BiFile = open("linCurves/passingParams_%s_Bi207.dat"%geometry, 'r')

    for line1 in CeFile:
        CeParams = line1.split()
        for line2 in SnFile:
            SnParams = line2.split()
            for line3 in BiFile:
                SnParams = line2.split()

makeLinearityParamFile()
runAllSourceSims(geometry="2010")
exit(0)




makeLinearityParamFile()




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
    parser.add_option("--runBetaBase",dest="runBetaBase",action="store_true",default=False,
                      help="Runs the betas (both little b on and off) without Linearity twiddles")
    parser.add_option("--runBetaLinCorr",dest="runBetaLinCorr",action="store_true",default=False,
                      help="Runs the betas for all sets of parameters deemed \"good\"")
    #Option to tell whatever you are running how many events to pass
    #Option to set geometry
    #option to set which set of "good" parameters to use
    
    parser.add_option("--runGainBismuth",dest="runGainBismuth",action="store_true",default=False,
                      help="Calculate and plot the Bi pulser spectra and track the gain relative to the reference run")
    parser.add_option("--makeDirectories",dest="makeDirectories",action="store_true",default=False,
                      help="Makes all the analysis directories.")
                      

    options, args = parser.parse_args()


    if options.makeLinCorrFile:
        makeLinearityParamFile()


 
    if options.createOctetLists:
        beta = BetaReplayManager()
        beta.createOctetLists()
    
    if options.findPedestals:
        beta = BetaReplayManager()
        #for octet in range(0,60,1):
        beta.findPedestals(5)

    if 0: 
        asymm = BetaAsymmetryManager()
        asymm.makeOctetAnalysisDirectories()

    if 0:
        beta = BetaReplayManager()
        for octet in range(0,59,1):
            beta.makeBasicHistograms(octet)


    if 1:
        octet_range = [0,59]#[20,28]#[45,50]#[38,40]#[0,59];
        beta = BetaReplayManager()
        for octet in range(octet_range[0],octet_range[1]+1,1):
            #beta.findPedestals(octet)
            #beta.runReplayPass1(octet)
            #beta.runGainBismuth(octet)
            #beta.runReplayPass2(octet)
            #beta.runReplayPass3(octet)
            beta.runReplayPass4(octet)


    #Running reverse calibrations
    if 0:
        octet_range = [0,2];
        beta = BetaReplayManager()
        for octet in range(octet_range[0],octet_range[1]+1,1):
            beta.runReverseCalibration(octet)
    
