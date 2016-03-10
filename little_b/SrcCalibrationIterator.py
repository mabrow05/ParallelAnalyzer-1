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


##### Make linearity curve files for order n polynomial
def makeLinearityParamFile():
    polOrder = 3 # cubic polynomial is default
    paramDeltas = {"p0":.0, "p1":0., "p2":0., "p3":0.}
    paramDeltaRanges = {"p0":(-0.,0.), "p1":(-0.2,0.2), "p2":(-0.0005,0.0005), "p3":(-0.000003,0.00003)}
    paramNSteps = {"p0":1, "p1":5, "p2":11, "p3":15}
    
    paramFile = open("linCurves/parameters.dat",'w');
    
    for p0 in np.linspace(paramDeltaRanges["p0"][0],paramDeltaRanges["p0"][1], paramNSteps["p0"]):
        for p1 in np.linspace(paramDeltaRanges["p1"][0],paramDeltaRanges["p1"][1], paramNSteps["p1"]):
            for p2 in np.linspace(paramDeltaRanges["p2"][0],paramDeltaRanges["p2"][1], paramNSteps["p2"]):
                for p3 in np.linspace(paramDeltaRanges["p3"][0],paramDeltaRanges["p3"][1], paramNSteps["p3"]):
                    paramDeltas["p0"] = p0;
                    paramDeltas["p1"] = p1;
                    paramDeltas["p2"] = p2;
                    paramDeltas["p3"] = p3;
                    paramFile.write("%f\t%f\t%f\t%f\n"%(p0,p1,p2,p3))
                    
    exit(0)

def runAllSourceSims(numEvents=100000):
    srcs = ["Sn113"]

    infile = open("linCurves/parameters.dat",'r')
    for line in infile:
        params = line.split()
        for src in srcs:
            os.system("./SimulationAnalyzer %s %i true %s %s %s %s"%(src,numEvents,params[0],params[1],params[2],params[3]))
            exit(0)

runAllSourceSims()




makeLinearityParamFile()




if __name__ == "__main__":

    parser = OptionParser()
    parser.add_option("--createOctetLists",dest="createOctetLists",action="store_true",default=False,
                      help="Creates Octet lists for individual octets stored in {dataPath}/octet_lists_MB")
    parser.add_option("--findPedestals",dest="findPedestals",action="store_true",default=False,
                      help="Finds the pedestals!")  
    parser.add_option("--replayPass1",dest="replayPass1",action="store_true",default=False,
                      help="Run replayPass1 on specified runs or octets (Applies cuts, reconstructs event positions, determines PID,Type,Side,etc.)")
    parser.add_option("--replayPass2",dest="replayPass2",action="store_true",default=False,
                      help="Run replayPass2 on specified runs or octets (Applies Bi pulser gain corrections")
    parser.add_option("--replayPass3",dest="replayPass3",action="store_true",default=False,
                      help="Run replayPass3 on specified runs or octets (Applies Xe position maps to scintillator response)")
    parser.add_option("--replayPass4",dest="replayPass4",action="store_true",default=False,
                      help="Run replayPass4 on specified runs or octets (Applies Energy Calibration)")
    parser.add_option("--runGainBismuth",dest="runGainBismuth",action="store_true",default=False,
                      help="Calculate and plot the Bi pulser spectra and track the gain relative to the reference run")
    parser.add_option("--makeDirectories",dest="makeDirectories",action="store_true",default=False,
                      help="Makes all the analysis directories.")
                      

    options, args = parser.parse_args()


    if options.makeDirectories:
        beta = BetaReplayManager()
        cal.makeAllDirectories()


 
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
    
