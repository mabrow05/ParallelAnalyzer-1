#!/usr/bin/python

import os
import sys
from optparse import OptionParser
from math import *
import MButils

##### Set up list of runs which are to be omitted from the Energy Calibration
omittedRuns = [17381,17383,17385,17382,17886,17912,19232]
omittedRanges = [(19347,19364),(18020,18055)]

for Range in omittedRanges:
    for run in range(Range[0],Range[1]+1,1):
        omittedRuns.append(run)


class CalibrationManager:
    
    def __init__(self,AnalysisType="MB"):
        self.AnalyzerPath = "../"
        self.runListPath = self.AnalyzerPath + "run_lists/"
        self.AnalysisDataPath = os.getenv("PARALLEL_DATA_PATH")
        self.srcPositionsPath = os.getenv("SOURCE_POSITIONS")
        self.srcPeakPath = os.getenv("SOURCE_PEAKS")
        self.replayPass3 = os.getenv("REPLAY_PASS3")
        self.srcListPath = os.getenv("SOURCE_LIST")

    def fitSourcePositions(self,srcRunPeriod=1, overwrite=False):
        filename = "Source_Calibration_Run_Period_%i.dat"%srcRunPeriod
        infile = open(self.runListPath + filename,'r')
        runs = []
        for line in infile:
        #checking if positions have already been fit
            filepath = self.srcPositionsPath + "source_positions_%i.dat"%int(line)
            print filepath
            if not MButils.fileExistsAndNotEmpty(filepath) or overwrite:
                runs.append(int(line))

        for run in runs:
            os.system("root -b -q '../source_positions/fit_source_positions.C(\"%i\")'"%run)
            print "Running fit_source_positions.C on run %i"%run
        

    def fitSourcePeaks(self,srcRunPeriod=1, overwrite=True):
        filename = "Source_Calibration_Run_Period_%i.dat"%srcRunPeriod
        infile = open(self.runListPath+filename,'r')
        runs = []
        for line in infile:       
            #checking if peaks have already been fit
            filepath = self.srcPeakPath +"source_peaks_%i.dat"%int(line)
            print filepath
            if not MButils.fileExistsAndNotEmpty(filepath) or overwrite:
                #if int(line)>18734 and int(line)<18756: #TAKE THIS OUT AFTER DOING COMPARISON W/ BRADS ENVELOPE
                runs.append(int(line))

        for run in runs:
            os.system("cd ../source_peaks; ./source_peaks.exe %i"%run)
            os.system("root -b -q '../source_peaks/plot_source_peaks.C(\"%i\")'"%run)
            print "Ran fit_source_peaks.C on run %i"%run


    def makeSourceCalibrationFile(self,CalibrationPeriod=1):

        outputFile = "../residuals/source_runs_RunPeriod_%i.dat"%(CalibrationPeriod)
        runList = []

        with open("../run_lists/Source_Calibration_Run_Period_%i.dat"%CalibrationPeriod) as runlist:
            for run in runlist:
                if os.path.isfile(self.srcListPath+"source_list_%i.dat"%int(run)):
                    srcList = open(self.srcListPath+"source_list_%i.dat"%int(run))
                    lines = []
                    for line in srcList:
                        lines.append(line)
                    if int(lines[0])>0:
                        runList.append(int(run))

        #print runList

        outfile = open(outputFile,'w')

        for run in runList:
            src_file = self.srcPeakPath + "source_peaks_%i.dat"%run
            if MButils.fileExistsAndNotEmpty(src_file) and run not in omittedRuns:
                infile = open(src_file,'r')
                for line in infile:
                    outfile.write(line)
                infile.close()

        outfile.close()
        print "Made combined source peak file for Calibration Period %i"%CalibrationPeriod

    
    #Calculates residuals for each PMT and as a whole for given run period
    def calculateResiduals(self, CalibrationPeriod=1):
        filename = "../residuals/source_runs_RunPeriod_%i.dat"%CalibrationPeriod
        if MButils.fileExistsAndNotEmpty(filename):
            os.system("root -l -b -q 'MB_calc_residuals.C (%i)'"%CalibrationPeriod)
            print "Calculated residuals for Calibration Period %i"%CalibrationPeriod
        else:
            print "No peak file to calculate residuals"
            sys.exit

    #Combines residuals for the calPeriods in the list given, for the PMT chosen (0 is PMT as a whole), and for the side
    # which can be "East", "West", or "Both"
    def makeGlobalResiduals(self,CalPeriods=[1], PMT=0, Side="Both"):
        CalPeriods.sort()
        periodLow = CalPeriods[0]
        periodHigh = CalPeriods[len(CalPeriods)-1]

        #for CalibrationPeriod in CalPeriods:

        sides = [] #holds what sides will be run
        if Side=="Both":
            sides = ["East","West"]
        else:
            sides=[Side]

        for side in sides:
            outfile=None
            if (PMT==0):
                outfile = open("../residuals/residuals_global_%s_periods_%i-%i.dat"%(side,periodLow,periodHigh),"w")
            elif (side=="East" and PMT>0):
                outfile = open("../residuals/residuals_global_%s_periods_%i-%i_PMTE%i.dat"%(side,periodLow,periodHigh,PMT),"w")
            elif (side=="West" and PMT>0):
                outfile = open("../residuals/residuals_global_%s_periods_%i-%i_PMTW%i.dat"%(side,periodLow,periodHigh,PMT),"w")

            for period in CalPeriods:
                filename=None
                if (PMT==0):
                    filename = "../residuals/residuals_%s_runPeriod_%i.dat"%(side,period)
                elif (side=="East" and PMT>0):
                    filename = "../residuals/residuals_%s_runPeriod_%i_PMTE%i.dat"%(side,period,PMT)
                elif (side=="West" and PMT>0):
                    filename = "../residuals/residuals_%s_runPeriod_%i_PMTW%i.dat"%(side,period,PMT)

                if os.path.isfile(filename):
                    resid = open(filename)
                    for line in resid:
                        outfile.write(line)

            outfile.close()
    

    def plotErrorEnvelope(self, calPeriodLow=2, calPeriodHigh=10, PMT=0):
        ## This runs code which calculates the mean and RMS of the global residual file from 
        ## from Calibration run periods CalPeriodLow to CalPeriodHigh. It prints out the mean and RMS
        ## for each source. Later this will be input into code which actually plots the error envelope.
        
        filenameEast = None
        filenameWest = None

        if calPeriodLow!=calPeriodHigh and not PMT:
            filenameEast = "../residuals/residuals_global_East_periods_%i-%i.dat"%(calPeriodLow, calPeriodHigh)
            filenameWest = "../residuals/residuals_global_West_periods_%i-%i.dat"%(calPeriodLow, calPeriodHigh)
        elif calPeriodLow==calPeriodHigh and not PMT:
            filenameEast = "../residuals/residuals_East_runPeriod_%i.dat"%calPeriodLow
            filenameWest = "../residuals/residuals_West_runPeriod_%i.dat"%calPeriodLow
        elif calPeriodLow!=calPeriodHigh and PMT:
            filenameEast = "../residuals/residuals_global_East_periods_%i-%i_PMTE%i.dat"%(calPeriodLow, calPeriodHigh, PMT)
            filenameWest = "../residuals/residuals_global_West_periods_%i-%i_PMTW%i.dat"%(calPeriodLow, calPeriodHigh, PMT)
        elif calPeriodLow==calPeriodHigh and PMT:
            filenameEast = "../residuals/residuals_East_runPeriod_%i_PMTE%i.dat"%(calPeriodLow,PMT)
            filenameWest = "../residuals/residuals_West_runPeriod_%i_PMTW%i.dat"%(calPeriodLow,PMT)
        
        if MButils.fileExistsAndNotEmpty(filenameEast) and MButils.fileExistsAndNotEmpty(filenameWest):
            print "Making Error Envelope for Run Periods %i to %i"%(calPeriodLow,calPeriodHigh)
            os.system("root -l -b -q 'MB_errorEnvelope.C (%i,%i,%i)'"%(calPeriodLow,calPeriodHigh,PMT))

        



if __name__ == "__main__":

    parser = OptionParser()
    parser.add_option("--fitSrcPeaks",dest="fitSrcPeaks",action="store_true",default=False,
                      help="Fit Source peaks for all runs unless otherwise stated")
    parser.add_option("--fitSrcPositions",dest="fitSrcPositions",action="store_true",default=False,
                      help="Fit for Source Positions for all runs. Only to be done if new replay_pass3 has been done on sources!")
    parser.add_option("--makeAllCalFiles",dest="makeAllCalFiles",action="store_true",default=False,
                      help="Recombines all source peaks for each Calibration period, then calculates residuals, and finally makes global residual file")
    parser.add_option("--makePeakFiles",dest="makePeakFiles",action="store_true",default=False,
                      help="Make combined source peak files for individual source calibration periods.")
    parser.add_option("--calcResiduals",dest="calcResiduals",action="store_true",default=False,
                      help="Calculate the residuals for the given run periods")
    parser.add_option("--makeGlobalResiduals",dest="makeGlobalResiduals",action="store_true",default=False,
                      help="Combine all residuals into one file for drawing error envelope")
    parser.add_option("--ErrorEnvelope",dest="ErrorEnvelope",action="store_true",default=False,
                      help="Make error envelope and save mean and sigma to file in ../error_envelope.")

    options, args = parser.parse_args()


    ### This will fit all the source peaks for any runs in the runPeriods list below
    if options.fitSrcPeaks:
        
        runPeriods = [2,3,4,5,6,7,8,9,10,11]
        cal = CalibrationManager()
        for period in runPeriods:
            cal.fitSourcePeaks(period)


    ### This will re-fit all the source positions. Note that by default, even with this option as true, 
    ### the overwrite option is False. These source positions have been checked by hand, and the only
    ### reason to redo them is if the position maps change.
    if options.fitSrcPositions:
        
        runPeriods = [2,3,4,5,6,7,8,9,10,11]
        cal = CalibrationManager()
        for period in runPeriods:
            cal.fitSourcePositions(period,False)
        

    ### If you have made changes to the runs which are to be ignored at the top of this script, you should run this 
    if options.makeAllCalFiles:
    
        runPeriods = [2,3,4,5,6,7,8,10]
        cal = CalibrationManager()
        for period in runPeriods:
            cal.makeSourceCalibrationFile(period)
            cal.calculateResiduals(period)

        cal.makeGlobalResiduals(runPeriods,PMT=1,Side="Both")


    ### useful if you are going to look at the linearity curves and residuals by eye separately
    if options.makePeakFiles:
        runPeriods = [7]
        cal=CalibrationManager()
        for period in runPeriods:
            cal.makeSourceCalibrationFile(period)


    ### Saves the results of calculating the mean and RMS of all the global residuals for the 
    ### given combination of calibration periods and for a certain PMT (PMT=0 is for all 4 combined)
    if options.ErrorEnvelope:
        cal = CalibrationManager()
        cal.plotErrorEnvelope(calPeriodLow=2,calPeriodHigh=10,PMT=1)
