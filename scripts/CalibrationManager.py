#!/usr/bin/python

import os
import sys
from optparse import OptionParser
from math import *
import MButils

##### Set up list of runs which are to be omitted from the Energy Calibration
omittedRuns = [19232]
omittedRanges = [(17925,18055)] #These runs are from Run period 4 and include very long runs and runs with no Sn or Bi

for Range in omittedRanges:
    for run in range(Range[0],Range[1]+1,1):
        omittedRuns.append(run)

##### Set up which runs to omit which PMTs.
#### This will be done by having a file with all source runs, where the PMT value
#### is set to 0 or 1 to represent false (don't use) and true (do use)

EPMT1 = [] #These hold individual runs where PMT was flaky or Bi pulser was not working
EPMT2 = []
EPMT3 = []
EPMT4 = []
WPMT1 = []
WPMT2 = []
WPMT3 = []
WPMT4 = []

EPMT1_runRanges = [] #These hold chunks of runs where PMT is dead or Bi pulser is not working.
EPMT2_runRanges = []
EPMT3_runRanges = []
EPMT4_runRanges = []
WPMT1_runRanges = [(17359,18055)]
WPMT2_runRanges = [(16983,17297)]
WPMT3_runRanges = []
WPMT4_runRanges = [(18745,18768),(19347,19960)]

for Range in EPMT1_runRanges:
    for run in range(Range[0],Range[1]+1,1):
        EPMT1.append(run)
for Range in EPMT2_runRanges:
    for run in range(Range[0],Range[1]+1,1):
        EPMT2.append(run)
for Range in EPMT3_runRanges:
    for run in range(Range[0],Range[1]+1,1):
        EPMT3.append(run)
for Range in EPMT4_runRanges:
    for run in range(Range[0],Range[1]+1,1):
        EPMT4.append(run)
for Range in WPMT1_runRanges:
    for run in range(Range[0],Range[1]+1,1):
        WPMT1.append(run)
for Range in WPMT2_runRanges:
    for run in range(Range[0],Range[1]+1,1):
        WPMT2.append(run)
for Range in WPMT3_runRanges:
    for run in range(Range[0],Range[1]+1,1):
        WPMT3.append(run)
for Range in WPMT4_runRanges:
    for run in range(Range[0],Range[1]+1,1):
        WPMT4.append(run)


class CalReplayManager:
    
    def __init__(self,AnalysisType="MB"):
        self.AnalyzerPath = "../"
        self.runListPath = self.AnalyzerPath + "run_lists/"
        self.AnalysisDataPath = os.getenv("PARALLEL_DATA_PATH")
        self.srcPositionsPath = os.getenv("SOURCE_POSITIONS")
        self.srcPeakPath = os.getenv("SOURCE_PEAKS")
        self.replayPass3 = os.getenv("REPLAY_PASS3")
        self.srcListPath = os.getenv("SOURCE_LIST")

    def runReplayPass1(self,srcRunPeriod=1):
        print "Running replay_pass1 for run period %i"%srcRunPeriod
        filename = "Source_Calibration_Run_Period_%i.dat"%srcRunPeriod
        infile = open(self.runListPath+filename,'r')
        runs = []
        for line in infile:      
            runs.append(int(line))
        
        for run in runs:
            os.system("cd ../replay_pass1/; ./replay_pass1.exe %i"%run)
        print "DONE"
        
    def runReplayPass2(self,srcRunPeriod=1):
        print "Running replay_pass2 for run period %i"%srcRunPeriod
        filename = "Source_Calibration_Run_Period_%i.dat"%srcRunPeriod
        infile = open(self.runListPath+filename,'r')
        runs = []
        for line in infile:      
            runs.append(int(line))
        
        for run in runs:
            os.system("cd ../replay_pass2/; ./replay_pass2.exe %i"%run)
        print "DONE"
        
    def runReplayPass3(self,srcRunPeriod=1):
        print "Running replay_pass3 for run period %i"%srcRunPeriod
        filename = "Source_Calibration_Run_Period_%i.dat"%srcRunPeriod
        infile = open(self.runListPath+filename,'r')
        runs = []
        for line in infile:      
            runs.append(int(line))
        
        for run in runs:
            os.system("cd ../replay_pass3/; ./replay_pass3.exe %i"%run)
        print "DONE"

    def runReplayPass4(self,srcRunPeriod=1):
        print "Running replay_pass4 for run period %i"%srcRunPeriod
        os.system("mkdir -s %s"%(os.getenv("REPLAY_PASS4")))
        filename = "Source_Calibration_Run_Period_%i.dat"%srcRunPeriod
        infile = open(self.runListPath+filename,'r')
        runs = []
        for line in infile:      
            runs.append(int(line))
        
        for run in runs:
            os.system("cd ../replay_pass4/; ./replay_pass4.exe %i"%run)
        print "DONE"
        
    def runGainBismuth(self,srcRunPeriod=1):
        print "Running gain_bismuth for run period %i"%srcRunPeriod
        filename = "Source_Calibration_Run_Period_%i.dat"%srcRunPeriod
        infile = open(self.runListPath+filename,'r')
        runs = []
        for line in infile:      
            runs.append(int(line))
        
        for run in runs:
            os.system("cd ../gain_bismuth/; ./gain_bismuth.exe %i"%run)
            os.system("root -l -b -q '../gain_bismuth/plot_gain_bismuth.C(\"%i\")'"%run)
        print "DONE"


        
class CalibrationManager:
    
    def __init__(self,AnalysisType="MB"):
        self.AnalyzerPath = "../"
        self.runListPath = self.AnalyzerPath + "run_lists/"
        self.AnalysisDataPath = os.getenv("PARALLEL_DATA_PATH")
        self.srcPositionsPath = os.getenv("SOURCE_POSITIONS")
        self.srcPeakPath = os.getenv("SOURCE_PEAKS")
        self.replayPass3 = os.getenv("REPLAY_PASS3")
        self.srcListPath = os.getenv("SOURCE_LIST")
        self.nPEcountPath = os.getenv("NPE_WEIGHTS")


    def calc_nPE_per_PMT(self, runAllRefRun=False, run=19359):

        if runAllRefRun:
            srcRunPeriod=1
            srcRefRuns = [17238,17370,17521,17925,18361,18621,18749,19232,19359,19511,19857,19899]
            for r in srcRefRuns:
                os.system("cd ../calc_nPE/; ./calc_nPE.exe %i"%r)
                print "Ran calc_nPE.exe for run %i"%r
                runs=[]
                filename = "Source_Calibration_Run_Period_%i.dat"%srcRunPeriod
                infile = open(self.runListPath+filename,'r')
                for line in infile:
                    runs.append(int(line))
                for rn in runs:
                    os.system("cp %s/nPE_weights_%i.dat %s/nPE_weights_%i.dat"%(self.nPEcountPath,r,self.nPEcountPath,rn))
                srcRunPeriod+=1
                infile.close()
                            
                    
        else:
            os.system("cd ../calc_nPE/; ./calc_nPE.exe %i"%run)
            print "Ran calc_nPE.exe for run %i"%run
        
        
        

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


    def fitSourcePeaksInEnergy(self,srcRunPeriod=1, overwrite=True):
        filename = "Source_Calibration_Run_Period_%i.dat"%srcRunPeriod
        infile = open(self.runListPath+filename,'r')
        runs = []
        for line in infile:       
            #checking if peaks have already been fit
            filepath = self.srcPeakPath +"source_peaks_EnergyPeak_%i.dat"%int(line)
            print filepath
            if not MButils.fileExistsAndNotEmpty(filepath) or overwrite:
                #if int(line)>18734 and int(line)<18756: #TAKE THIS OUT AFTER DOING COMPARISON W/ BRADS ENVELOPE
                runs.append(int(line))

        for run in runs:
            os.system("cd ../source_peaks; ./source_peaks_EnergyPeak.exe %i"%run)
            os.system("root -b -q '../source_peaks/plot_source_peaks_Energy.C(\"%i\")'"%run)
            print "Ran fit_source_peaks_Energy.C on run %i"%run


    def makePMTrunFile(self,CalibrationPeriod=1, master=False):
        if not master:
            outputFile = "../residuals/PMT_runQuality_SrcPeriod_%i.dat"%(CalibrationPeriod)
            runList = []

            with open("../run_lists/Source_Calibration_Run_Period_%i.dat"%CalibrationPeriod) as runlist:
                for run in runlist:
                    if os.path.isfile(self.srcListPath+"source_list_%i.dat"%int(run)) and int(run) not in omittedRuns:
                        srcList = open(self.srcListPath+"source_list_%i.dat"%int(run))
                        lines = []
                        for line in srcList:
                            lines.append(line)
                            if int(lines[0])>0:
                                runList.append(int(run))

            outfile = open(outputFile,'w')

            for run in runList:
                pmtList = [1,1,1,1,1,1,1,1]
                if run in EPMT1:
                    pmtList[0]=0
                if run in EPMT2:
                    pmtList[1]=0
                if run in EPMT3:
                    pmtList[2]=0
                if run in EPMT4:
                    pmtList[3]=0
                if run in WPMT1:
                    pmtList[4]=0
                if run in WPMT2:
                    pmtList[5]=0
                if run in WPMT3:
                    pmtList[6]=0
                if run in WPMT4:
                    pmtList[7]=0

                outfile.write("%i %i %i %i %i %i %i %i %i\n"%(run,pmtList[0],pmtList[1],pmtList[2],pmtList[3],
                                                          pmtList[4],pmtList[5],pmtList[6],pmtList[7]))

            outfile.close()
            print "Done writing PMT file for Source Period %i"%CalibrationPeriod

        #Update the master list of PMT quality
        if master:
            masterFile = open("../residuals/PMT_runQuality_master.dat",'w')
            
            for run in range(16983,20000,1):
                pmtList = [1,1,1,1,1,1,1,1]
                if run in EPMT1:
                    pmtList[0]=0
                if run in EPMT2:
                    pmtList[1]=0
                if run in EPMT3:
                    pmtList[2]=0
                if run in EPMT4:
                    pmtList[3]=0
                if run in WPMT1:
                    pmtList[4]=0
                if run in WPMT2:
                    pmtList[5]=0
                if run in WPMT3:
                    pmtList[6]=0
                if run in WPMT4:
                    pmtList[7]=0

                masterFile.write("%i %i %i %i %i %i %i %i %i\n"%(run,pmtList[0],pmtList[1],pmtList[2],pmtList[3],
                                                          pmtList[4],pmtList[5],pmtList[6],pmtList[7]))
            masterFile.close()
            print "Updated master list of PMT run quality"



    def makeSourceCalibrationFile(self,CalibrationPeriod=1, PeaksInEnergy=False):
        #This utilizes the omittedRuns and removes them from the calibration. Any time you make a change
        # to the runs which are to be omitted, you shoud rerun this!

        outputFile = None
        src_file_base = None
        if PeaksInEnergy:
            outputFile = "../residuals/source_runs_EnergyPeaks_RunPeriod_%i.dat"%(CalibrationPeriod)
            src_file_base = self.srcPeakPath + "source_peaks_EnergyPeak_"
        else:
            outputFile = "../residuals/source_runs_RunPeriod_%i.dat"%(CalibrationPeriod)
            src_file_base = self.srcPeakPath + "source_peaks_"
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
            src_file = src_file_base+"%i.dat"%run
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
        if os.path.isfile(filename):
            os.system("root -b -q 'MB_calc_residuals.C (%i)'"%CalibrationPeriod)
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
                        if not line[0:3]=="PMT":
                            outfile.write(line)

            outfile.close()
        if PMT==0:
            print "Produced global residual file for weighted average of all PMTs, run periods %i-%i, and sides "%(periodLow,periodHigh), sides
        else:
            print "Produced global residual file for PMT %i, run periods %i-%i, and sides "%(PMT,periodLow,periodHigh), sides
    

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

class BetaDecayDataManager:
    def __init__(self):
        self.AnalyzerPath = "../"
        self.runListPath = self.AnalyzerPath + "run_lists/"
        self.AnalysisDataPath = os.getenv("PARALLEL_DATA_PATH")
        self.srcPositionsPath = os.getenv("SOURCE_POSITIONS")
        self.srcPeakPath = os.getenv("SOURCE_PEAKS")
        self.replayPass3 = os.getenv("REPLAY_PASS3")
        self.srcListPath = os.getenv("SOURCE_LIST")

    def runReplayPass4(self, run=None):
        return 0



if __name__ == "__main__":

    parser = OptionParser()
    parser.add_option("--fitSrcPeaks",dest="fitSrcPeaks",action="store_true",default=False,
                      help="Fit Source peaks for all runs unless otherwise stated")
    parser.add_option("--fitSrcPeaksInEnergy",dest="fitSrcPeaksInEnergy",action="store_true",default=False,
                      help="Fit Source peaks post replay_pass4 for all runs unless otherwise stated")
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
    parser.add_option("--makePMTrunFile",dest="makePMTrunFile",action="store_true",default=False,
                      help="Make file with booleans for whether to use each PMT for each run.")
                      

    options, args = parser.parse_args()


    ### This will fit all the source peaks for any runs in the runPeriods list below
    if options.fitSrcPeaks:
        
        runPeriods = [1,2,3,4,5,6,7,8,9,10,11,12]
        cal = CalibrationManager()
        for period in runPeriods:
            cal.fitSourcePeaks(period)

    if options.fitSrcPeaksInEnergy:
        
        runPeriods = [1]#,2,3,4,5,6,7,8,9,10,11,12]
        cal = CalibrationManager()
        for period in runPeriods:
            cal.fitSourcePeaksInEnergy(period)


    ### This will re-fit all the source positions. Note that by default, even with this option as true, 
    ### the overwrite option is False. These source positions have been checked by hand, and the only
    ### reason to redo them is if the position maps change.
    if options.fitSrcPositions:
        
        runPeriods = [1,2,3,4,5,6,7,8,9,10,11,12]
        cal = CalibrationManager()
        for period in runPeriods:
            cal.fitSourcePositions(period,False)
        
    ### Makes a file with each run followed by a boolean (0,1) for whether each PMT should be used or not
    if options.makePMTrunFile:
        runPeriods = [1,2,3,4,5,6,7,8,9,10,11,12]
        cal=CalibrationManager()
        for period in runPeriods:
            cal.makePMTrunFile(period)


    ### If you have made changes to the runs which are to be ignored at the top of this script, you should run this 
    if options.makeAllCalFiles:
    
        runPeriods = [1,2,3,4,5,6,7,8,9,10,11,12]
        cal = CalibrationManager()
        for period in runPeriods:
            cal.makeSourceCalibrationFile(period)
            cal.makePMTrunFile(period)
            #cal.calculateResiduals(period)
            

        #cal.makeGlobalResiduals(runPeriods,PMT=1,Side="Both")


    ### useful if you are going to look at the linearity curves and residuals by eye separately
    if options.makePeakFiles:
        runPeriods = [1,2,3,4,5,6,7,8,9,10,11,12]
        cal=CalibrationManager()
        for period in runPeriods:
            cal.makeSourceCalibrationFile(period, True)


    ### Saves the results of calculating the mean and RMS of all the global residuals for the 
    ### given combination of calibration periods and for a certain PMT (PMT=0 is for all 4 combined)
    if options.ErrorEnvelope:
        cal = CalibrationManager()
        cal.plotErrorEnvelope(calPeriodLow=1,calPeriodHigh=12,PMT=0)

    

    ## Makes file holding all the residuals for each PMT for each run which is to be used
    if options.makeGlobalResiduals:
        cal = CalibrationManager()
        runPeriods = [1,2,3,4,5,6,7,8,9,10,11,12]
        pmts = [0]#[1,2,3,4] #PMT 0 is for the weighted average of all 4
        for pmt in pmts:
            cal.makeGlobalResiduals(runPeriods,PMT=pmt,Side="Both")


    if 0:
        rep = CalReplayManager()
        cal = CalibrationManager()
        runPeriods = [1,2,3,4,5,6,7,8,9,10,11,12]
        for runPeriod in runPeriods:
            #rep.runReplayPass1(runPeriod)
            #rep.runGainBismuth(runPeriod)
            #rep.runReplayPass2(runPeriod)
            #rep.runReplayPass3(runPeriod)
            #cal.fitSourcePeaks(runPeriod)
            rep.runReplayPass4(runPeriod)

    if 0:
        cal = CalibrationManager()
        #cal.calc_nPE_per_PMT(True)
        cal.makePMTrunFile(master=True)


    #Trying to figure out why the east side isn't reconstructed as well after replay pass 4
    if 1: 
        runPeriods =  [1,2,3,4,5,6,7,8,9,10,11,12]
        rep = CalReplayManager()
        cal = CalibrationManager()
        for runPeriod in runPeriods:
            rep.runReplayPass4(runPeriod)
            cal.makeSourceCalibrationFile(runPeriod, True)
