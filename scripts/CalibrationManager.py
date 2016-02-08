#!/usr/bin/python

########################################################################
# This code is meant to control the replaying and energy calibration of 
# the source data for a given data set. It is currently set up for the 
# 2011/2012 data set and needs to be adapted for 2012/2013 unless we 
# enumerate these source runs starting with the next number in line after
# 2011/2012...
########################################################################

import os
import sys
from optparse import OptionParser
from math import *
import MButils

##### Set up list of runs which are to be omitted from the Energy Calibration
omittedRuns = [19232, 20529, 20530, 20531, 20823, 20824, 20825, 20826, 20827, 21097]
# 20529 - very low statistics
# 20530,20531,20823-20827 - lost West event triggers
omittedRanges = [(17942,18055), (20901, 20917)] #These runs are from Run period 4 and include very long runs and runs with no Sn or Bi
                                               # And also run period 15 which is useless

for Range in omittedRanges:
    for run in range(Range[0],Range[1]+1,1):
        omittedRuns.append(run)

##### Set up which runs to omit which PMTs.
#### This will be done by having a file with all source runs, where the PMT value
#### is set to 0 or 1 to represent false (don't use) and true (do use)

EPMT1 = [] #These hold individual runs where PMT was flaky or Bi pulser was not working. 
EPMT2 = []
EPMT3 = []
EPMT4 = [20517,20519,20821,20822]
WPMT1 = []
WPMT2 = []
WPMT3 = []
WPMT4 = []

EPMT1_runRanges = [] #These hold chunks of runs where PMT is dead or Bi pulser is not working.
EPMT2_runRanges = []
EPMT3_runRanges = []
EPMT4_runRanges = [(17233,18055), (20121,23173)]
WPMT1_runRanges = [(17359,18055)]
WPMT2_runRanges = [(16983,17297)]
WPMT3_runRanges = []
WPMT4_runRanges = [(18370,18386),(18745,18768),(19347,19960),(20000,23000)]

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
        self.replayPass1 = os.getenv("REPLAY_PASS1")
        self.replayPass2 = os.getenv("REPLAY_PASS2")
        self.replayPass3 = os.getenv("REPLAY_PASS3")
        self.replayPass4 = os.getenv("REPLAY_PASS4")
        self.srcListPath = os.getenv("SOURCE_LIST")
        self.gainBismuthPath = os.getenv("GAIN_BISMUTH")
        self.nPEweightsPath = os.getenv("NPE_WEIGHTS")
        self.octetListPath = os.getenv("OCTET_LIST")
        self.triggerFuncPath = os.getenv("TRIGGER_FUNC")
        self.revCalSimPath = os.getenv("REVCALSIM")
        self.UKspecReplayPath = os.getenv("UK_SPEC_REPLAY")
        self.AnalysisResultsPath = os.getenv("ANALYSIS_RESULTS")

    def makeAllDirectories(self):
        #os.system("mkdir -p %s"%self.AnalysisDataPath)
        os.system("mkdir -p %s"%self.srcPositionsPath)
        os.system("mkdir -p %s"%os.getenv("BASIC_HISTOGRAMS"))
        os.system("mkdir -p %s"%os.getenv("PEDESTALS"))
        os.system("mkdir -p %s"%self.srcPeakPath)
        os.system("mkdir -p %s"%self.srcListPath)
        os.system("mkdir -p %s"%self.replayPass1)
        os.system("mkdir -p %s"%self.replayPass2)
        os.system("mkdir -p %s"%self.replayPass3)
        os.system("mkdir -p %s"%self.replayPass4)
        os.system("mkdir -p %s"%self.gainBismuthPath)
        os.system("mkdir -p %s"%self.nPEweightsPath)
        os.system("mkdir -p %s"%self.octetListPath)
        os.system("mkdir -p %s"%self.triggerFuncPath)
        os.system("mkdir -p %s"%self.revCalSimPath)
        os.system("mkdir -p %s"%self.UKspecReplayPath)
        os.system("mkdir -p %s"%self.AnalysisResultsPath)

    def runReplayPass1(self,srcRunPeriod=1, sourceORxenon="source"):
        print "Running replay_pass1 for %s run period %i"%(sourceORxenon,srcRunPeriod)
        filename=None
        if sourceORxenon=="source":
            filename = "Source_Calibration_Run_Period_%i.dat"%srcRunPeriod
        elif sourceORxenon=="xenon":
            filename = "Xenon_Calibration_Run_Period_%i.dat"%srcRunPeriod
        else:
            print "Not a valid source type!! Options: source or xenon"
            exit();
        infile = open(self.runListPath+filename,'r')
        runs = []
        for line in infile:      
            runs.append(int(line))
        
        for run in runs:
            os.system("cd ../replay_pass1/; ./replay_pass1.exe %i"%run)
        print "DONE"
        
    def runReplayPass2(self,srcRunPeriod=1, sourceORxenon="source"):
        print "Running replay_pass2 for %s run period %i"%(sourceORxenon,srcRunPeriod)
        filename=None
        if sourceORxenon=="source":
            filename = "Source_Calibration_Run_Period_%i.dat"%srcRunPeriod
        elif sourceORxenon=="xenon":
            filename = "Xenon_Calibration_Run_Period_%i.dat"%srcRunPeriod
        else:
            print "Not a valid source type!! Options: source or xenon"
            exit();
        infile = open(self.runListPath+filename,'r')
        runs = []
        for line in infile:      
            runs.append(int(line))
        
        for run in runs:
            os.system("cd ../replay_pass2/; ./replay_pass2.exe %i"%run)
        print "DONE"
        
    def runReplayPass3(self,srcRunPeriod=1, sourceORxenon="source"):
        print "Running replay_pass3 for %s run period %i"%(sourceORxenon,srcRunPeriod)
        filename=None
        if sourceORxenon=="source":
            filename = "Source_Calibration_Run_Period_%i.dat"%srcRunPeriod
        elif sourceORxenon=="xenon":
            filename = "Xenon_Calibration_Run_Period_%i.dat"%srcRunPeriod
        else:
            print "Not a valid source type!! Options: source or xenon"
            exit();
        infile = open(self.runListPath+filename,'r')
        runs = []
        for line in infile:      
            runs.append(int(line))
        
        for run in runs:
            os.system("cd ../replay_pass3/; ./replay_pass3.exe %i"%run)
        print "DONE"

    def runReplayPass4(self,srcRunPeriod=1, sourceORxenon="source"):
        print "Running replay_pass4 for %s run period %i"%(sourceORxenon,srcRunPeriod)
        os.system("mkdir -p %s"%(os.getenv("REPLAY_PASS4")))
        filename=None
        if sourceORxenon=="source":
            filename = "Source_Calibration_Run_Period_%i.dat"%srcRunPeriod
        elif sourceORxenon=="xenon":
            filename = "Xenon_Calibration_Run_Period_%i.dat"%srcRunPeriod
        else:
            print "Not a valid source type!! Options: source or xenon"
            exit();
        infile = open(self.runListPath+filename,'r')
        runs = []
        for line in infile:      
            runs.append(int(line))
        
        for run in runs:
            os.system("cd ../replay_pass4/; ./replay_pass4.exe %i"%run)
        print "DONE"
        
    def runGainBismuth(self,srcRunPeriod=1, sourceORxenon="source"):
        print "Running gain_bismuth for %s run period %i"%(sourceORxenon,srcRunPeriod)
        filename=None
        if sourceORxenon=="source":
            filename = "Source_Calibration_Run_Period_%i.dat"%srcRunPeriod
        elif sourceORxenon=="xenon":
            filename = "Xenon_Calibration_Run_Period_%i.dat"%srcRunPeriod
        else:
            print "Not a valid source type!! Options: source or xenon"
            exit();
        infile = open(self.runListPath+filename,'r')
        runs = []
        for line in infile:      
            runs.append(int(line))
        
        for run in runs:
            os.system("cd ../gain_bismuth/; ./gain_bismuth.exe %i"%run)
            os.system("root -l -b -q '../gain_bismuth/plot_gain_bismuth.C(\"%i\")'"%run)
        print "DONE"


    def makeBasicHistograms(self, srcRunPeriod=1, sourceORxenon="source"):
        print "Making Basic Histograms for %s run period %i"%(sourceORxenon,srcRunPeriod)
        filename=None
        if sourceORxenon=="source":
            filename = "Source_Calibration_Run_Period_%i.dat"%srcRunPeriod
        elif sourceORxenon=="xenon":
            filename = "Xenon_Calibration_Run_Period_%i.dat"%srcRunPeriod
        else:
            print "Not a valid source type!! Options: source or xenon"
            exit();
        infile = open(self.runListPath+filename,'r')
        runs = []
        for line in infile:      
            runs.append(int(line))
        
        for run in runs:
            print "Making basic histograms for run %i"%run
            os.system("cd ../basic_histograms/; ./basic_histograms.exe %i"%run)
            os.system("root -l -b -q '../basic_histograms/plot_basic_histograms.C(\"%i\")'"%run)
            
        print "DONE"


    def findPedestals(self, srcRunPeriod=1, sourceORxenon="source"):
        print "Running pedestals for %s run period %i"%(sourceORxenon,srcRunPeriod)
        filename=None
        if sourceORxenon=="source":
            filename = "Source_Calibration_Run_Period_%i.dat"%srcRunPeriod
        elif sourceORxenon=="xenon":
            filename = "Xenon_Calibration_Run_Period_%i.dat"%srcRunPeriod
        else:
            print "Not a valid source type!! Options: source or xenon"
            exit();
        infile = open(self.runListPath+filename,'r')
        runs = []
        for line in infile:      
            runs.append(int(line))
        
        for run in runs:
            os.system("cd ../pedestals/; ./pedestals.exe %i"%run)
            
        print "DONE"


    def runReverseCalibration(self, srcRunPeriod=1, sourceORxenon="source"):
        print "Running reverse calibration for %s run period %i"%(sourceORxenon,srcRunPeriod)
        filename=None
        if sourceORxenon=="source":
            filename = "Source_Calibration_Run_Period_%i.dat"%srcRunPeriod
        elif sourceORxenon=="xenon":
            filename = "Xenon_Calibration_Run_Period_%i.dat"%srcRunPeriod
        else:
            print "Not a valid source type!! Options: source or xenon"
            exit();
        infile = open(self.runListPath+filename,'r')
        runs = []
        for line in infile:      
            runs.append(int(line))

        for run in runs:
            filename = self.srcListPath+"source_list_%i.dat"%run
            if not MButils.fileExistsAndNotEmpty(filename):
                continue
            srcFile = open(filename)
            srcFileInput = []
            for line in srcFile:
                srcFileInput.append(line)
            for src in range(1,int(srcFileInput[0])+1):
                source = srcFileInput[src][0:2]
                if source=="Ce":
                    source = source+"139"
                elif source=='Sn':
                    source = source + "113"
                elif source=="Bi":
                    source = source + "207"
                elif source == "In":
                    source = source + "114"
                elif source == "Cd":
                    source = source + "109"
                elif source == "Cs":
                    source = source+"137"

                os.system("./../simulation_comparison/revCalSim.exe %i %s"%(run, source))
                #print "./../simulation_comparison/revCalSim.exe %i %s"%(run, source)

        print "Finished reverse calibration for " + sourceORxenon + "run period %i"%srcRunPeriod


        
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
        self.revCalSimPath = os.getenv("REVCALSIM")

    def findCalibrationPeriod(self, runNumber):
        if runNumber <= 17297:
            return 1
        elif runNumber <= 17439:
            return 2
        elif runNumber <= 17734:
            return 3
        elif runNumber <= 17955:
            return 4
        elif runNumber <= 18386:
            return 5
        elif runNumber <= 18683:
            return 6
        elif runNumber <= 18994:
            return 7
        elif runNumber <= 19239:
            return 8
        elif runNumber <= 19544:
            return 9
        elif runNumber <= 20000:
            return 11
        elif runNumber <= 20741:
            return 13
        elif runNumber <= 20837:
            return 14
        elif runNumber <= 21237:
            return 16
        elif runNumber <= 21605:
            return 17
        elif runNumber <= 21863:
            return 18
        elif runNumber <= 22118:
            return 19
        elif runNumber <= 22238:
            return 20
        elif runNumber <= 22630:
            return 22
        elif runNumber <= 23173:
            return 23
        else:
            print "Bad Run Number!!! No Calibration Period..."
            exit



    def calc_nPE_per_PMT(self, runAllRefRun=False, run=19359, writeNPEforAllRuns=False, year="2011-2012"):

        srcSn_nPE_Runs = [17238,17370,17521,17925,18361,18621,18749,19232,19359,19511,19857,19899,20520,20823,20905,21091,21315,21683,21918,22219,22298,22441,22771,22925] #These are runs for which the nPE per channel are calculated
        # May or may not match the reference runs

        srcRunPeriodRange = None
        runRange = None

        if year=="2011-2012":
            srcRunPeriodRange = [1,12]
            runRange = [16000,20000]

        elif year=="2012-2013":
            srcRunPeriodRange = [13,24]
            runRange = [20000,23173]

        if runAllRefRun:
            
            for srcRunPeriod in range(srcRunPeriodRange[0],srcRunPeriodRange[1]+1,1):
                r = srcSn_nPE_Runs[srcRunPeriod-1]
                os.system("cd ../calc_nPE/; ./calc_nPE.exe %i"%r)
                print "Ran calc_nPE.exe for run %i"%r
                #runs=[]
                #filename = "Source_Calibration_Run_Period_%i.dat"%srcRunPeriod
                #infile = open(self.runListPath+filename,'r')
                #for line in infile:
                #    runs.append(int(line))
                #for rn in runs:
                #    os.system("cp %s/nPE_weights_%i.dat %s/nPE_weights_%i.dat"%(self.nPEcountPath,r,self.nPEcountPath,rn))
                #infile.close()

        elif writeNPEforAllRuns:
            for rn in range(runRange[0],runRange[1]+1,1):
                calPeriod = self.findCalibrationPeriod(rn)
                os.system("cp %s/nPE_weights_%i.dat %s/nPE_weights_%i.dat"%(self.nPEcountPath,srcSn_nPE_Runs[calPeriod-1],self.nPEcountPath,rn))

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
                runs.append(int(line))

        for run in runs:
            os.system("cd ../source_peaks; ./source_peaks.exe %i"%run)
            os.system("root -b -q '../source_peaks/plot_source_peaks.C(\"%i\")'"%run)
            print "Ran fit_source_peaks.C on run %i"%run


    def fitSourcePeaksInEnergy(self,srcRunPeriod=1, PMTbyPMT=False, Simulation=False, overwrite=True):
        filename = "Source_Calibration_Run_Period_%i.dat"%srcRunPeriod
        infile = open(self.runListPath+filename,'r')
        runs = []
        for line in infile:       
            #checking if peaks have already been fit
            filepath = self.srcPeakPath +"source_peaks_EvisPMTbyPMT_%i.dat"%int(line)
            if not PMTbyPMT:
                filepath = self.srcPeakPath +"source_peaks_EnergyPeak_%i.dat"%int(line)
            
            print filepath
            if not MButils.fileExistsAndNotEmpty(filepath) or overwrite:
                runs.append(int(line))

        for run in runs:
            if PMTbyPMT:
                if Simulation:
                     os.system("cd ../source_peaks; ./sim_source_peaks_EvisPMTbyPMT.exe %i"%run)
                     os.system("root -b -q '../source_peaks/plot_sim_source_peaks_Evis.C(\"%i\")'"%run)
                else:   
                    os.system("cd ../source_peaks; ./source_peaks_EvisPMTbyPMT.exe %i"%run)
                    os.system("root -b -q '../source_peaks/plot_source_peaks_Evis.C(\"%i\")'"%run)
            else:
                if Simulation:
                    os.system("cd ../source_peaks; ./sim_source_peaks_EnergyPeak.exe %i"%run)
                    #os.system("root -b -q '../source_peaks/plot_source_peaks_Energy.C(\"%i\")'"%run)
                else:
                    os.system("cd ../source_peaks; ./source_peaks_EnergyPeak.exe %i"%run)
                    os.system("root -b -q '../source_peaks/plot_source_peaks_Energy.C(\"%i\")'"%run)
            print "Ran plot_source_peaks.C on run %i"%run


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
            
            for run in range(16983,23174,1):
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



    def makeSourceCalibrationFile(self,CalibrationPeriod=1, PeaksInEnergy=False, PMTbyPMT=False, Simulation=False):
        #This utilizes the omittedRuns and removes them from the calibration. Any time you make a change
        # to the runs which are to be omitted, you shoud rerun this!

        outputFile = None
        src_file_base = None
        if PeaksInEnergy and PMTbyPMT:
            if Simulation:
                outputFile = "../residuals/SIM_source_runs_EvisPMTbyPMT_RunPeriod_%i.dat"%(CalibrationPeriod)
                src_file_base = self.revCalSimPath + "/source_peaks/source_peaks_EvisPMTbyPMT_"
            else:
                outputFile = "../residuals/source_runs_EvisPMTbyPMT_RunPeriod_%i.dat"%(CalibrationPeriod)
                src_file_base = self.srcPeakPath + "source_peaks_EvisPMTbyPMT_"               
        elif PeaksInEnergy and not PMTbyPMT:
            if Simulation:
                outputFile = "../residuals/SIM_source_runs_EnergyPeaks_RunPeriod_%i.dat"%(CalibrationPeriod)
                src_file_base = self.revCalSimPath + "/source_peaks/source_peaks_EnergyPeak_"
            else:
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


    #### Calculate the linearity curves for each PMT from a source run period
    def LinearityCurves(self,CalibrationPeriod=1):
        filename = "../residuals/source_runs_RunPeriod_%i.dat"%CalibrationPeriod #Files which hold fitted source peaks
        if os.path.isfile(filename):
            os.system("root -b -q 'LinearityCurves.C (%i)'"%CalibrationPeriod)
            print "Calculated linearity curves for Calibration Period %i"%CalibrationPeriod
        else:
            print "No peak file to calculate linearity curves"
            sys.exit


    #Calculates residuals for each PMT and as a whole for given run period
    def calculateResiduals(self, CalibrationPeriod=1, PMTbyPMT=False):
        filename = None
        if PMTbyPMT: #Residuals are calculated for each PMT, then the weighted average of the residuals is plotted
            filename = "../residuals/source_runs_EvisPMTbyPMT_RunPeriod_%i.dat"%CalibrationPeriod
        else: #Residuals are calculated based on the final weighted energy peak
            filename = "../residuals/source_runs_EnergyPeaks_RunPeriod_%i.dat"%CalibrationPeriod

        if os.path.isfile(filename):
            if PMTbyPMT:
                os.system("root -b -q 'MB_calc_residuals_EvisPMTbyPMT.C (%i)'"%CalibrationPeriod)
            else:
                os.system("root -b -q 'MB_calc_residuals_finalEnergyFits.C (%i)'"%CalibrationPeriod)
            
            print "Calculated residuals for Calibration Period %i"%CalibrationPeriod
        else:
            print "No peak file to calculate residuals"
            sys.exit

    #Combines residuals for the calPeriods in the list given, for the PMT chosen (0 is PMT as a whole), and for the side
    # which can be "East", "West", or "Both"
    def makeGlobalResiduals(self,CalPeriods=[1], PMT=0, Side="Both",InEnergy=True, PMTbyPMT=False):
        CalPeriods.sort()
        periodLow = CalPeriods[0]
        periodHigh = CalPeriods[len(CalPeriods)-1]

        sides = [] #holds what sides will be run
        if Side=="Both":
            sides = ["East","West"]
        else:
            sides=[Side]

        for side in sides:
            outfile=None
            if (PMT==0):
                if InEnergy and PMTbyPMT:
                    outfile = open("../residuals/residuals_global_EvisPMTbyPMT_%s_periods_%i-%i.dat"%(side,periodLow,periodHigh),"w")
                elif InEnergy and not PMTbyPMT:
                    outfile = open("../residuals/residuals_global_EnergyPeaks_%s_periods_%i-%i.dat"%(side,periodLow,periodHigh),"w")
                else:
                    outfile = open("../residuals/residuals_global_%s_periods_%i-%i.dat"%(side,periodLow,periodHigh),"w")
            elif (side=="East" and PMT>0):
                outfile = open("../residuals/residuals_global_%s_periods_%i-%i_PMTE%i.dat"%(side,periodLow,periodHigh,PMT),"w")
            elif (side=="West" and PMT>0):
                outfile = open("../residuals/residuals_global_%s_periods_%i-%i_PMTW%i.dat"%(side,periodLow,periodHigh,PMT),"w")

            for period in CalPeriods:
                filename=None
                if (PMT==0):
                    if InEnergy and PMTbyPMT:
                        filename = "../residuals/residuals_EvisPMTbyPMT_%s_runPeriod_%i.dat"%(side,period)
                    elif InEnergy and not PMTbyPMT:
                        filename = "../residuals/residuals_%s_EnergyPeaks_runPeriod_%i.dat"%(side,period)
                    else:
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
    

    def plotErrorEnvelope(self, calPeriodLow=2, calPeriodHigh=10, PMT=0, InEnergy=False):
        ## This runs code which calculates the mean and RMS of the global residual file
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
    parser.add_option("--LinearityCurves",dest="LinearityCurves",action="store_true",default=False,
                      help="Run script which calculates linearity curves for PMTs.")
    parser.add_option("--makeDirectories",dest="makeDirectories",action="store_true",default=False,
                      help="Makes all the analysis directories.")
                      

    options, args = parser.parse_args()


    if options.makeDirectories:
        cal = CalReplayManager()
        cal.makeAllDirectories()


    ### This will fit all the source peaks for any runs in the runPeriods list below
    if options.fitSrcPeaks:
        
        runPeriods = [1,2,3,4,5,6,7,8,9,10,11,12]
        cal = CalibrationManager()
        for period in runPeriods:
            cal.fitSourcePeaks(period)

    if options.fitSrcPeaksInEnergy:
        
        runPeriods = [1,2,3,4,5,6,7,8,9,10,11,12]
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
        runPeriods = [13,14,15,16,17,18,19,20,21,22,23,24]#[1,2,3,4,5,6,7,8,9,10,11,12]
        cal=CalibrationManager()
        for period in runPeriods:
            cal.makePMTrunFile(period)
            
        cal.makePMTrunFile(1,True) #Updates the master list of PMT quality over all runs


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
        runPeriods = [13,14,16,17,18,19,20,21,22,23]#[1,2,3,4,5,6,7,8,9,10,11,12]
        cal=CalibrationManager()
        for period in runPeriods:
            cal.makeSourceCalibrationFile(period, False)


    ### Does the linearity curves for all the source calibration periods
    if options.LinearityCurves:
        runPeriods =[1,2,3,4,5,6,7,8,9,10,11,12]
        cal = CalibrationManager()
        for runPeriod in runPeriods:
            cal.LinearityCurves(runPeriod)


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
            cal.makeGlobalResiduals(runPeriods,PMT=pmt,Side="Both", InEnergy=True, PMTbyPMT=False)


    #### All the steps for completely replaying runs (without doing a new calibration or new position maps along the way)
    if 0:
        rep = CalReplayManager()
        cal = CalibrationManager()
        runPeriods = [14]#[13,14,16,17,18,19,20,21,22,23,24]# [1,2,3,4,5,6,7,8,9,10,11,12]
        for runPeriod in runPeriods:
            #rep.makeBasicHistograms(runPeriod, sourceORxenon="source")
            #rep.findPedestals(runPeriod)
            #rep.runReplayPass1(runPeriod)
            #rep.runGainBismuth(runPeriod)
            #rep.runReplayPass2(runPeriod)
            #rep.runReplayPass3(runPeriod)
            #cal.fitSourcePositions(runPeriod)
            #cal.fitSourcePeaks(runPeriod)
            cal.makeSourceCalibrationFile(runPeriod, False)
            #rep.runReplayPass4(runPeriod)

    ### Making the files which hold the PMT quality
    if 0:
        cal = CalibrationManager()
        #cal.calc_nPE_per_PMT(runAllRefRun=False,writeNPEforAllRuns=True,year="2012-2013")
        cal.makePMTrunFile(master=True)

    ### Simulation reverse calibration procedure
    if 0: 
        runPeriods =[13,14,16,17,18,19,20,21,22,23,24]#[1,2,3,4,5,6,7,8,9,10,11,12]
        rep = CalReplayManager()
        cal = CalibrationManager()
        
        for runPeriod in runPeriods:
            rep.runReverseCalibration(runPeriod)
            cal.fitSourcePeaksInEnergy(runPeriod, True, Simulation=True)
            cal.makeSourceCalibrationFile(runPeriod, PeaksInEnergy=True, PMTbyPMT=True, Simulation=True)

    ### Source Run Calibration Steps...
    if 1: 
        runPeriods = [13,14,16,17,18,19,20,21,22,23,24]#[1,2,3,4,5,6,7,8,9,10,11,12]#[5,6,7,8,9,10,11]#
        rep = CalReplayManager()
        cal = CalibrationManager()
        
        #for runPeriod in runPeriods:
            #cal.makeSourceCalibrationFile(period, False)
            #cal.makeSourceCalibrationFile(runPeriod, PeaksInEnergy=True, PMTbyPMT=True, Simulation=True)
            #cal.LinearityCurves(runPeriod)
            #rep.runReplayPass4(runPeriod)
            #cal.fitSourcePeaksInEnergy(runPeriod, PMTbyPMT=True, Simulation=False)
            #cal.makeSourceCalibrationFile(runPeriod, PeaksInEnergy=True, PMTbyPMT=True, Simulation=False)
            #cal.calculateResiduals(runPeriod, PMTbyPMT=True)
            
        cal.makeGlobalResiduals(runPeriods,PMT=0,Side="Both",InEnergy=True, PMTbyPMT=True)

    ### Replaying Xe Runs. Note that the position maps are calculated post replayPass2 and only need to
    ### be done once unless fundamental changes to the code are made upstream
    if 0: 
        runPeriods = [8,9,10]#[2,3,4,5,7] #### 1-7 are from 2011/2012, while 8-10 are from 2012/2013
        rep = CalReplayManager()
        cal = CalibrationManager()
        #cal.calc_nPE_per_PMT(runAllRefRun=False,writeNPEforAllRuns=True)
        for runPeriod in runPeriods:    
            #rep.makeBasicHistograms(runPeriod, sourceORxenon="xenon")
            #rep.findPedestals(runPeriod, sourceORxenon="xenon")
            #rep.runReplayPass1(runPeriod, sourceORxenon="xenon")
            #rep.runGainBismuth(runPeriod, sourceORxenon="xenon")
            #rep.runReplayPass2(runPeriod, sourceORxenon="xenon")
            #rep.runReplayPass3(runPeriod, sourceORxenon="xenon")
            rep.runReplayPass4(runPeriod, sourceORxenon="xenon")

            
    
