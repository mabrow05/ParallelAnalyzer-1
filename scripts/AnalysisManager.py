#!/usr/bin/python

import os
from math import *
import MButils.py

class CalibrationManager:
    
    def __init__(self,AnalysisType="MB"):
        self.AnalyzerPath = "../"
        self.AnalysisDataPath = "/extern/UCNA/"
        self.runListPath = self.AnalyzerPath + "run_lists/"
        self.srcPositionsPath = self.AnalysisDataPath + "source_positions_MB/"
        self.srcPeakPath = self.AnalysisDataPath + "source_peaks_MB/"


    def fitSourcePositions(srcRunPeriod=1, overwrite=False):
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
        

    def fitSourcePeaks(srcRunPeriod=1, overwrite=True):
        filename = "Source_Calibration_Run_Period_%i.dat"%srcRunPeriod
        infile = open(runlist_path+filename,'r')
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
            print "Running fit_source_peaks.C on run %i"%run


    def makeSourceCalibrationFile(rmin=0, rmax=0):
        outputFile=AnalysisDataPath+"CalibrationPlots_MB/source_runs_%i-%i.dat"%(rmin,rmax)
        outfile = open(outputFile,'w')
        for run in range(rmin,rmax+1,1):
            src_file = srcPeakPath + "source_peaks_%i.dat"%run
            if MButils.fileExistsAndNotEmpty(src_file):
                infile = open(src_file,'r')
                for line in infile:
                    outfile.write(line)
                infile.close()
        outfile.close()

    def calculateResiduals(rmin=0, rmax=0):
        ## need to make MB_calc_residuals.C into a function which takes arguments.
        ## for now, this needs to be opened and run by hand
        return 0
    
    def plotErrorEnvelope(rmin=0, rmax=0):
        ## need to make MB_errorEnvelope.C into a function which takes arguments
        return 0

        



if __name__ == "__main__":

