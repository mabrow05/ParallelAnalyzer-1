#!/usr/bin/python

#Runs source peak fitting routine. Set the overwrite var to False if you
# do not want to overwrite previously fitted peaks

import os
import MButils


overwrite = True ### whether to overwrite previous source peak data files
src_peak_path = os.environ["SOURCE_PEAKS"]
runlist_path = "../run_lists/"

for p in range(1,12,1):
    filename = "Source_Calibration_Run_Period_%i.dat"%p
    infile = open(runlist_path+filename,'r')
    runs = []
    for line in infile:
        
        #checking if peaks have already been fit
        filepath = src_peak_path +"source_peaks_%i.dat"%int(line)
        print filepath
        if not MButils.fileExistsAndNotEmpty(filepath) or overwrite:
            #if int(line)>17949 and int(line)<20000: #TAKE THIS OUT AFTER DOING COMPARISON W/ BRADS ENVELOPE
                runs.append(int(line))

    for run in runs:
        #os.system("cd ../source_peaks; ./source_peaks.exe %i"%run)
        os.system("root -b -q '../source_peaks/plot_source_peaks.C(\"%i\")'"%run)
        print "Running fit_source_peaks.C on run %i"%run
