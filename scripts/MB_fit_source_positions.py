#!/usr/bin/python

#runs the position maps if they do not already exist in "filepath"

import os
import MButils

overwrite = False ### whether to overwrite previous source position data files
                 ### This shouldn't be necessary unless changes are made
                 ### to the Wirechamber position reconstruction routine
src_positions_path = os.environ("SOURCE_POSITIONS")
runlist_path = "../run_lists/"

for p in range(1,12,1):
    filename = "Source_Calibration_Run_Period_%i.dat"%p
    infile = open(runlist_path+filename,'r')
    runs = []
    for line in infile:
        #checking if positions have already been fit
        filepath = src_positions_path + "source_positions_%i.dat"%int(line)
        print filepath
        if not MButils.fileExistsAndNotEmpty(filepath) or overwrite:
            runs.append(int(line))

    for run in runs:
        os.system("root -b -q '../source_positions/fit_source_positions.C(\"%i\")'"%run)
        print "Running fit_source_positions.C on run %i"%run




    

