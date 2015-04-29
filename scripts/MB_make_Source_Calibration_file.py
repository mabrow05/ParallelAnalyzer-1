#!/usr/bin/python

### Script to generate files with all source data included to produce 
### error envelope type plot

import os
import MButils

#runRange = [17359,19959]
CalibrationPeriod = 11
omittedRuns = [17383,17385,17382,17381,17521,17876,17886,17909,17912,17950,18749,18024,19859,19364,19363,19362,19361,19360,19359,19358,19357,19356,19355,19347,19239]

src_list_path = "/extern/UCNA/source_list_MB/"
src_peak_path = "/extern/UCNA/source_peaks_MB/"
outputFile = "/extern/UCNA/CalibrationPlots_MB/source_runs_RunPeriod_%i.dat"%(CalibrationPeriod)

runList = []

with open("../run_lists/Source_Calibration_Run_Period_%i.dat"%CalibrationPeriod) as runlist:
    for run in runlist:
        if os.path.isfile(src_list_path+"source_list_%i.dat"%int(run)):
            srcList = open(src_list_path+"source_list_%i.dat"%int(run))
            lines = []
            for line in srcList:
                lines.append(line)
            if int(lines[0])>0:
                runList.append(int(run))

print runList

outfile = open(outputFile,'w')

for run in runList:
    src_file = src_peak_path + "source_peaks_%i.dat"%run
    if MButils.fileExistsAndNotEmpty(src_file) and run not in omittedRuns:
        infile = open(src_file,'r')
        for line in infile:
            outfile.write(line)
        infile.close()

outfile.close()
    
