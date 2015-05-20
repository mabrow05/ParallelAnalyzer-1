#!/usr/bin/python

### Script to generate files with all source data included to produce 
### error envelope type plot

import os
import MButils

#runRange = [17359,19959]
CalibrationPeriods = [2]#,3,4,5,6,7,8,9,10,11] 
omittedRuns = [17381,17383,17385,17382,17886,17912,19232]
omittedRanges = [(19347,19364),(18020,18055)]

for Range in omittedRanges:
    for run in range(Range[0],Range[1]+1,1):
        omittedRuns.append(run)

print omittedRuns

src_list_path = os.getenv("SOURCE_LIST")
src_peak_path = os.getenv("SOURCE_PEAKS")

for CalibrationPeriod in CalibrationPeriods:
    outputFile = "../residuals/source_runs_RunPeriod_%i.dat"%(CalibrationPeriod)
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
    
