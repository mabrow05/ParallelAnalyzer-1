#!/usr/bin/python

import os

runs = [17238,17370,17521,17925,18361,18621,18749,19232,19359,19511,19857,19899]

sources = ["Ce139", "Sn113", "Bi207"]

runPeriod=1

for run in runs:
    outfile = open("weightedSimPeaks_runPeriod_%i.dat"%runPeriod,'w')
    for source in sources:
        #os.system("root -l -b -q 'weightPeaks.C(%i,\"%s\")'"%(run,source))

        # to combine the peaks from a run into a source run period file
        if 1:           
            filename = "%i_%s_weightedSimPeaks.dat"%(run,source)
            peakFile = open(filename,'r')
            lines = []
            for line in peakFile:
                lines.append(line)
                #print line
            peakFile.close()
            for i in range(len(lines),0,-1):
                outfile.write(lines[i-1]) # to print lower Bi first
                #print lines[i-1][:-1]
                if i>1: 
                    outfile.write("\n")
    runPeriod+=1
    outfile.close()
