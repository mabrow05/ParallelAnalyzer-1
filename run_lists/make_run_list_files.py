#!/usr/bin/python

#produce the files which hold the number of sources and what type for each run

import os

runType = "Source" #"Xenon"
year = "2012-2013"
numStart = None
logfile = None

if year[:4]=="2011":
    numStart = 1
    logfile = "../Aux/UCNA Run Log.txt"
elif year[:4]=="2012":
    numStart = 13
    logfile = "../Aux/UCNA Run Log 2012.txt"

runPeriod_runs = []
numPeriods=0
with open(logfile,"r") as log:
    srcPeriodContinue=False # Whether or not we are in a source period
    for line in log:
        if line[0]=="*" and line[7:9]==runType[0:2]:
            if not srcPeriodContinue:
                numPeriods = numPeriods+1
                runPeriod_runs.append([])
                srcPeriodContinue = True
            runPeriod_runs[numPeriods-1].append(int(line[1:6]))
            #print "%i %i"%(numPeriods-1+numStart,int(line[1:6]))

        elif line[0]=="*" and not line[7:9]==runType[0:2]:
            srcPeriodContinue = False

        elif line[0:10]=="##########":
            srcPeriodContinue = False
        
for i in range(0,numPeriods,1):
    with open("%s_Calibration_Run_Period_%i.dat"%(runType,i+numStart),"w") as output:
        for j in range(0,len(runPeriod_runs[i]),1):
            print "%i %i"%(i+numStart,runPeriod_runs[i][j])
            output.write("%i\n"%runPeriod_runs[i][j])


