#!/usr/bin/python

import os

runPeriods = [1,2,3,4,5,6,7,8,9,10,11,12]
InEnergy = False;
filename = ""
filein = None

for runPeriod in runPeriods:
    input_hold = []
    if InEnergy:
        filename = "../residuals/source_runs_EnergyPeaks_RunPeriod_%i_hold.dat"%runPeriod
        filein = open(filename, "r")
        for line in filein:
            input_hold.append(line.split())
        filein.close()
    else:
        filename = "../residuals/source_runs_RunPeriod_%i_hold.dat"%runPeriod
        filein = open(filename, "r")
        for line in filein:
            input_hold.append(line.split())
        filein.close()
    
    for line in input_hold:
        print '%s %s %s %s %s %s %s %s %s %s'%(line[0], line[1], line[2], line[3], line[4], line[5], line[6], line[7], line[8], line[9])

    fileout = open("../residuals/source_runs_RunPeriod_%i.dat"%runPeriod, "w")
    eQ = ""
    for line in input_hold:
        #print line[1]
        if float(line[1])==98.2:
            eQ="Ce"
        elif float(line[1])==331.2:
            eQ="Sn"#317.8
        elif float(line[1])==443.0:
            eQ="Bi2" #448.8
        elif float(line[1])==928.0:
            eQ="Bi1" #926.0
        fileout.write('%s %s %s %s %s %s %s %s %s %s\n'%(line[0], eQ, line[2], line[3], line[4], line[5], line[6], line[7], line[8], line[9]))

    fileout.close()
