#!/usr/bin/python

### Script to generate files with all source data included to produce 
### error envelope type plot

import os
import MButils

PMT=1 # 0 -> PMTs averaged over a side; 1,2,3,4 -> Individual PMT
Side="Both" #options: "East", "West", "Both"
periodLow=2
periodHigh=10
skipPeriods=[1,9]

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

    for period in range(periodLow,periodHigh+1,1):
        if period not in skipPeriods:
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
                    outfile.write(line)

    outfile.close()
