#!/usr/bin/python

### Script to generate files with all source data included to produce 
### error envelope type plot

import os
import MButils

PMT=1 # 0 -> PMTs averaged over a side; 1,2,3,4 -> Individual PMT
periodLow=2
periodHigh=10
Periods = [1,2,3,4,5,6,7,8,9,10,11]
skipPeriods=[1,4,9]

for side in ["East","West"]:
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
