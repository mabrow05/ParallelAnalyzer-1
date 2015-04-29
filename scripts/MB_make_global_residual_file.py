#!/usr/bin/python

### Script to generate files with all source data included to produce 
### error envelope type plot

import os
import MButils

periodLow=4
periodHigh=11

for side in ["East","West"]:
    outfile = open("residuals_global_%s_periods_%i-%i.dat"%(side,periodLow,periodHigh),"w")

    for period in range(periodLow,periodHigh+1,1):
        if os.path.isfile("residuals_%s_runPeriod_%i.dat"%(side,period)):
            resid = open("residuals_%s_runPeriod_%i.dat"%(side,period))
            for line in resid:
                outfile.write(line)

    outfile.close()
