#!/usr/bin/python
import sys
import csv
import os
from math import *
sys.path.append("../../UNBLINDING")
from asymmetry import *
        
    
if __name__=="__main__":
    year = 2011
    emin = 220.
    emax = 670.
    delta0fracShift = 0.
    delta1fracShift = 0.
    delta2fracShift = 0.
    delta3fracShift = 0.

    # create objects for each event type and total asymmetry
    A = uncertaintyHandler(year,"C")
    type0A = uncertaintyHandler(year,"D")
    type1A = uncertaintyHandler(year,"F")
    type2A = uncertaintyHandler(year,"J")
    type3A = uncertaintyHandler(year,"K")

    # Load each of their statistical uncertainties and UnCorr A values
    A.statUncertainties()
    type0A.statUncertainties()
    type1A.statUncertainties()
    type2A.statUncertainties()
    type3A.statUncertainties()

    # 
    
    # store fractions, delta_i values, and A_i values with errors
    frac0 = [ A.stat_percent_err[i]**2/typ0A.stat_percent_err[i]**2 for i in range(0,len(A.stat_percent_err)) ]
    delta0 = readAngleCorr(year,"D") 
