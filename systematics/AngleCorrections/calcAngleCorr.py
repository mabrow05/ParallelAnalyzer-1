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
    delta0fracShift = 0. # 20% = 0.20
    delta1fracShift = 0.
    delta2fracShift = 0.
    delta3fracShift = 0.

    # create objects for each event type and total asymmetry
    A = uncertaintyHandler(year,"C",sim=True)
    type0A = uncertaintyHandler(year,"D",sim=True)
    type1A = uncertaintyHandler(year,"F",sim=True)
    type2A = uncertaintyHandler(year,"J",sim=True)
    type3A = uncertaintyHandler(year,"K",sim=True)

    # Load each of their statistical uncertainties and UnCorr A values
    A.statUncertainties()
    type0A.statUncertainties()
    type1A.statUncertainties()
    type2A.statUncertainties()
    type3A.statUncertainties()

    # 
    
    # store fractions, delta_i values, and A_i values with errors
    frac0 = [ A.stat_percent_err[i]**2/typ0A.stat_percent_err[i]**2 for i in range(0,len(A.stat_percent_err)) ]
    delta0 = readAngleCorr(year,"D",delta0fracShift) 
    frac1 = [ A.stat_percent_err[i]**2/typ1A.stat_percent_err[i]**2 for i in range(1,len(A.stat_percent_err)) ]
    delta1 = readAngleCorr(year,"F",delta1fracShift) 
    frac2 = [ A.stat_percent_err[i]**2/typ2A.stat_percent_err[i]**2 for i in range(2,len(A.stat_percent_err)) ]
    delta2 = readAngleCorr(year,"J",delta2fracShift) 
    frac3 = [ A.stat_percent_err[i]**2/typ3A.stat_percent_err[i]**2 for i in range(3,len(A.stat_percent_err)) ]
    delta3 = readAngleCorr(year,"K",,delta3fracShift)

    # now create the total delta_3i corrections

    delta_30 = [ [frac0[i]*delta0[0][i]*type0A.realA[i] for i in range(0,len(type0A.realA))],
                 [sqrt(delta0fracShift**2+type0A.stat_percent_err[i]**2) for i in range(0,len(type0A.realA))] ]
    delta_31 = [ [frac1[i]*delta1[0][i]*type1A.realA[i] for i in range(0,len(type1A.realA))],
                 [sqrt(delta1fracShift**2+type1A.stat_percent_err[i]**2) for i in range(0,len(type1A.realA))] ]
    delta_32 = [ [frac2[i]*delta2[0][i]*type2A.realA[i] for i in range(0,len(type2A.realA))],
                 [sqrt(delta2fracShift**2+type2A.stat_percent_err[i]**2) for i in range(0,len(type2A.realA))] ]
    delta_33 = [ [frac3[i]*delta3[0][i]*type3A.realA[i] for i in range(0,len(type3A.realA))],
                 [sqrt(delta3fracShift**2+type3A.stat_percent_err[i]**2) for i in range(0,len(type3A.realA))] ]

    # create the total delta_3 correction from the above corrections

    delta_3 = [ [(1.+delta_30[0][i])*(1.+delta_31[0][i])*(1.+delta_32[0][i])*(1.+delta_33[0][i]) for i in range(0,len(delta_30))],
                [sqrt( (delta_30[1][i]/(1.+delta_30[0][i]))**2 + (delta_31[1][i]/(1.+delta_31[0][i]))**2
                       + (delta_32[1][i]/(1.+delta_32[0][i]))**2 + (delta_33[1][i]/(1.+delta_33[0][i]))**2 )for i in range(0,len(delta_30))] ]

    for i in range(0,len(delta_3[0])):
        print("%f %f %f"%(10.*i+5.,delta_3[0][i],delta_3[1][i]))
