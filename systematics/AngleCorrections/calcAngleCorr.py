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
    delta0fracShift = 0.25 # 20% = 0.20
    delta1fracShift = 0.25
    delta2fracShift = 0.25
    delta3fracShift = 0.25

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

    #print(totalStatErr(A.realA,A.realAerr,220., 670.))
    #print(totalStatErr(type3A.realA,type3A.realAerr,220., 670.))
    #exit(0)

    # 
    
    # store fractions, delta_i values, and A_i values with errors
    totaldenom = [ ( ( 1./type0A.realAerr[i]**2 if type0A.realAerr[i]>0. else 0.)
                     +( 1./type1A.realAerr[i]**2 if type1A.realAerr[i]>0. else 0.)
                     +( 1./type2A.realAerr[i]**2 if type2A.realAerr[i]>0. else 0.)
                     +( 1./type3A.realAerr[i]**2 if type3A.realAerr[i]>0. else 0.) )
                   for i in range(0,len(type0A.realAerr)) ]
    
    #calculate the weighted average A
    weightedA = [ [( ((1./type0A.realAerr[i])**2*type0A.realA[i] if type0A.realAerr[i]>0. else 0.)+
                    ((1./type1A.realAerr[i])**2*type1A.realA[i] if type1A.realAerr[i]>0. else 0.)+
                    ((1./type2A.realAerr[i])**2*type2A.realA[i] if type2A.realAerr[i]>0. else 0.)+
                    ((1./type3A.realAerr[i])**2*type3A.realA[i] if type3A.realAerr[i]>0. else 0.))/(totaldenom[i] if totaldenom[i]!=0. else 1.)
                   for i in range(0,len(type0A.realAerr))],
                  [ 1./totaldenom[i]**2 if totaldenom[i]!=0. else 0. for i in range(0,len(type0A.realAerr))] ]
                  
                  


    frac0 = [ ((1./type0A.realAerr[i])**2/totaldenom[i] if type0A.realAerr[i]>0. else 0.) for i in range(0,len(type0A.realAerr)) ]
    delta0 = readAngleCorr(year,"D",delta0fracShift) 
    
    frac1 = [ ((1./type1A.realAerr[i])**2/totaldenom[i] if type1A.realAerr[i]>0. else 0.) for i in range(0,len(type1A.realAerr)) ]
    delta1 = readAngleCorr(year,"F",delta1fracShift) 
    
    frac2 = [ ((1./type2A.realAerr[i])**2/totaldenom[i] if type2A.realAerr[i]>0. else 0.) for i in range(0,len(type2A.realAerr)) ]
    delta2 = readAngleCorr(year,"J",delta2fracShift)
    
    frac3 = [ ((1./type3A.realAerr[i])**2/totaldenom[i] if type3A.realAerr[i]>0. else 0.) for i in range(0,len(type3A.realAerr)) ]
    delta3 = readAngleCorr(year,"K",delta3fracShift)


    # now create the total delta_3i corrections and their fractional uncertainties

    delta_30 = [ [(frac0[i]*delta0[0][i]*fabs(type0A.realA[i]/weightedA[0][i]) if weightedA[0][i]!=0. else 0.) for i in range(0,len(type0A.realA))],
                 [delta0fracShift for i in range(0,len(type0A.realA))] ]
    
    delta_31 = [ [(frac1[i]*delta1[0][i]*fabs(type1A.realA[i]/weightedA[0][i]) if weightedA[0][i]!=0. else 0.) for i in range(0,len(type1A.realA))],
                 [delta1fracShift for i in range(0,len(type1A.realA))] ]

    delta_32 = [ [(frac2[i]*delta2[0][i]*fabs(type2A.realA[i]/weightedA[0][i]) if weightedA[0][i]!=0. else 0.) for i in range(0,len(type2A.realA))],
                 [delta2fracShift for i in range(0,len(type2A.realA))] ]

    delta_33 = [ [(frac3[i]*delta3[0][i]*fabs(type3A.realA[i]/weightedA[0][i]) if weightedA[0][i]!=0. else 0.) for i in range(0,len(type3A.realA))],
                 [delta3fracShift for i in range(0,len(type3A.realA))] ]

    ##### including statistical err on correction
    #delta_30 = [ [(frac0[i]*delta0[0][i]*fabs(type0A.realA[i]/weightedA[0][i]) if weightedA[0][i]!=0. else 0.) for i in range(0,len(type0A.realA))],
    #             [sqrt(delta0fracShift**2+type0A.stat_percent_err[i]**2) for i in range(0,len(type0A.realA))] ]
    
    #delta_31 = [ [(frac1[i]*delta1[0][i]*fabs(type1A.realA[i]/weightedA[0][i]) if weightedA[0][i]!=0. else 0.) for i in range(0,len(type1A.realA))],
    #             [sqrt(delta1fracShift**2+type1A.stat_percent_err[i]**2) for i in range(0,len(type1A.realA))] ]
    #delta_32 = [ [(frac2[i]*delta2[0][i]*fabs(type2A.realA[i]/weightedA[0][i]) if weightedA[0][i]!=0. else 0.) for i in range(0,len(type2A.realA))],
    #             [sqrt(delta2fracShift**2+type2A.stat_percent_err[i]**2) for i in range(0,len(type2A.realA))] ]
    #delta_33 = [ [(frac3[i]*delta3[0][i]*fabs(type3A.realA[i]/weightedA[0][i]) if weightedA[0][i]!=0. else 0.) for i in range(0,len(type3A.realA))],
    #             [sqrt(delta3fracShift**2+type3A.stat_percent_err[i]**2) for i in range(0,len(type3A.realA))] ]

    # create the total delta_3 correction from the above corrections with the actual uncertainty

    delta_3 = [ [(1.+delta_30[0][i])*(1.+delta_31[0][i])*(1.+delta_32[0][i])*(1.+delta_33[0][i])-1. for i in range(0,len(delta_30[0]))],
                [((1.+delta_30[0][i])*(1.+delta_31[0][i])*(1.+delta_32[0][i])*(1.+delta_33[0][i]))*
                 sqrt( (delta_30[1][i]*delta_30[0][i]/(1.+delta_30[0][i]))**2 + (delta_31[1][i]*delta_31[0][i]/(1.+delta_31[0][i]))**2
                       + (delta_32[1][i]*delta_32[0][i]/(1.+delta_32[0][i]))**2 + (delta_33[1][i]*delta_33[0][i]/(1.+delta_33[0][i]))**2 )for i in range(0,len(delta_30[0]))] ]

    #for i in range(0,len(delta_3[0])):
    #    print("%0.0f %0.7f %0.7f"%(10.*i+5.,delta_3[0][i],(delta_3[1][i]/delta_3[0][i] if delta_3[0][i]!=0. else 0.)))


    statUncert = [fabs(weightedA[1][i]/weightedA[0][i]) if weightedA[0][i]!=0. else 0. for i in range(0,len(delta_3[0]))]

    totalFracUncertA = weightRealStats([delta_3[1][i]/(1.+delta_3[0][i]) for i in range(0,len(delta_3[0]))],statUncert,220., 670.)
    totalFracCorrA = weightRealStats([delta_3[0][i]/(1+delta_3[0][i]) for i in range(0,len(delta_3[0]))],statUncert,220., 670.)
    print
    print
    print("Total fractional Uncert: %f"%totalFracUncertA)
    print("Total fractional Corr: %f"%totalFracCorrA)
    print("WeightedA: %f"%weightRealStats(weightedA[0],weightedA[1],220.,670.))
    print("RealA: %f"%weightRealStats(A.realA,A.realAerr,220.,670.))


    # calculate the contributions to the uncertainty from each type
    total_delta_30 = weightRealStats([delta_30[0][i]/(1+delta_30[0][i]) for i in range(0,len(delta_30[0]))],statUncert,220.,670.)
    total_delta_31 = weightRealStats([delta_31[0][i]/(1+delta_31[0][i]) for i in range(0,len(delta_31[0]))],statUncert,220.,670.)
    total_delta_32 = weightRealStats([delta_32[0][i]/(1+delta_32[0][i]) for i in range(0,len(delta_32[0]))],statUncert,220.,670.)
    total_delta_33 = weightRealStats([delta_33[0][i]/(1+delta_33[0][i]) for i in range(0,len(delta_33[0]))],statUncert,220.,670.)
    
    print
    print("Idividual fractional Correction on A")
    print("Type 0: %f"%(total_delta_30/(1.+total_delta_30)))
    print("Type 1: %f"%(total_delta_31/(1.+total_delta_31)))
    print("Type 2: %f"%(total_delta_32/(1.+total_delta_32)))
    print("Type 3: %f"%(total_delta_33/(1.+total_delta_33)))

    print
    print("Idividual fractional uncertainties on A")
    print("Type 0: %f"%fabs(delta0fracShift*total_delta_30/(1.+total_delta_30)))
    print("Type 1: %f"%fabs(delta1fracShift*total_delta_31/(1.+total_delta_31)))
    print("Type 2: %f"%fabs(delta2fracShift*total_delta_32/(1.+total_delta_32)))
    print("Type 3: %f"%fabs(delta3fracShift*total_delta_33/(1.+total_delta_33)))
    print
    print("fractional Contribution to the total uncertainty as ratio of Delta_3i/Delta_30")
    print("Type 1: %f"%fabs(delta1fracShift*total_delta_31/(1.+total_delta_31)/(delta0fracShift*total_delta_30/(1.+total_delta_30))))
    print("Type 2: %f"%fabs(delta1fracShift*total_delta_32/(1.+total_delta_32)/(delta0fracShift*total_delta_30/(1.+total_delta_30))))
    print("Type 3: %f"%fabs(delta1fracShift*total_delta_33/(1.+total_delta_33)/(delta0fracShift*total_delta_30/(1.+total_delta_30))))
    

