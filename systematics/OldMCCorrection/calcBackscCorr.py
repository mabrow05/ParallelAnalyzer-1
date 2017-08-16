#!/usr/bin/python
import sys
import csv
import os
from math import *
sys.path.append("../../UNBLINDING")
from asymmetry import *
        
delta0fracShift = 0.10 # 20% = 0.20
delta1fracShift = 0.30
delta2fracShift = 0.40
delta3fracShift = 0.20

useEffStatErr=True
effStatErr0 = 0.0025
effStatErr1 = 0.00074
effStatErr2 = 0.00062
effStatErr3 = 0.00059

incStatErr = True

def doBackscCorr(year,anaCh,emin,emax):
    
    Type0=False
    Type1=False
    Type2=False
    Type3=False
    if anaCh=="A" or anaCh=="B" or anaCh=="C" or anaCh=="D" or anaCh=="E":
        Type0 = True
    if anaCh=="A" or anaCh=="B" or anaCh=="C" or anaCh=="E" or anaCh=="F":
        Type1 = True
    if anaCh=="A" or anaCh=="C" or anaCh=="E" or anaCh=="G" or anaCh=="H" or anaCh=="J":
        Type2 = True
    if anaCh=="A" or anaCh=="C" or anaCh=="E" or anaCh=="G" or anaCh=="H" or anaCh=="K":
        Type3 = True

     # create objects for each event type and total asymmetry
    A = uncertaintyHandler(year,anaCh,sim=True)
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

    #for i in range(0,len(type2A.realA)):
    #    print("%f %f %f %f"%((i*10.+5.),type2A.realA[i],type2A.realAerr[i],type2A.stat_percent_err[i]))

    #exit(0)

    #print(totalStatErr(A.realA,A.realAerr,emin, emax))
    #print(totalStatErr(type3A.realA,type3A.realAerr,emin, emax))
    #exit(0)

    # 
    
    # store fractions, delta_i values, and A_i values with errors
    totaldenom = [ ( ( 1./type0A.realAerr[i]**2 if type0A.realAerr[i]>0. and Type0 else 0.)
                     +( 1./type1A.realAerr[i]**2 if type1A.realAerr[i]>0. and Type1 else 0.)
                     +( 1./type2A.realAerr[i]**2 if type2A.realAerr[i]>0. and Type2 else 0.)
                     +( 1./type3A.realAerr[i]**2 if type3A.realAerr[i]>0. and Type3 else 0.) )
                   for i in range(0,len(type0A.realAerr)) ]
    
    #calculate the weighted average A
    weightedA = [ [( ((1./type0A.realAerr[i])**2*type0A.realA[i] if type0A.realAerr[i]>0. and Type0 else 0.)+
                     ((1./type1A.realAerr[i])**2*type1A.realA[i] if type1A.realAerr[i]>0. and Type1 else 0.)+
                     ((1./type2A.realAerr[i])**2*type2A.realA[i] if type2A.realAerr[i]>0. and Type2 else 0.)+
                     ((1./type3A.realAerr[i])**2*type3A.realA[i] if type3A.realAerr[i]>0. and Type3 else 0.))/(totaldenom[i] if totaldenom[i]!=0. else 1.)
                   for i in range(0,len(type0A.realAerr))],
                  [ 1./totaldenom[i]**2 if totaldenom[i]!=0. else 0. for i in range(0,len(type0A.realAerr))] ]

                  


    frac0 = [ ((1./type0A.realAerr[i])**2/totaldenom[i] if type0A.realAerr[i]>0. and Type0 else 0.) for i in range(0,len(type0A.realAerr)) ]
    delta0 = readOldBackscCorr(year,anaCh,delta0fracShift,"0")     
    
    frac1 = [ ((1./type1A.realAerr[i])**2/totaldenom[i] if type1A.realAerr[i]>0. and Type1 else 0.) for i in range(0,len(type1A.realAerr)) ]
    delta1 = readOldBackscCorr(year,anaCh,delta1fracShift,"1") 
    
    frac2 = [ ((1./type2A.realAerr[i])**2/totaldenom[i] if type2A.realAerr[i]>0. and Type2 else 0.) for i in range(0,len(type2A.realAerr)) ]
    delta2 = readOldBackscCorr(year,anaCh,delta2fracShift,"2")
    
    frac3 = [ ((1./type3A.realAerr[i])**2/totaldenom[i] if type3A.realAerr[i]>0. and Type3 else 0.) for i in range(0,len(type3A.realAerr)) ]
    delta3 = readOldBackscCorr(year,anaCh,delta3fracShift,"3")
    #for i in range(0,len(frac3)):
    #    print("%f %f"%(delta3[0][i],frac3[i]))

    # now create the total delta_2i corrections and their fractional uncertainties

    #delta_20 = [ [(frac0[i]*delta0[0][i]*fabs(type0A.realA[i]/weightedA[0][i]) if weightedA[0][i]!=0. else 0.) for i in range(0,len(type0A.realA))],
    #             [sqrt(delta0fracShift**2+(type0A.stat_percent_err[i]**2 if incStatErr else 0.)) for i in range(0,len(type0A.realA))] ]
    
    #delta_21 = [ [(frac1[i]*delta1[0][i]*fabs(type1A.realA[i]/weightedA[0][i]) if weightedA[0][i]!=0. else 0.) for i in range(0,len(type1A.realA))],
    #             [sqrt(delta1fracShift**2+(type1A.stat_percent_err[i]**2 if incStatErr else 0.)) for i in range(0,len(type1A.realA))] ]

    #delta_22 = [ [(frac2[i]*delta2[0][i]*fabs(type2A.realA[i]/weightedA[0][i]) if weightedA[0][i]!=0. else 0.) for i in range(0,len(type2A.realA))],
     #            [sqrt(delta2fracShift**2+(type2A.stat_percent_err[i]**2 if incStatErr else 0.)) for i in range(0,len(type2A.realA))] ]

    #delta_23 = [ [(frac3[i]*delta3[0][i]*fabs(type3A.realA[i]/weightedA[0][i]) if weightedA[0][i]!=0. else 0.) for i in range(0,len(type3A.realA))],
    #             [sqrt(delta3fracShift**2+(type3A.stat_percent_err[i]**2 if incStatErr else 0.)) for i in range(0,len(type3A.realA))] ]

    ##### including statistical err on correction
    delta_20 = [ [delta0[0][i] if frac0[i]>0. and delta0[0][i]!=-1. else 0. for i in range(0,len(type0A.realA))],
                 [sqrt((delta0fracShift*delta0[0][i])**2
                       +((type0A.stat_percent_err[i]*delta0[0][i])**2 if incStatErr else 0.)
                       +(effStatErr0**2 if useEffStatErr else 0.)) for i in range(0,len(type0A.realA))] ]
    
    delta_21 = [ [delta1[0][i] if frac1[i]>0. and delta1[0][i]!=-1. else 0. for i in range(0,len(type1A.realA))],
                 [sqrt((delta1fracShift*delta1[0][i])**2
                       +((type1A.stat_percent_err[i]*delta1[0][i])**2 if incStatErr else 0.)
                       +(effStatErr1**2 if useEffStatErr else 0.)) for i in range(0,len(type1A.realA))] ]

    delta_22 = [ [delta2[0][i] if frac2[i]>0. and delta2[0][i]!=-1. else 0. for i in range(0,len(type2A.realA))],
                 [sqrt((delta2fracShift*delta2[0][i])**2
                       +((type2A.stat_percent_err[i]*delta2[0][i])**2 if incStatErr else 0.)
                       +(effStatErr2**2 if useEffStatErr else 0.)) for i in range(0,len(type2A.realA))] ]

    delta_23 = [ [delta3[0][i] if frac3[i]>0. and delta3[0][i]!=-1. else 0. for i in range(0,len(type3A.realA))],
                 [sqrt((delta3fracShift*delta3[0][i])**2
                       +((type3A.stat_percent_err[i]*delta2[0][i])**2 if incStatErr else 0.)
                       +(effStatErr3**2 if useEffStatErr else 0.)) for i in range(0,len(type3A.realA))] ]

    # create the total delta_2 correction from the above corrections with the actual uncertainty

    delta_2 = [ [(1.+delta_20[0][i])*(1.+delta_21[0][i])*(1.+delta_22[0][i])*(1.+delta_23[0][i])-1. for i in range(0,len(delta_20[0]))],
                [((1.+delta_20[0][i])*(1.+delta_21[0][i])*(1.+delta_22[0][i])*(1.+delta_23[0][i]))*
                 sqrt( (delta_20[1][i]/(1.+delta_20[0][i]))**2 + (delta_21[1][i]/(1.+delta_21[0][i]))**2
                       + (delta_22[1][i]/(1.+delta_22[0][i]))**2 + (delta_23[1][i]/(1.+delta_23[0][i]))**2 )for i in range(0,len(delta_20[0]))] ]

    #for i in range(0,len(delta_2[0])):
    #    print("%0.0f %0.7f %0.7f"%(10.*i+5.,delta_2[0][i],(delta_2[1][i]/delta_2[0][i] if delta_2[0][i]!=0. else 0.)))


    statUncert = [fabs(weightedA[1][i]/weightedA[0][i]) if weightedA[0][i]!=0. else 0. for i in range(0,len(delta_2[0]))]

    totalFracUncertA = weightRealStats([delta_2[1][i]/(1.+delta_2[0][i]) for i in range(0,len(delta_2[0]))],statUncert,emin, emax)
    totalFracCorrA = weightRealStats([delta_2[0][i]/(1+delta_2[0][i]) for i in range(0,len(delta_2[0]))],statUncert,emin, emax)
    print
    print
    print("Total fractional Uncert: %f"%totalFracUncertA)
    print("Total fractional Corr: %f"%totalFracCorrA)
    print("WeightedA: %f"%weightRealStats(weightedA[0],weightedA[1],emin,emax))
    print("RealA: %f"%weightRealStats(A.realA,A.realAerr,emin,emax))


    # calculate the contributions to the uncertainty from each type
    total_delta_20 = weightRealStats([delta_20[0][i]/(1+delta_20[0][i]) for i in range(0,len(delta_20[0]))],statUncert,emin,emax)
    total_delta_21 = weightRealStats([delta_21[0][i]/(1+delta_21[0][i]) for i in range(0,len(delta_21[0]))],statUncert,emin,emax)
    total_delta_22 = weightRealStats([delta_22[0][i]/(1+delta_22[0][i]) for i in range(0,len(delta_22[0]))],statUncert,emin,emax)
    total_delta_23 = weightRealStats([delta_23[0][i]/(1+delta_23[0][i]) for i in range(0,len(delta_23[0]))],statUncert,emin,emax)
    
    total_delta_20err = weightRealStats([delta_20[1][i]/(1+delta_20[0][i]) for i in range(0,len(delta_20[0]))],statUncert,emin,emax)
    total_delta_21err = weightRealStats([delta_21[1][i]/(1+delta_21[0][i]) for i in range(0,len(delta_21[0]))],statUncert,emin,emax)
    total_delta_22err = weightRealStats([delta_22[1][i]/(1+delta_22[0][i]) for i in range(0,len(delta_22[0]))],statUncert,emin,emax)
    total_delta_23err = weightRealStats([delta_23[1][i]/(1+delta_23[0][i]) for i in range(0,len(delta_23[0]))],statUncert,emin,emax)
    

    print
    print("Idividual fractional Correction on A")
    print("Type 0: %f"%(total_delta_20))
    print("Type 1: %f"%(total_delta_21))
    print("Type 2: %f"%(total_delta_22))
    print("Type 3: %f"%(total_delta_23))

    print
    print("Idividual fractional uncertainties on A")
    print("Type 0: %f"%fabs(total_delta_20err))
    print("Type 1: %f"%fabs(total_delta_21err))
    print("Type 2: %f"%fabs(total_delta_22err))
    print("Type 3: %f"%fabs(total_delta_23err))
    print
    #print("fractional Contribution to the total uncertainty as ratio of delta_2i/delta_20")
    #print("Type 1: %f"%fabs(delta1fracShift*total_delta_21/(1.+total_delta_21)/(delta0fracShift*total_delta_20/(1.+total_delta_20))))
    #print("Type 2: %f"%fabs(delta1fracShift*total_delta_22/(1.+total_delta_22)/(delta0fracShift*total_delta_20/(1.+total_delta_20))))
    #print("Type 3: %f"%fabs(delta1fracShift*total_delta_23/(1.+total_delta_23)/(delta0fracShift*total_delta_20/(1.+total_delta_20))))

    with open("%s_delta_20_anaCh%s.txt"%("2011-2012" if year==2011 else "2012-2013",anaCh),"w") as f:
        for i in range(0,len(delta_20[0])):
            f.write("%f\t%0.7f\t%0.7f\n"%((i*10.+5.),delta_20[0][i],fabs(delta_20[1][i])))

    with open("%s_delta_21_anaCh%s.txt"%("2011-2012" if year==2011 else "2012-2013",anaCh),"w") as f:
        for i in range(0,len(delta_21[0])):
            f.write("%f\t%0.7f\t%0.7f\n"%((i*10.+5.),delta_21[0][i],fabs(delta_21[1][i])))

    with open("%s_delta_22_anaCh%s.txt"%("2011-2012" if year==2011 else "2012-2013",anaCh),"w") as f:
        for i in range(0,len(delta_22[0])):
            f.write("%f\t%0.7f\t%0.7f\n"%((i*10.+5.),delta_22[0][i],fabs(delta_22[1][i])))
            
    with open("%s_delta_23_anaCh%s.txt"%("2011-2012" if year==2011 else "2012-2013",anaCh),"w") as f:
        for i in range(0,len(delta_23[0])):
            f.write("%f\t%0.7f\t%0.7f\n"%((i*10.+5.),delta_23[0][i],fabs(delta_23[1][i])))

            
    with open("%s_delta_2_anaCh%s.txt"%("2011-2012" if year==2011 else "2012-2013",anaCh),"w") as f:
        for i in range(0,len(delta_2[0])):
            f.write("%f\t%0.7f\t%0.7f\n"%((i*10.+5.),delta_2[0][i],fabs(delta_2[1][i])))



if __name__=="__main__":
    
    Year = [2011,2012]
    anaChoices = ["A","B","C","D","F","G","H","J","K"]
    
    emin = 190.
    emax = 750.
    
    #doBackscCorr(2012,"C",emin,emax)
    for year in Year:
        for anaCh in anaChoices:
            doBackscCorr(year,anaCh,emin,emax)
    
    

