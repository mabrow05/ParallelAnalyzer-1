#!/usr/bin/python
import sys
import csv
import os
from math import *

def readAsymm(year,eastEff,westEff):
    A = []
    with open('%s_asymmErecon_E%0.2f_W%0.2f.txt'%(year,eastEff,westEff),'rb') as tsvin:
    #with open('%s_BHasymm_%sField_polW.txt'%(year,field),'rb') as tsvin:
        tsvin = csv.reader(tsvin, delimiter='\t')
        for row in tsvin:
            A.append([float(row[0]),float(row[1]),float(row[2])])
            #print("%s\t%s\t%s"%(row[0],row[1],row[2]))

    return A


def weightAsymm(data,emin,emax):
    numer = sum([ x[1]/x[2]**2 for x in data if emin<x[0]<emax ])
    denom = sum([ 1./x[2]**2 for x in data if emin<x[0]<emax ])
    return [numer/denom,sqrt(1./denom)]

def writeCorrectionByBin(uncorr,corr,year):
    with open("efficiencyCorr_%s.txt"%year,'w') as fout:
        for i in range(0,len(uncorr)):
            c = corr[i][1]/uncorr[i][1] if fabs(uncorr[i][1])>0. else 1.
            fout.write("%f\t%0.10f\n"%(uncorr[i][0],c))
            #print("%f\t%0.10f"%(uncorr[i][0],c))
                



    
if __name__=="__main__":
    year = "2011-2012"
    emin = 230.
    emax = 750.
    filemin = 0
    filemax = 200
    
    A_eff = readAsymm(year,1.19,1.09)
    A_eff_int = weightAsymm(A_eff,emin,emax)

    A = readAsymm(year,0.,0.)
    A_int = weightAsymm(A,emin,emax)

    writeCorrectionByBin(A_eff,A,year)
    
    print("A = %f +/- %f"%(A_int[0],A_int[1]))
    print("A_eff = %f +/- %f"%(A_eff_int[0],A_eff_int[1]))
    print("DeltaFieldDip (DeltaA/A) = %f"%(A_int[0]/A_eff_int[0] - 1))
    print("UnCorrelated err = %f"%(fabs(A_int[0]/A_eff_int[0])*
                      sqrt( (A_int[1]/A_int[0])**2 +(A_eff_int[1]/A_eff_int[0])**2 ) ) )
    
