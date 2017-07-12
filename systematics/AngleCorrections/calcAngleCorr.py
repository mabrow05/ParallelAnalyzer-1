#!/usr/bin/python
import sys
import csv
import os
from math import *

def getEnergyBin(energy):
    return int(energy)/10
    
def getBinEnergyMid(b):
    return b*10. + 5.

def getBinEnergyLowEdge(b):
    return b*10.

def getBinEnergyUpperEdge(b):
    return (b+1.)*10.

def weightRealStats(data,stats,e0,e1):
        """Weight the given data array by the real experimental statistics"""

        numer = sum([(1./stats[i])**2 * data[i] 
                     for i in range(getEnergyBin(e0),getEnergyBin(e1)) if e0<=getBinEnergyMid(i)<=e1])
        denom = sum([(1./stats[i])**2  
                     for i in range(getEnergyBin(e0),getEnergyBin(e1)) if e0<=getBinEnergyMid(i)<=e1])
        
        return numer/denom


def statUncertainties(year="2011-2012"):
    """Reads in the Statistical Uncertainty of every bin, calculates and 
        stores the total fitted asymmetry and the total percent error
        from statistics alone"""
    octLow=0
    octHigh=59
    if year!="2011-2012":
        octLow=60
        octHigh=121
    
    infile = open( os.environ["ANALYSIS_RESULTS"]+"/Asymmetries/"+
                   "UnCorr_OctetAsymmetries_AnaChC_Octets_%i-%i_BinByBin.txt"
                   %(octLow,octHigh),'r' )
    A_en = []
    A = []
    Aerr = []
    percErr = []
    
    if infile:
        for line in infile:
            l = line.split()
            if float(l[0])<1000.:
                A_en.append(float(l[0]))
                A.append(float(l[1]))
                Aerr.append(float(l[2]))
                if float(l[1])!=0.: 
                    percErr.append( fabs( float(l[2])/float(l[1]) ) )
                    
                else:
                    percErr.append( 1. )

    else: 
        print("Couldn't open file for statistical uncertainties")
        exit(0)


    return percErr
        

        

def readAngleCorr(year="2011-2012",percErr=0.2):
    A = [[],[]]
    with open('angleCorr_%s.txt'%(year),'rb') as tsvin:
        tsvin = csv.reader(tsvin, delimiter='\t')
        for row in tsvin:
            A[0].append(float(row[3]))
            A[1].append(fabs(float(row[3])*percErr))
            #print("%s\t%s\t%s"%(row[0],row[1],row[2]))

    return A



    
if __name__=="__main__":
    year = "2011-2012"
    emin = 230.
    emax = 750.
    perc_Err=0.25

    angleCorr = readAngleCorr(year,percErr=perc_Err)
    statErr = statUncertainties(year)
    print(weightRealStats(angleCorr[0],statErr,emin,emax))
    print(weightRealStats(angleCorr[1],statErr,emin,emax))
