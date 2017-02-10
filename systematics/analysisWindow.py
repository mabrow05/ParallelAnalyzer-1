#!/usr/bin/python

import sys
import os
from math import *
#sys.path.append(".")
from EnergyUncertainty.EnergyErrors import *

# I need to read in statistical uncertainties to store in the 



def getEnergyBin(energy):
    return int(energy)/10
    
def getBinEnergyMid(b):
    return b*10. + 5.

def weightRealStats(data,stats,e0,e1):
        """Weight the given data array by the real experimental statistics"""
        
        numer = sum([(1./stats[i])**2 * data[i] 
                     for i in range(getEnergyBin(e0),getEnergyBin(e1)) if e0<=getBinEnergyMid(i)<=e1])
        denom = sum([(1./stats[i])**2  
                     for i in range(getEnergyBin(e0),getEnergyBin(e1)) if e0<=getBinEnergyMid(i)<=e1])

        return numer/denom

def totalStatErr(data,dataErr,e0,e1):
    numer = sum([(1./dataErr[i])**2 * data[i] 
                 for i in range(getEnergyBin(e0),getEnergyBin(e1)) if e0<=getBinEnergyMid(i)<=e1])
    weightsum = sum([(1./dataErr[i])**2 
                     for i in range(getEnergyBin(e0),getEnergyBin(e1)) if e0<=getBinEnergyMid(i)<=e1])
    A = numer/weightsum;
    return fabs( ( 1./sqrt(weightsum) ) / A )

def sumErrors(err):
    return sqrt(sum([er**2 for er in err]))


class analysisWindow:
    
    def __init__(self,year):
        self.stats = []
        self.A0 = 0.
        self.statUncert = 0.
        self.year = 0
        self.octLow = 0
        self.octHigh = 0
        self.realA = []
        self.realAerr = []
        self.stat_percent_err = []
        self.energy_err = []

        self.year = year
        print(self.year)
        if self.year == 2011:
            self.octLow = 0
            self.octHigh = 59
            
        elif self.year == 2012:
            self.octLow = 60
            self.octHigh = 121
        else:
            print('bad year given!!')
            exit(0)

    def statUncertainties(self):
        """Reads in the Statistical Uncertainty of every bin, calculates and 
        stores the total fitted asymmetry and the total percent error
        from statistics alone"""
        
        infile = open( os.environ["ANALYSIS_RESULTS"]+"/Asymmetries/"+
                       "UnCorr_OctetAsymmetries_AnaChA_180-780_Octets_%i-%i_BinByBin.txt"
                       %(self.octLow,self.octHigh),'r' )

        A = []
        Aerr = []
        percErr = []

        if infile:
            for line in infile:
                l = line.split()
                A.append(float(l[1]))
                Aerr.append(float(l[2]))
                if float(l[1])!=0.: 
                    percErr.append( fabs( float(l[2])/float(l[1]) ) )
                    
                else:
                    percErr.append( 1. )

        else: 
            print("Couldn't open file for statistical uncertainties")
            exit(0)
        
        self.realA = A
        self.realAerr = Aerr 
        self.stat_percent_err = percErr

    

    def readEnergyUncertainties(self):
        """Read in the Energy uncertainties and store the percent 
        error from every bin."""
        
        infile = open( os.environ["ANALYSIS_RESULTS"]+"/Corrections/"+
                       "EnergyLinearityUncertainty_%i.txt"%self.year)

        percErr = []

        if infile:
            for line in infile:
                if line[0]!="#":
                    l = line.split()
                    percErr.append( fabs(float(l[3])) )

        else: 
            print("Couldn't open file for Energy uncertainties")
            exit(0)
        

        self.energy_err = percErr

    def readSystematicUncertainties(self,errFac=0.2):
        """Read in the systematic correction bin-by-bin for final
        asymmetries, calculates delta_A/A, and then multiplies by the 
        error factor before filling the percent error for every bin."""
        
        
        percErr = []

    def makeTheoryUncertainties(self,lambda_err=0.):
        """Calculate the theory uncertainties given the error in lambda provided"""
        percErr = []
        
    

    def minimizer(self):
        """Minimizes the error and outputs all errors in order or largest to smallest"""
        
        self.statUncertainties()
        self.readEnergyUncertainties()
        
        for x in self.energy_err:
            print(x)        

        minErr = 100.
        enBinLow = None
        enBinHigh = None

        for lowBin in range(5,120-1):
            for highBin in range(lowBin+1,120):

                Errors = []

                Errors.append( weightRealStats(self.energy_err,self.stat_percent_err,220.,670.) )
                Errors.append( totalStatErr(self.realA,self.realAerr,220.,670.) )

                errSum = sumErrors(Errors)
                
                if errSum < minErr:
                    minErr = errSum
                    enBinLow = lowBin
                    enBinHigh = highBin


        #print("Weighted Energy Error: %f"%weightRealStats(self.energy_err,self.stat_percent_err,220.,670.))
        #print("Total Statistical Error: %f"%totalStatErr(self.realA,self.realAerr,220.,670.))
        print("Min total Err = %f"%minErr)
        print("Energy Range = %f-%f"%(getBinEnergyMid(enBinLow),getBinEnergyMid(enBinHigh)))

if __name__ == "__main__":
    
    minim = analysisWindow(2011)
    minim.minimizer()
