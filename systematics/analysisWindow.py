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

    if len(data)!=len(dataErr):
        return 0

    numer = sum([(1./dataErr[i])**2 * data[i] 
                 for i in range(getEnergyBin(e0),getEnergyBin(e1)) if e0<=getBinEnergyMid(i)<=e1])
    weightsum = sum([(1./dataErr[i])**2 
                     for i in range(getEnergyBin(e0),getEnergyBin(e1)) if e0<=getBinEnergyMid(i)<=e1])
    A = numer/weightsum if weightsum>0. else 1.
    return fabs( ( 1./sqrt(weightsum) ) / A ) if weightsum>0. and fabs(A)>0. else 1.

def sumErrors(err):
    return sqrt(sum([er**2 for er in err]))


class uncertaintyHandler:
    
    def __init__(self,year,anaChoice="C"):
        self.anaChoice = anaChoice
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
        self.theory_err = []
        self.theory_corr = []
        self.systA = []
        self.syst_err = []
        self.syst_corr = []

        self.year = year
        
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
                       "UnCorr_OctetAsymmetries_AnaCh%s_Octets_%i-%i_BinByBin.txt"
                       %(self.anaChoice,self.octLow,self.octHigh),'r' )
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
        
        self.realA = A
        self.realAerr = Aerr 
        self.stat_percent_err = percErr


    def readMCSystematicCorrections(self):
        """Read in the systematic correction bin-by-bin for final
        asymmetries, calculates delta_A/A"""

        A_corr = []

        # begin with making sure uncorrected asymm is loaded
        if len(self.stat_percent_err)==0: 
            self.statUncertainties() 

         
        # Now read in corrected asymmetry
        infile = open( os.environ["ANALYSIS_RESULTS"]+"/Asymmetries/"+
                       "DeltaExpOnly_OctetAsymmetries_AnaCh%s_Octets_%i-%i_BinByBin.txt"
                       %(self.anaChoice,self.octLow,self.octHigh),'r' )

        if infile:
            for line in infile:
                l = line.split()
                if float(l[0])<1000.:
                    A_corr.append(float(l[1]))
                    
        else: 
            print("Couldn't open file for statistical uncertainties")
            exit(0)

        infile.close()

        
        # Calculate the percent error on the correction
        corr = []
        percErr = []
        
        for i in range( 0, len(A_corr) ):
            correction = 100.
            if fabs(self.realA[i])>0.:
                correction = A_corr[i]/self.realA[i] - 1. 
            
            corr.append(correction)
            percErr.append( fabs(correction*0.2) )
            
        self.systA = A_corr
        self.syst_err = percErr
        self.syst_corr = corr

    def calcMCSystematicError(self,elow,ehigh):

        if len(self.systA)==0:
            self.readMCSystematicCorrections()

         # Now Simulation
        infile = open( os.environ["SIM_ANALYSIS_RESULTS"]+"/Asymmetries/"+
                       "DeltaExpOnly_OctetAsymmetries_AnaCh%s_Octets_%i-%i_BinByBin.txt"
                       %(self.anaChoice,self.octLow,self.octHigh),'r' )
        sim_A_en = []
        sim_A = []
        sim_Aerr = []
        sim_percErr = []

        if infile:
            for line in infile:
                l = line.split()
                if float(l[0])<1000.:
                    sim_A_en.append(float(l[0]))
                    sim_A.append(float(l[1]))
                    sim_Aerr.append(float(l[2]))
                    if float(l[1])!=0.: 
                        sim_percErr.append( fabs( float(l[2])/float(l[1]) ) )
                        
                    else:
                        sim_percErr.append( 1. )
                    

        else: 
            print("Couldn't open file for statistical uncertainties")
            exit(0)

    ### here is the crucial part. We take the anaChoice D sim and data to calculate 
        ### delta A diff to subtract off of the difference between sim and data from 
        ### the analysis choice of interest.

        # Data first
        anaD = uncertaintyHandler(self.year,"D")
        anaD.readMCSystematicCorrections()

        # Now Simulation
        infile = open( os.environ["SIM_ANALYSIS_RESULTS"]+"/Asymmetries/"+
                       "DeltaExpOnly_OctetAsymmetries_AnaCh%s_Octets_%i-%i_BinByBin.txt"
                       %(self.anaChoice,self.octLow,self.octHigh),'r' )
        simD_A_en = []
        simD_A = []
        simD_Aerr = []
        simD_percErr = []

        if infile:
            for line in infile:
                l = line.split()
                if float(l[0])<1000.:
                    simD_A_en.append(float(l[0]))
                    simD_A.append(float(l[1]))
                    simD_Aerr.append(float(l[2]))
                    if float(l[1])!=0.: 
                        simD_percErr.append( fabs( float(l[2])/float(l[1]) ) )
                        #print(fabs( float(l[2])/float(l[1]) ) , l[1] )
                    else:
                        simD_percErr.append( 1. )

        else: 
            print("Couldn't open file for statistical uncertainties")
            exit(0)
        
       
        ### Now calculate asymmetry over bins
        dataA_anaD = weightRealStats(anaD.systA,anaD.stat_percent_err,elow,ehigh)
        simA_anaD = weightRealStats(simD_A,simD_percErr,elow,ehigh)

        dataA = weightRealStats(self.systA,self.stat_percent_err,elow,ehigh)
        simA = weightRealStats(sim_A,sim_percErr,elow,ehigh)

        deltaAOffset = dataA_anaD-simA_anaD # the reference delta

        deltaAOffsetPerc = deltaAOffset / dataA_anaD if dataA_anaD!=0 else 1000.
        #return fabs(deltaAOffsetPerc)
        

       # print("( (%f-%f)/%f ) - %f = %f"%(dataA,simA,dataA,deltaAOffsetPerc,fabs( ( (dataA-simA) / dataA ) - deltaAOffsetPerc ) if fabs(dataA)>0. else 100000.))
        #return fabs( ( (dataA-simA) / dataA ) - deltaAOffsetPerc ) if fabs(dataA)>0. else 100000.

        print("( ( %f-%f ) - %f ) / %f = %f"%(dataA,simA,deltaAOffset,dataA,fabs( ( (dataA-simA)  - deltaAOffset) / dataA ) if dataA!=0. else 100000. ))
        return fabs( ( (dataA-simA) -deltaAOffset) / dataA ) if fabs(dataA)>0. else 100000. 
        
    

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

    def makeSystematicCorrections(self,errFac=0.2):
        """Read in the systematic correction bin-by-bin for final
        asymmetries, calculates delta_A/A, and then multiplies by the 
        error factor before filling the percent error for every bin."""

        A_uncorr = []
        A_corr = []

        # begin with uncorrected asymmetry
        infile = open( os.environ["ANALYSIS_RESULTS"]+"/Asymmetries/"+
                       "UnCorr_OctetAsymmetries_AnaCh%s_Octets_%i-%i_BinByBin.txt"
                       %(self.anaChoice,self.octLow,self.octHigh),'r' )

        if infile:
            for line in infile:
                l = line.split()
                if float(l[0])<1000.:
                    A_uncorr.append(float(l[1]))

        else: 
            print("Couldn't open file for statistical uncertainties")
            exit(0)

         

        infile.close()

        # Now read in corrected asymmetry
        infile = open( os.environ["ANALYSIS_RESULTS"]+"/Asymmetries/"+
                       "DeltaExpOnly_OctetAsymmetries_AnaCh%s_Octets_%i-%i_BinByBin.txt"
                       %(self.anaChoice,self.octLow,self.octHigh),'r' )

        if infile:
            for line in infile:
                l = line.split()
                if float(l[0])<1000.:
                    A_corr.append(float(l[1]))

        else: 
            print("Couldn't open file for statistical uncertainties")
            exit(0)

        infile.close()

        #for x in A_corr:
        #    print(x)

        # Calculate the percent error on the correction
        corr = []
        percErr = []
        
        for i in range( 0, len(A_corr) ):
            correction = 100.
            if fabs(A_uncorr[i])>0.:
                correction = A_corr[i]/A_uncorr[i] - 1. 
            
            corr.append(correction)
            percErr.append( fabs(correction*errFac) )

        self.syst_err = percErr
        self.syst_corr = corr


    


    def makeTheoryUncertainties(self,lambda_err=0.):
        """Calculate the theory uncertainties given the error in lambda provided"""
        percErr = []
        
    def calcEnergyUncert(self,emin,emax):

        if len(self.energy_err)==0: 
            self.readEnergyUncertainties() 
            
        if len(self.stat_percent_err)==0: 
            self.statUncertainties() 
            
        return weightRealStats(self.energy_err,self.stat_percent_err,emin,emax)

    def calcStatUncert(self,emin,emax):

        if len(self.stat_percent_err)==0: 
            self.statUncertainties() 
            
        return totalStatErr(self.realA,self.realAerr,emin,emax) 
        

    def calcSystematicUncert(self,emin,emax):

        if len(self.syst_err)==0: 
            self.makeSystematicCorrections() 
            
        if len(self.stat_percent_err)==0: 
            self.statUncertainties() 
            
        return weightRealStats(self.syst_err,self.stat_percent_err,emin,emax)


    def minimizer(self):
        """Minimizes the error and outputs all errors in order or largest to smallest"""
        
        self.statUncertainties()
        self.readEnergyUncertainties()
        #self.makeSystematicCorrections()
        self.readMCSystematicCorrections()
        

        minErr = 100.
        minEnergyErr = 0.
        minStatisticalErr = 0.
        minSystematicErr = 0.
        enBinLow = None
        enBinHigh = None

        for lowBin in range(0,79-1):
            for highBin in range(lowBin+1,79):

                Errors = []

                Errors.append( self.calcEnergyUncert(     getBinEnergyMid(lowBin),getBinEnergyMid(highBin) ) )
                Errors.append( self.calcStatUncert(       getBinEnergyMid(lowBin),getBinEnergyMid(highBin) ) )

                ##### Systematic error handling

                # calcSystematicUncert applies a uniform percent uncertainty to every bin correction
                #Errors.append( self.calcSystematicUncert( getBinEnergyMid(lowBin),getBinEnergyMid(highBin) ) )
                
                # Ignores systematics in minimization
                #Errors.append(0.)

                # Multiplies final correction by some fractional percent uncertainty ( 20% in 2010 )
                #Errors.append( fabs( weightRealStats(self.syst_corr,self.stat_percent_err,getBinEnergyMid(lowBin),getBinEnergyMid(highBin)) ) * 0.2)

                Errors.append( self.calcMCSystematicError(getBinEnergyMid(lowBin),getBinEnergyMid(highBin)) );

                errSum = sumErrors(Errors)
                
                if errSum < minErr:
                    minErr = errSum
                    enBinLow = lowBin
                    enBinHigh = highBin
                    minEnergyErr = Errors[0]
                    minStatisticalErr = Errors[1]
                    minSystematicErr = Errors[2]


        minSystCorrTotal = weightRealStats(self.syst_corr,self.stat_percent_err,getBinEnergyMid(enBinLow),getBinEnergyMid(enBinHigh))
        print
        print("Systematic Correction: %f"%minSystCorrTotal)
        print
        print("Weighted Energy Error: %f"%minEnergyErr)
        print("Total Statistical Error: %f"%minStatisticalErr)
        print("Total Systematic Error: %f"%minSystematicErr)
        print("Min total Err = %f"%minErr)
        print("Energy Range = %f-%f"%(getBinEnergyMid(enBinLow),getBinEnergyMid(enBinHigh)))

if __name__ == "__main__":
    
    uncert = uncertaintyHandler(2011,"C")
    uncert.minimizer()
    print(uncert.calcEnergyUncert(220,680))
    print(uncert.calcStatUncert(220,680))
    print(uncert.calcSystematicUncert(220,680))
    print(uncert.calcMCSystematicError(220,680) )
