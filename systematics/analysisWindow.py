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

def getAsymmetry(anaCh,octLow,octHigh,elow,ehigh,sim=True):
    os.system("./MC_Corrections/AsymmCalculator.exe %s %i %i %f %f %i false > log.txt"%(anaCh,octLow, octHigh, elow,ehigh,int(sim)))
    infile = open("asymm_hold.txt",'r')
    asymm=[0,0]
    if infile:
        for line in infile:
            entry = line.split()
            asymm[0] = float(entry[1])
            asymm[1] = float(entry[2])
                #print("type0sim %f"%type0sim)
        infile.close()
    return asymm

def readAngleCorr(year=2011,anaCh="C",percErr=0.2):
    A = [[],[]]
    yearString=None
    if year==2011:
        yearString = "2011-2012"
    else:
        yearString = "2012-2013"
    with open(os.environ["ANALYSIS_CODE"]+'systematics/AngleCorrections/%s_DeltaAngle_anaCh%s.txt'%(yearString,anaCh),'rb') as tsvin:
        tsvin = csv.reader(tsvin, delimiter='\t')
        for row in tsvin:
            A[0].append(float(row[2]))
            A[1].append(fabs(float(row[2])*percErr))
            #print("%s\t%s\t%s"%(row[0],row[1],row[2]))

    return A

def readBackscCorr(year=2011,anaCh="C",percErr=0.2):
    A = [[],[]]
    yearString=None
    if year==2011:
        yearString = "2011-2012"
    else:
        yearString = "2012-2013"
        
    with open(os.environ["ANALYSIS_CODE"]+"/systematics/OldMCCorrection/"+
              "deltaBSALL_anaCh%s_%s.txt"%(anaCh,yearString),'rb') as tsvin:
        tsvin = csv.reader(tsvin, delimiter='\t')
        for row in tsvin:
            A[0].append(float(row[1]))
            A[1].append(fabs(float(row[1])*percErr))
            #print("%s\t%s\t%s"%(row[0],row[1],row[2]))

    return A


class uncertaintyHandler:
    
    def __init__(self,year,anaChoice="C"):
        self.geom = None
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
            self.geom = "2011-2012"
            self.octLow = 0
            self.octHigh = 59
            
        elif self.year == 2012:
            self.geom = "2012-2013"
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
            print(correction)
            corr.append(correction)
            percErr.append( fabs(correction*0.2) )
            
        self.systA = A_corr
        self.syst_err = percErr
        self.syst_corr = corr

    


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

    def makeSystematicCorrections(self,errFacDeltaBS=0.2,errFacDeltaAngle=0.25):
        """Read in the systematic correction bin-by-bin for final
        asymmetries, calculates delta_A/A, and then multiplies by the 
        error factor before filling the percent error for every bin."""

        BS_corr = readBackscCorr(self.year,self.anaChoice,errFacDeltaBS)
        Angle_corr = readAngleCorr(self.year,self.anaChoice,errFacDeltaAngle)

        MC_corr = []
        MC_corrErr = []
        
        for i in range( 0, len(BS_corr[0]) ):
       
            MC_corr.append( (1+BS_corr[0][i])*(1+Angle_corr[0][i])-1. )
            MC_corrErr.append( sqrt( ((1+BS_corr[0][i])*Angle_corr[1][i])**2 + ((1+Angle_corr[0][i])*BS_corr[1][i])**2) )
            #if fabs(A_uncorr[i])>0.:
            #     correction = A_BScorr[i]/A_uncorr[i] - 1. 
            
        self.syst_err = MC_corrErr
        self.syst_corr = MC_corr


    


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

        #if len(self.stat_percent_err)==0: 
        #    self.statUncertainties() 
            
        #return totalStatErr(self.realA,self.realAerr,emin,emax) 
        asymm = getAsymmetry(self.anaChoice,self.octLow,self.octHigh,emin,
                             emax,False)
        self.A0 = asymm[0]
        return fabs( asymm[1]/asymm[0] )

    def calcSystematicUncert(self,emin,emax):

        if len(self.syst_err)==0: 
            self.makeSystematicCorrections() 
            
        if len(self.stat_percent_err)==0: 
            self.statUncertainties() 
            
        return weightRealStats(self.syst_err,self.stat_percent_err,emin,emax)

    def gainUncert(self,emin,emax):

        if len(self.stat_percent_err)==0: 
            self.statUncertainties() 

        infile = open( os.environ["ANALYSIS_RESULTS"]+"/Corrections/"+
                       "gainUncertainty_%i.txt"%self.year)

        gainErr = []

        if infile:
            for line in infile:
                if line[0]!="#":
                    l = line.split()
                    gainErr.append( fabs(float(l[3])) )
                    #print(float(l[3]))

        else: 
            print("Couldn't open file for gain uncertainties")
            exit(0)
        
            
        return weightRealStats(gainErr,self.stat_percent_err,emin,emax)



    def minimizer(self,errFacDeltaBS=0.2,errFacDeltaAngle=0.2):
        """Minimizes the error and outputs all errors in order or largest to smallest"""

        self.statUncertainties()
        self.readEnergyUncertainties()
        self.makeSystematicCorrections(errFacDeltaBS,errFacDeltaAngle)
        #self.readMCSystematicCorrections()

        ofile = open("%i_Uncertainties.txt"%self.year, 'w')
        

        minErr = 100.
        minEnergyErr = 0.
        minStatisticalErr = 0.
        minSystematicErr = 0.
        enBinLow = None
        enBinHigh = None

        minNumBins = 30
        lowbound = 10
        upperbound = 75
        for lowBin in range(lowbound,upperbound+1-minNumBins):
            for highBin in range(lowBin+minNumBins,upperbound+1):

                Errors = []

                Errors.append( self.calcEnergyUncert(     getBinEnergyMid(lowBin),getBinEnergyMid(highBin) ) )
                Errors.append( self.calcStatUncert(       getBinEnergyLowEdge(lowBin),getBinEnergyUpperEdge(highBin) ) )
                Errors.append( self.calcSystematicUncert( getBinEnergyMid(lowBin),getBinEnergyMid(highBin) ) )

                errSum = sumErrors(Errors)
                


                print("****************************")
                print("%i-%i keV total Error = %f\n"%(getBinEnergyLowEdge(lowBin),getBinEnergyUpperEdge(highBin),errSum))
                print("Asymm = %f\n"%self.A0)
                print("Weighted Energy Error: %f"%Errors[0])
                print("Total Statistical Error: %f"%Errors[1])
                print("Total Systematic Error: %f"%Errors[2])

                
                ofile.write("****************************\n")
                ofile.write("Asymm= %f\n"%self.A0)
                ofile.write("%i-%i keV_total_Error= %f\n\n"%(getBinEnergyLowEdge(lowBin),getBinEnergyUpperEdge(highBin),errSum))
                ofile.write("Weighted_Energy_Error: %f\n"%Errors[0])
                ofile.write("Total_Statistical_Error: %f\n"%Errors[1])
                ofile.write("Total_Systematic_Error: %f\n"%Errors[2])

                
                if errSum < minErr:
                    minErr = errSum
                    enBinLow = lowBin
                    enBinHigh = highBin
                    minEnergyErr = Errors[0]
                    minStatisticalErr = Errors[1]
                    minSystematicErr = Errors[2]


        ofile.close()

       # minSystCorrTotal = weightRealStats(self.syst_corr,self.stat_percent_err,getBinEnergyMid(enBinLow),getBinEnergyMid(enBinHigh))
        print
        #print("Systematic Correction: %f"%minSystCorrTotal)
        print
        print("Weighted Energy Error: %f"%minEnergyErr)
        print("Total Statistical Error: %f"%minStatisticalErr)
        print("Total Systematic Error: %f"%minSystematicErr)
        print("Min total Err = %f"%minErr)
        print("Energy Range = %f-%f"%(getBinEnergyLowEdge(enBinLow),getBinEnergyUpperEdge(enBinHigh)))

if __name__ == "__main__":
    
    year=2011
    uncert = uncertaintyHandler(year,"C")
    #uncert.minimizer(errFacDeltaBS=0.25,errFacDeltaAngle=0.25)

    #2011: Energy Range = 210.-760. for bsErr=0.25 and angleErr=0.25 : 0.0074
    #2012: Energy Range = 210.-760. for bsErr=0.25 and angleErr=0.25 : 0.0084

    #2011: Energy Range = 210.-750. for bsErr=0.2 and angleErr=0.25 : 0.006945
    #2012: Energy Range = 200.-740. for bsErr=0.2 and angleErr=0.25 : 0.00811
    # BUT!!! 2012 gives same uncertainty of 2011 energy range

    #2011: Energy Range = 230.-740. for bsErr=0.2 and angleErr=0.2
    #2012: Energy Range = 190.-760. for bsErr=0.2 and angleErr=0.2

    #also consider analysis choice "B" which seems to have the lowers Total uncert

    if 1:
    
        lowBin = 21
        highBin = 74
        
        uncert.statUncertainties()
        uncert.readEnergyUncertainties()
        uncert.makeSystematicCorrections(errFacDeltaBS=0.20,errFacDeltaAngle=0.25)
        
        Errors = []
        
        Errors.append( uncert.calcEnergyUncert(     getBinEnergyMid(lowBin),getBinEnergyMid(highBin) ) )
        Errors.append( uncert.calcStatUncert(       getBinEnergyLowEdge(lowBin),getBinEnergyUpperEdge(highBin) ) )
        Errors.append( uncert.calcSystematicUncert( getBinEnergyMid(lowBin),getBinEnergyMid(highBin) ) )
        
        errSum = sumErrors(Errors)
        
       
        print("\nAsymm = %f\n"%uncert.A0)
        print("Weighted Energy Error: %f"%Errors[0])
        print("Total Statistical Error: %f"%Errors[1])
        print("Total Systematic Error: %f"%Errors[2])
        print("%i-%i keV total Error = %f\n"%(getBinEnergyLowEdge(lowBin),getBinEnergyUpperEdge(highBin),errSum))
        

