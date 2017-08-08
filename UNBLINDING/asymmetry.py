#!/usr/bin/python

import sys
import os
from math import *
sys.path.append("%s/systematics"%os.environ["ANALYSIS_CODE"])
from EnergyUncertainty.EnergyErrors import *
from CombineCorrelated import *

# I need to read in statistical uncertainties to store in the 

UNBLIND = False


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

        numer = sum([ ( (1./stats[i])**2 * data[i] if stats[i]!=0. else 0. )
                     for i in range(getEnergyBin(e0),getEnergyBin(e1)) if e0<=getBinEnergyMid(i)<=e1])
        denom = sum([ ( (1./stats[i])**2 if stats[i]!=0. else 0. )
                     for i in range(getEnergyBin(e0),getEnergyBin(e1)) if e0<=getBinEnergyMid(i)<=e1])
        
        return numer/denom if denom!=0. else 0.

def totalStatErr(data,dataErr,e0,e1):

    if len(data)!=len(dataErr):
        return 0

    numer = sum([ ((1./dataErr[i])**2 * data[i] if dataErr[i]!=0. else 0. )
                  for i in range(getEnergyBin(e0),getEnergyBin(e1)) if e0<=getBinEnergyMid(i)<=e1])
    weightsum = sum([ ( (1./dataErr[i])**2 if dataErr[i]!=0. else 0. )
                      for i in range(getEnergyBin(e0),getEnergyBin(e1)) if e0<=getBinEnergyMid(i)<=e1])
    A = numer/weightsum if weightsum>0. else 1.
    return fabs( ( 1./sqrt(weightsum) ) / A ) if weightsum>0. and fabs(A)>0. else 1.

def sumErrors(err):
    return sqrt(sum([er**2 for er in err]))

def getAsymmetry(anaCh,octLow,octHigh,elow,ehigh,corr="UnCorr",sim=False,withPOL=True):
    os.system("cd %s/Asymmetry; ./MBAnalyzer.exe %s %i %i %f %f %s %s 0 false %s %s > log.txt"%(os.environ["ANALYSIS_CODE"],anaCh,octLow, octHigh, elow, ehigh, corr, "true" if sim else "false", "true" if withPOL else "false",
                                                                                                "true" if UNBLIND else "false"))
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

def readOldAngleCorr(year=2011,anaCh="C",percErr=0.2):
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

def readOldBackscCorr(year=2011,anaCh="C",percErr=0.2,bsType="ALL"):
    A = [[],[]]
    yearString=None
    if year==2011:
        yearString = "2011-2012"
    else:
        yearString = "2012-2013"
        
    with open(os.environ["ANALYSIS_CODE"]+"/systematics/OldMCCorrection/"+
              "deltaBS%s_anaCh%s_%s.txt"%(bsType,anaCh,yearString),'rb') as tsvin:
        tsvin = csv.reader(tsvin, delimiter='\t')
        for row in tsvin:
            A[0].append(float(row[1]) if row[1]!="nan" and row[1]!="-nan" 
                        and row[1]!="inf" and row[1]!="-inf" else 0.)
            A[1].append(fabs(float(row[1])*percErr))
            #print("%s\t%s\t%s"%(row[0],row[1],row[2]))

    return A


def readAngleCorr(year=2011,anaCh="C",bsType="ALL"):
    A = [[],[]]
    yearString=None
    if year==2011:
        yearString = "2011-2012"
    else:
        yearString = "2012-2013"
    with open(os.environ["ANALYSIS_CODE"]+"/systematics/AngleCorrections/" +
              "%s_delta_3%s_anaCh%s.txt"
              %(yearString,(bsType if bsType!="ALL" else ""),anaCh),'rb') as tsvin:
        tsvin = csv.reader(tsvin, delimiter='\t')
        for row in tsvin:
            A[0].append(float(row[1]) if row[1]!="nan" and row[1]!="-nan" 
                        and row[1]!="inf" and row[1]!="-inf" else 0.)
            A[1].append(fabs(float(row[2])))
            #print("%s\t%s\t%s"%(row[0],row[1],row[2]))

    return A

def readBackscCorr(year=2011,anaCh="C",bsType="ALL"):
    A = [[],[]]
    yearString=None
    if year==2011:
        yearString = "2011-2012"
    else:
        yearString = "2012-2013"
    with open(os.environ["ANALYSIS_CODE"]+"/systematics/OldMCCorrection/" +
              "%s_delta_2%s_anaCh%s.txt"
              %(yearString,(bsType if bsType!="ALL" else ""),anaCh),'rb') as tsvin:
        tsvin = csv.reader(tsvin, delimiter='\t')
        for row in tsvin:
            A[0].append(float(row[1]) if row[1]!="nan" and row[1]!="-nan" 
                        and row[1]!="inf" and row[1]!="-inf" else 0.)
            A[1].append(fabs(float(row[2])))
            #print("%s\t%s\t%s"%(row[0],row[1],row[2]))

    return A



class uncertaintyHandler:
    
    def __init__(self,year,anaChoice="C",sim=False):
        self.sim = sim
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
        self.rad_err = []
        self.recoil_err = []
        self.theory_corr = []
        self.systA = []
        self.BS_corr = []
        self.Angle_corr = []

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
        from statistics alone""" #
        os.system("cd %s/Asymmetry; ./MBAnalyzer.exe %s %i %i %f %f %s %s 0 false %s %s 2&>1 log.txt"%(os.environ["ANALYSIS_CODE"],self.anaChoice,self.octLow,self.octHigh, 220, 670, 
                                                                                                    "UnCorr","true" if self.sim else "false", "false", "false"))
        infile = open( os.environ["%sANALYSIS_RESULTS"%("SIM_" if self.sim else "")]+"/Asymmetries/"+
                       "UnCorr_OctetAsymmetries_AnaCh%s_Octets_%i-%i_BinByBin.txt"
                       %(self.anaChoice,self.octLow,self.octHigh),'r' )

        print(os.environ["%sANALYSIS_RESULTS"%("SIM_" if self.sim else "")]+"/Asymmetries/"+
              "UnCorr_OctetAsymmetries_AnaCh%s_Octets_%i-%i_BinByBin.txt"
              %(self.anaChoice,self.octLow,self.octHigh))
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
                        percErr.append( 0. )

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

    def readTheoryUncertainties(self):
        """Read in the Energy uncertainties and store the percent 
        error from every bin."""
        
        infile = open( os.environ["ANALYSIS_RESULTS"]+"/Corrections/"+
                       "RadiativeUncertainty_%i.txt"%self.year)

        percErr = []

        if infile:
            for line in infile:
                if line[0]!="#":
                    l = line.split()
                    percErr.append( fabs(float(l[3])) )

        else: 
            print("Couldn't open file for Energy uncertainties")
            exit(0)
        

        self.rad_err = percErr
        infile.close()

        infile = open( os.environ["ANALYSIS_RESULTS"]+"/Corrections/"+
                       "RecoilOrderUncertainty_%i.txt"%self.year)

        percErr = []

        if infile:
            for line in infile:
                if line[0]!="#":
                    l = line.split()
                    percErr.append( fabs(float(l[3])) )

        else: 
            print("Couldn't open file for Energy uncertainties")
            exit(0)
        

        self.recoil_err = percErr

    def makeSystematicCorrections(self):
        """Read in the systematic correction bin-by-bin for final
        asymmetries, calculates delta_A/A, and then multiplies by the 
        error factor before filling the percent error for every bin."""

        self.BS_corr = readBackscCorr(self.year,self.anaChoice)
        self.Angle_corr = readAngleCorr(self.year,self.anaChoice)
        


    def makeTheoryUncertainties(self,lambda_err=0.):
        """Calculate the theory uncertainties given the error in lambda provided"""
        percErr = []
        
    def calcEnergyUncert(self,emin,emax):

        if len(self.energy_err)==0: 
            self.readEnergyUncertainties() 
            
        if len(self.stat_percent_err)==0: 
            self.statUncertainties() 
            
        return weightRealStats(self.energy_err,self.realAerr,emin,emax)

    def calcRadiativeUncert(self,emin,emax):

        if len(self.rad_err)==0: 
            self.readTheoryUncertainties() 
            
        if len(self.stat_percent_err)==0: 
            self.statUncertainties() 
            
        return weightRealStats(self.rad_err,self.realAerr,emin,emax)

    def calcRecoilUncert(self,emin,emax):

        if len(self.recoil_err)==0: 
            self.readTheoryUncertainties() 
            
        if len(self.stat_percent_err)==0: 
            self.statUncertainties() 
            
        return weightRealStats(self.recoil_err,self.realAerr,emin,emax)

    def calcPolarimetryCorr(self,emin,emax):
        

        UnCorr = getAsymmetry(self.anaChoice,self.octLow,self.octHigh,emin,
                              emax,"AllCorr",self.sim,False)
        PolCorr = getAsymmetry(self.anaChoice,self.octLow,self.octHigh,emin,
                               emax,"AllCorr",self.sim,True)
        
        deltaPol = (PolCorr[0]/UnCorr[0]-1.)
        #print("deltaPol: %f"%(deltaPol/(1.+deltaPol)))

        infile = open( os.environ["%sANALYSIS_RESULTS"%("SIM_" if self.sim else "")]+"/Corrections/"+
                       "PolUncertOctetAsymmetries_AnaCh%s_Octets_%i-%i.txt"%(self.anaChoice,self.octLow,self.octHigh))
        
        percErr = []
        
        if infile:
            for line in infile:
                if line[0]!="#":
                    l = line.split()
                    percErr.append( fabs(float(l[1])) )

        else: 
            print("Couldn't open file for depol uncertainties")
            exit(0)
        
        uncert = weightRealStats(percErr,self.realAerr,emin,emax)

        return [deltaPol/(1.+deltaPol),uncert]

                
    def calcStatUncert(self,emin,emax):

        #if len(self.stat_percent_err)==0: 
        #    self.statUncertainties() 
            
        #return totalStatErr(self.realA,self.realAerr,emin,emax) 
        asymm = getAsymmetry(self.anaChoice,self.octLow,self.octHigh,emin,
                             emax,"AllCorr",self.sim)
        self.A0 = asymm[0]
        asymm = getAsymmetry(self.anaChoice,self.octLow,self.octHigh,emin,
                             emax,"UnCorr",self.sim,False)
        return fabs( asymm[1]/asymm[0] )

    def calcMCCorr(self,emin,emax,percErrBS=0.25,percErrAngle=0.25):

        if len(self.BS_corr[0])==0: 
            self.makeSystematicCorrections() 
            
        if len(self.stat_percent_err)==0: 
            self.statUncertainties() 
            
        #return weightRealStats([self.syst_err[i]/fabs(self.syst_corr[i]+1) for i in range(len(self.syst_err))],self.realAerr,emin,emax)
        #return weightRealStats(self.syst_err,self.realAerr,emin,emax)
        
        UnCorr = getAsymmetry(self.anaChoice,self.octLow,self.octHigh,emin,
                              emax,"UnCorr",self.sim)
        AngleCorr = getAsymmetry(self.anaChoice,self.octLow,self.octHigh,emin,
                                 emax,"DeltaAngleOnly",self.sim)
        BSCorr = getAsymmetry(self.anaChoice,self.octLow,self.octHigh,emin,
                              emax,"DeltaBackscatterOnly",self.sim)

        deltaAngle = (AngleCorr[0]/UnCorr[0]-1.)
        deltaBS = (BSCorr[0]/UnCorr[0]-1.)
        deltaCorr = (1.+deltaAngle)*(1+deltaBS)-1.
        deltaCorrErr = sqrt( (percErrAngle*deltaAngle/(1.+deltaAngle))**2 + (percErrBS*deltaBS/(1.+deltaBS))**2)

        #print("deltaBS: %f +/- %f"%(deltaBS,percErrBS*deltaBS))
        #print("deltaAngle: %f +/- %f"%(deltaAngle,percErrAngle*deltaAngle))
        return [deltaCorr/(1+deltaCorr),deltaCorrErr/(1+deltaCorr)]

    def calcAngleCorr(self,emin,emax,BStype="ALL"):

        if len(self.BS_corr[0])==0: 
            self.makeSystematicCorrections() 
            
        if len(self.stat_percent_err)==0: 
            self.statUncertainties() 
            
        #return weightRealStats([self.syst_err[i]/fabs(self.syst_corr[i]+1) for i in range(len(self.syst_err))],self.realAerr,emin,emax)
        #return weightRealStats(self.syst_err,self.realAerr,emin,emax)
        
        UnCorr = getAsymmetry(self.anaChoice,self.octLow,self.octHigh,emin,
                              emax,"UnCorr",self.sim)
        AngleCorr = getAsymmetry(self.anaChoice,self.octLow,self.octHigh,emin,
                                 emax,"DeltaAngle%s"%BStype,self.sim)
        
        deltaAngle = (AngleCorr[0]/UnCorr[0]-1.)
        #print("deltaAngle: %f +/- %f"%(deltaAngle/(1.+deltaAngle),fabs(percErrAngle*deltaAngle/(1.+deltaAngle))))
        return [deltaAngle/(1.+deltaAngle),weightRealStats([self.Angle_corr[1][i]/(1.+self.Angle_corr[0][i]) for i in range(0,len(self.Angle_corr[0]))],self.realAerr,emin,emax)]#fabs(percErrAngle*deltaAngle/(1.+deltaAngle))]


    def calcBackscCorr(self,emin,emax, BStype="ALL"):

        if len(self.BS_corr[0])==0: 
            self.makeSystematicCorrections() 
            
        if len(self.stat_percent_err)==0: 
            self.statUncertainties() 
            
        #return weightRealStats([self.syst_err[i]/fabs(self.syst_corr[i]+1) for i in range(len(self.syst_err))],self.realAerr,emin,emax)
        #return weightRealStats(self.syst_err,self.realAerr,emin,emax)
        
        UnCorr = getAsymmetry(self.anaChoice,self.octLow,self.octHigh,emin,
                              emax,"UnCorr",self.sim)
        BackscCorr = getAsymmetry(self.anaChoice,self.octLow,self.octHigh,emin,
                                  emax,"DeltaBacksc%s"%BStype,self.sim)
        
        deltaBacksc = (BackscCorr[0]/UnCorr[0]-1.)
        #print("deltaBacksc: %f +/- %f"%(deltaBacksc/(1.+deltaBacksc),fabs(percErrBacksc*deltaBacksc/(1.+deltaBacksc))))
        return [deltaBacksc/(1.+deltaBacksc),weightRealStats([self.BS_corr[1][i]/(1.+self.BS_corr[0][i]) for i in range(0,len(self.BS_corr[0]))],self.realAerr,emin,emax)]#fabs(percErrBacksc*deltaBacksc/(1.+deltaBacksc))]


    def calcRecoilOrderCorr(self,emin,emax,percErr=0.02):
            
        if len(self.stat_percent_err)==0: 
            self.statUncertainties() 
            
        #return weightRealStats([self.syst_err[i]/fabs(self.syst_corr[i]+1) for i in range(len(self.syst_err))],self.realAerr,emin,emax)
        #return weightRealStats(self.syst_err,self.realAerr,emin,emax)
        
        UnCorr = getAsymmetry(self.anaChoice,self.octLow,self.octHigh,emin,
                              emax,"UnCorr",self.sim)
        Corr = getAsymmetry(self.anaChoice,self.octLow,self.octHigh,emin,
                             emax,"DeltaTheoryRecoil",self.sim)
        
        delta = (Corr[0]/UnCorr[0]-1.)
        #print("delta: %f +/- %f"%(delta/(1.+delta),fabs(percErr*delta/(1.+delta))))
        return [delta/(1.+delta),fabs(percErr*delta/(1.+delta))]

    
    def calcRadiativeCorr(self,emin,emax,percErr=0.02):

        if len(self.stat_percent_err)==0: 
            self.statUncertainties() 
            
        #return weightRealStats([self.syst_err[i]/fabs(self.syst_corr[i]+1) for i in range(len(self.syst_err))],self.realAerr,emin,emax)
        #return weightRealStats(self.syst_err,self.realAerr,emin,emax)
        
        UnCorr = getAsymmetry(self.anaChoice,self.octLow,self.octHigh,emin,
                              emax,"UnCorr",self.sim)
        Corr = getAsymmetry(self.anaChoice,self.octLow,self.octHigh,emin,
                             emax,"DeltaTheoryRadiative",self.sim)
        
        delta = (Corr[0]/UnCorr[0]-1.)
        #print("delta: %f +/- %f"%(delta/(1.+delta),fabs(percErr*delta/(1.+delta))))
        return [delta/(1.+delta),fabs(percErr*delta/(1.+delta))]


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
        
            
        return weightRealStats(gainErr,self.realAerr,emin,emax)



    def minimizer(self):
        """Minimizes the error and outputs all errors in order or largest to smallest"""

        self.statUncertainties()
        self.readEnergyUncertainties()
        self.makeSystematicCorrections()
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
                Errors.append( self.calcAngleCorr( getBinEnergyMid(lowBin),getBinEnergyMid(highBin) )[1] )
                Errors.append( self.calcBackscCorr( getBinEnergyMid(lowBin),getBinEnergyMid(highBin) )[1] )

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

       # minSystCorrTotal = weightRealStats(self.syst_corr,self.realAerr,getBinEnergyMid(enBinLow),getBinEnergyMid(enBinHigh))
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

    #2011: Energy Range = 220.-750. for bsErr=0.25 and angleErr=0.25 : 0.0064, but can go as low as 190-750
    #2012: Energy Range = 180.-740. for bsErr=0.25 and angleErr=0.25 : 0.00765. but same for 190-750

    #2011: Energy Range = 210.-750. for bsErr=0.2 and angleErr=0.25 : 0.006945
    #2012: Energy Range = 200.-740. for bsErr=0.2 and angleErr=0.25 : 0.00811
    # BUT!!! 2012 gives same uncertainty of 2011 energy range

    #2011: Energy Range = 230.-740. for bsErr=0.2 and angleErr=0.2
    #2012: Energy Range = 190.-760. for bsErr=0.2 and angleErr=0.2

    #also consider analysis choice "B" which seems to have the lowers Total uncert

    if 0:
    
        lowBin = 19
        highBin = 74
        
        uncert.statUncertainties()
        uncert.readEnergyUncertainties()
        uncert.readTheoryUncertainties()
        uncert.makeSystematicCorrections()
        
        Errors = []
        
        Errors.append( uncert.calcEnergyUncert(     getBinEnergyMid(lowBin),getBinEnergyMid(highBin) ) )
        Errors.append( uncert.calcStatUncert(       getBinEnergyLowEdge(lowBin),getBinEnergyUpperEdge(highBin) ) )
        Errors.append( uncert.calcAngleCorr( getBinEnergyMid(lowBin),getBinEnergyMid(highBin) )[1] )
        Errors.append( uncert.calcBackscCorr( getBinEnergyMid(lowBin),getBinEnergyMid(highBin) )[1] )
                
        errSum = sumErrors(Errors)
        
       
        print("\nAsymm = %f\n"%uncert.A0)
        print("Weighted Energy Error: %f"%Errors[0])
        print("Total Statistical Error: %f"%Errors[1])
        print("Total Angle Corr Error: %f"%Errors[2])
        print("Total backsc Corr Error: %f"%Errors[3])
        print("%i-%i keV total Error = %f\n"%(getBinEnergyLowEdge(lowBin),getBinEnergyUpperEdge(highBin),errSum))
        
        print("Weighted Radiative Corr Error: %f"%uncert.calcRadiativeUncert(getBinEnergyMid(lowBin),getBinEnergyMid(highBin)))
        print("Weighted Recoil Corr Error: %f"%uncert.calcRecoilUncert(getBinEnergyMid(lowBin),getBinEnergyMid(highBin)))
        print
        print("Total delta_20: %0.6f +/- %0.6f"%(uncert.calcBackscCorr( getBinEnergyMid(lowBin),getBinEnergyMid(highBin),"0" )[0],
                                                 uncert.calcBackscCorr( getBinEnergyMid(lowBin),getBinEnergyMid(highBin),"0" )[1]))
        print("Total delta_21: %0.6f +/- %0.6f"%(uncert.calcBackscCorr( getBinEnergyMid(lowBin),getBinEnergyMid(highBin),"1" )[0],
                                                 uncert.calcBackscCorr( getBinEnergyMid(lowBin),getBinEnergyMid(highBin),"1" )[1]))
        print("Total delta_22: %0.6f +/- %0.6f"%(uncert.calcBackscCorr( getBinEnergyMid(lowBin),getBinEnergyMid(highBin),"2" )[0],
                                                 uncert.calcBackscCorr( getBinEnergyMid(lowBin),getBinEnergyMid(highBin),"2" )[1]))
        print("Total delta_23: %0.6f +/- %0.6f"%(uncert.calcBackscCorr( getBinEnergyMid(lowBin),getBinEnergyMid(highBin),"3" )[0],
                                                 uncert.calcBackscCorr( getBinEnergyMid(lowBin),getBinEnergyMid(highBin),"3" )[1]))
        print("Total delta_2: %0.6f +/- %0.6f"%(uncert.calcBackscCorr( getBinEnergyMid(lowBin),getBinEnergyMid(highBin),"ALL" )[0],
                                                 uncert.calcBackscCorr( getBinEnergyMid(lowBin),getBinEnergyMid(highBin),"ALL" )[1]))

        print
        print("Total delta_30: %0.6f +/- %0.6f"%(uncert.calcAngleCorr( getBinEnergyMid(lowBin),getBinEnergyMid(highBin),"0" )[0],
                                                 uncert.calcAngleCorr( getBinEnergyMid(lowBin),getBinEnergyMid(highBin),"0" )[1]))
        print("Total delta_31: %0.6f +/- %0.6f"%(uncert.calcAngleCorr( getBinEnergyMid(lowBin),getBinEnergyMid(highBin),"1" )[0],
                                                 uncert.calcAngleCorr( getBinEnergyMid(lowBin),getBinEnergyMid(highBin),"1" )[1]))
        print("Total delta_32: %0.6f +/- %0.6f"%(uncert.calcAngleCorr( getBinEnergyMid(lowBin),getBinEnergyMid(highBin),"2" )[0],
                                                 uncert.calcAngleCorr( getBinEnergyMid(lowBin),getBinEnergyMid(highBin),"2" )[1]))
        print("Total delta_33: %0.6f +/- %0.6f"%(uncert.calcAngleCorr( getBinEnergyMid(lowBin),getBinEnergyMid(highBin),"3" )[0],
                                                 uncert.calcAngleCorr( getBinEnergyMid(lowBin),getBinEnergyMid(highBin),"3" )[1]))
        print("Total delta_3: %0.6f +/- %0.6f"%(uncert.calcAngleCorr( getBinEnergyMid(lowBin),getBinEnergyMid(highBin),"ALL" )[0],
                                                 uncert.calcAngleCorr( getBinEnergyMid(lowBin),getBinEnergyMid(highBin),"ALL" )[1]))
        print
        print("Total delta_MC: %0.6f +/- %0.6f"%(uncert.calcMCCorr( getBinEnergyMid(lowBin),getBinEnergyMid(highBin),0.25,0.25)[0],
                                                 uncert.calcMCCorr( getBinEnergyMid(lowBin),getBinEnergyMid(highBin),0.25,0.25)[1]))
             
        print
        delta20 = uncert.calcBackscCorr(getBinEnergyMid(lowBin),getBinEnergyMid(highBin),"0")
        delta21 = uncert.calcBackscCorr(getBinEnergyMid(lowBin),getBinEnergyMid(highBin),"1")
        delta22 = uncert.calcBackscCorr(getBinEnergyMid(lowBin),getBinEnergyMid(highBin),"2")
        delta23 = uncert.calcBackscCorr(getBinEnergyMid(lowBin),getBinEnergyMid(highBin),"3") 
        delta2 = uncert.calcBackscCorr(getBinEnergyMid(lowBin),getBinEnergyMid(highBin),"ALL") 
        delta30 = uncert.calcAngleCorr(getBinEnergyMid(lowBin),getBinEnergyMid(highBin),"0")
        delta31 = uncert.calcAngleCorr(getBinEnergyMid(lowBin),getBinEnergyMid(highBin),"1")
        delta32 = uncert.calcAngleCorr(getBinEnergyMid(lowBin),getBinEnergyMid(highBin),"2")
        delta33 = uncert.calcAngleCorr(getBinEnergyMid(lowBin),getBinEnergyMid(highBin),"3") 
        delta3 = uncert.calcAngleCorr(getBinEnergyMid(lowBin),getBinEnergyMid(highBin),"ALL") 


        indErrors = sumErrors([delta2[1],delta3[1]])
        print("Total Individual MC Uncert: %f +/- %f"%(uncert.calcMCCorr( getBinEnergyMid(lowBin),getBinEnergyMid(highBin))[0],
                                                       indErrors))
        
    if 1:
        anaChoice = "C"
        lowBin = 19
        highBin = 74
        
        errDeltaRecoil = 0.005 # These need to be fixed
        errDeltaRadiative = 0.005

        #### ************ make sure these are correct values
        deltaGainErr = 0.0018
        deltaMuonVetoErr = 0.0003
        deltaNeutronBG = 0.0001
        deltaNeutronBGErr = 0.0002
        
        ############### 2011-2012 #####################
        deltaEff2011 = 0.0013
        deltaEff2011Err = 0.0001

        deltaField2011 = 0.
        deltaField2011Err = 0.0012
        
        A2011 = uncertaintyHandler(2011,anaChoice)
        A2011.statUncertainties()
        A2011.readEnergyUncertainties()
        A2011.makeSystematicCorrections()

        statErr2011 = A2011.calcStatUncert(getBinEnergyMid(lowBin),getBinEnergyMid(highBin))
        energyErr2011 = A2011.calcEnergyUncert(getBinEnergyMid(lowBin),getBinEnergyMid(highBin))
        deltaBacksc2011 = A2011.calcBackscCorr(getBinEnergyLowEdge(lowBin),getBinEnergyUpperEdge(highBin))
        deltaAngle2011 = A2011.calcAngleCorr(getBinEnergyLowEdge(lowBin),getBinEnergyUpperEdge(highBin))
        deltaRecoil2011 = A2011.calcRecoilOrderCorr(getBinEnergyLowEdge(lowBin),getBinEnergyUpperEdge(highBin),errDeltaRecoil)
        deltaRadiative2011 = A2011.calcRadiativeCorr(getBinEnergyLowEdge(lowBin),getBinEnergyUpperEdge(highBin),errDeltaRadiative)
        
        A0_2011 = A2011.A0*(1.+deltaEff2011)*(1.+deltaField2011)*(1.+deltaNeutronBG/(1.-deltaNeutronBG))

        enUncert2011 = fabs(energyErr2011*A0_2011)
        statUncert2011 = fabs(statErr2011*A0_2011)
        angleUncert2011 = fabs(deltaAngle2011[1]*A0_2011)
        backscUncert2011 = fabs(deltaBacksc2011[1]*A0_2011)
        radiativeUncert2011 = fabs(deltaRadiative2011[1]*A0_2011)
        recoilUncert2011 = fabs(deltaRecoil2011[1]*A0_2011)
        fieldUncert2011 = fabs(deltaField2011Err/(1.+deltaField2011)*A0_2011)
        effUncert2011 = fabs(deltaEff2011Err/(1.+deltaEff2011)*A0_2011)
         
        depolCorr2011 = A2011.calcPolarimetryCorr(getBinEnergyMid(lowBin),getBinEnergyMid(highBin))


        #### 2011 uncertainty table
        print("\n\n****************** 2011-2012 ****************\n\n")
        print("A0 = %0.6f +/- %0.6f"%(A0_2011,sumErrors([enUncert2011,statUncert2011,angleUncert2011,backscUncert2011,radiativeUncert2011,
                                                         recoilUncert2011,fieldUncert2011,effUncert2011,depolCorr2011[1]*A0_2011,
                                                         deltaGainErr*A0_2011,deltaMuonVetoErr*A0_2011,deltaNeutronBGErr*A0_2011])))
        print("")

        print("\t\t% Corr\t\t% Unc.")
        print("depol\t\t%0.2f\t\t%0.2f"%(depolCorr2011[0]*100.,depolCorr2011[1]*100.))
        print("Energy\t\t\t\t%0.2f"%(energyErr2011*100.))
        print("BS\t\t%0.2f\t\t%0.2f"%(deltaBacksc2011[0]*100.,deltaBacksc2011[1]*100.))
        print("Angle\t\t%0.2f\t\t%0.2f"%(deltaAngle2011[0]*100.,deltaAngle2011[1]*100.))
        print("gain\t\t\t\t%0.2f"%(deltaGainErr*100.))
        print("MWPC eff\t%0.2f\t\t%0.2f"%(deltaEff2011/(1+deltaEff2011)*100.,deltaEff2011Err/(1.+deltaEff2011)*100.))
        print("field\t\t\t\t%0.2f"%(deltaField2011Err*100.))
        print("UCN Bg\t\t%0.2f\t\t%0.2f"%(deltaNeutronBG*100.,deltaNeutronBGErr*100.))
        print("muon\t\t\t\t%0.2f\n"%(deltaMuonVetoErr*100.))
        print("-----------------------------------------------")
        print("Stats\t\t\t\t%0.2f"%(statErr2011*100.))
        print("-----------------------------------------------")
        print("\tTheory Corrections")
        print("R.O.\t\t%0.2f\t\t%0.2f"%(deltaRecoil2011[0]*100.,0.03))#deltaRecoil[1]*100.))
        print("Rad\t\t%0.2f\t\t%0.2f"%(deltaRadiative2011[0]*100.,0.05) )#deltaRadiative[1]*100.))
        
        MC = MeasCombiner()
        
        mu0 = MC.new_meas_mu(A0_2011)
        stat0 = MC.new_meas_err(mu0,statUncert2011)	# stat
        err01 = MC.new_meas_err(mu0,enUncert2011) # syst - energy
        err02 = MC.new_meas_err(mu0,angleUncert2011) # syst - MC Angle
        err03 = MC.new_meas_err(mu0,backscUncert2011) # syst - MC Backsc
        #dst0 = MC.new_meas_err(mu0,0.00*A0_2011) # depol stat
        #err03 = MC.new_meas_err(mu0,0.000*A0_2011) # depol syst
        err04 = MC.new_meas_err(mu0,depolCorr2011[1]*A0_2011) # depol total
        err05 = MC.new_meas_err(mu0,effUncert2011) # MWPC eff
        err06 = MC.new_meas_err(mu0,fieldUncert2011) # Field Dip
        err07 = MC.new_meas_err(mu0,fabs(deltaGainErr*A0_2011)) # gain
        err08 = MC.new_meas_err(mu0,fabs(deltaMuonVetoErr*A0_2011)) # MuonVeto
        err09 = MC.new_meas_err(mu0,fabs(deltaNeutronBGErr*A0_2011)) # NeutronBG
        err010 = MC.new_meas_err(mu0,0.0003*A0_2011)#recoilUncert2011) # R.O. Corr
        err011 = MC.new_meas_err(mu0,0.0005*A0_2011)#radiativeUncert2011) # Radiative Corr
        
        

        ############### 2012-2013 #####################
        deltaEff2012 = 0.0012
        deltaEff2012Err = 0.0001

        deltaField2012 = 0.
        deltaField2012Err = 0.001
        
        A2012 = uncertaintyHandler(2012,anaChoice)
        A2012.statUncertainties()
        A2012.readEnergyUncertainties()
        A2012.makeSystematicCorrections()

        statErr2012 = A2012.calcStatUncert(getBinEnergyMid(lowBin),getBinEnergyMid(highBin))
        energyErr2012 = A2012.calcEnergyUncert(getBinEnergyMid(lowBin),getBinEnergyMid(highBin))
        deltaBacksc2012 = A2012.calcBackscCorr(getBinEnergyLowEdge(lowBin),getBinEnergyUpperEdge(highBin))
        deltaAngle2012 = A2012.calcAngleCorr(getBinEnergyLowEdge(lowBin),getBinEnergyUpperEdge(highBin))
        deltaRecoil2012 = A2012.calcRecoilOrderCorr(getBinEnergyLowEdge(lowBin),getBinEnergyUpperEdge(highBin),errDeltaRecoil)
        deltaRadiative2012 = A2012.calcRadiativeCorr(getBinEnergyLowEdge(lowBin),getBinEnergyUpperEdge(highBin),errDeltaRadiative)
        
        A0_2012 = A2012.A0*(1.+deltaEff2012)*(1.+deltaField2012)*(1.+deltaNeutronBG/(1.-deltaNeutronBG))

        enUncert2012 = fabs(energyErr2012*A0_2012)
        statUncert2012 = fabs(statErr2012*A0_2012)
        angleUncert2012 = fabs(deltaAngle2012[1]*A0_2012)
        backscUncert2012 = fabs(deltaBacksc2012[1]*A0_2012)
        radiativeUncert2012 = fabs(deltaRadiative2012[1]*A0_2012)
        recoilUncert2012 = fabs(deltaRecoil2012[1]*A0_2012)
        fieldUncert2012 = fabs(deltaField2012Err/(1.+deltaField2012)*A0_2012)
        effUncert2012 = fabs(deltaEff2012Err/(1.+deltaEff2012)*A0_2012)
         
        depolCorr2012 = A2012.calcPolarimetryCorr(getBinEnergyMid(lowBin),getBinEnergyMid(highBin))
                


        #### 2012 uncertainty table
        print("\n\n****************** 2012-2013 ****************\n\n")
        print("A0 = %0.6f +/- %0.6f"%(A0_2012,sumErrors([enUncert2012,statUncert2012,angleUncert2012,backscUncert2012,radiativeUncert2012,
                                                         recoilUncert2012,fieldUncert2012,effUncert2012,depolCorr2012[1]*A0_2012,
                                                         deltaGainErr*A0_2012,deltaMuonVetoErr*A0_2012,deltaNeutronBGErr*A0_2012])))
        print("")

        print("\t\t% Corr\t\t% Unc.")
        print("depol\t\t%0.2f\t\t%0.2f"%(depolCorr2012[0]*100.,depolCorr2012[1]*100.))
        print("Energy\t\t\t\t%0.2f"%(energyErr2012*100.))
        print("BS\t\t%0.2f\t\t%0.2f"%(deltaBacksc2012[0]*100.,deltaBacksc2012[1]*100.))
        print("Angle\t\t%0.2f\t\t%0.2f"%(deltaAngle2012[0]*100.,deltaAngle2012[1]*100.))
        print("gain\t\t\t\t%0.2f"%(deltaGainErr*100.))
        print("MWPC eff\t%0.2f\t\t%0.2f"%(deltaEff2012/(1+deltaEff2012)*100.,deltaEff2012Err/(1.+deltaEff2012)*100.))
        print("field\t\t\t\t%0.2f"%(deltaField2012Err*100.))
        print("UCN Bg\t\t%0.2f\t\t%0.2f"%(deltaNeutronBG*100.,deltaNeutronBGErr*100.))
        print("muon\t\t\t\t%0.2f\n"%(deltaMuonVetoErr*100.))
        print("-----------------------------------------------")
        print("Stats\t\t\t\t%0.2f"%(statErr2012*100.))
        print("-----------------------------------------------")
        print("\tTheory Corrections")
        print("R.O.\t\t%0.2f\t\t%0.2f"%(deltaRecoil2012[0]*100.,0.03))#deltaRecoil[1]*100.))
        print("Rad\t\t%0.2f\t\t%0.2f"%(deltaRadiative2012[0]*100.,0.05))#deltaRadiative[1]*100.))
        

        mu1 = MC.new_meas_mu(A0_2012)
        stat1 = MC.new_meas_err(mu1,statUncert2012)	# stat
        err11 = MC.new_meas_err(mu1,enUncert2012) # syst - energy
        err12 = MC.new_meas_err(mu1,angleUncert2012) # syst - MC Angle
        err13 = MC.new_meas_err(mu1,backscUncert2012) # syst - MC Backsc
        #dst1 = MC.new_meas_err(mu1,0.00*A0_2012) # depol stat
        #err13 = MC.new_meas_err(mu1,0.000*A0_2012) # depol syst
        err14 = MC.new_meas_err(mu1,depolCorr2012[1]*A0_2012) # depol total
        err15 = MC.new_meas_err(mu1,effUncert2012) # MWPC eff
        err16 = MC.new_meas_err(mu1,fieldUncert2012) # Field Dip
        err17 = MC.new_meas_err(mu1,fabs(deltaGainErr*A0_2012)) # gain
        err18 = MC.new_meas_err(mu1,fabs(deltaMuonVetoErr*A0_2012)) # MuonVeto
        err19 = MC.new_meas_err(mu1,fabs(deltaNeutronBGErr*A0_2012)) # NeutronBG
        err110 = MC.new_meas_err(mu1,0.0003*A0_2012)#recoilUncert2012) # R.O. Corr
        err111 = MC.new_meas_err(mu1,0.0005*A0_2012)#radiativeUncert2012) # Radiative Corr
        
    
	
        MC.add_correlation(err01,err11,1.0)
        MC.add_correlation(err02,err12,1.)
        MC.add_correlation(err03,err13,1.)
        MC.add_correlation(err04,err14,1.0)
        MC.add_correlation(err05,err15,1.0)
        MC.add_correlation(err06,err16,1.0)
        MC.add_correlation(err07,err17,1.0)
        MC.add_correlation(err08,err18,1.0)
        MC.add_correlation(err09,err19,1.0)
        MC.add_correlation(err010,err110,1.0)
        MC.add_correlation(err011,err111,1.0)

        result = MC.calc_combo()
        

        #### 2012 uncertainty table
        print("\n\n****************** FINAL RESULT ****************\n\n")
        print("A0 = %0.6f +/- %0.6f"%(result[0],result[1]))
        print "stat =",MC.errcombo([stat0,stat1])
        print "syst =",MC.errcombo([err01,err11,err02,err12,err03,err13,err04,err14,err05,err15,err06,err16,err07,err17,err08,err18,err09,err19,err010,err110,err011,err111])
        print("")

        print("\t\t% Corr2011\t% Corr2012\t% Unc.")
        print("depol\t\t%0.2f\t\t%0.2f\t\t%0.2f"%(depolCorr2011[0]*100.,depolCorr2012[0]*100.,fabs(MC.errcombo([err04,err14])/result[0])*100.))
        print("Energy\t\t\t\t\t\t%0.2f"%(fabs(MC.errcombo([err01,err11])/result[0])*100.))
        print("BS\t\t%0.2f\t\t%0.2f\t\t%0.2f"%(deltaBacksc2011[0]*100.,deltaBacksc2012[0]*100.,fabs(MC.errcombo([err03,err13])/result[0])*100.))
        print("Angle\t\t%0.2f\t\t%0.2f\t\t%0.2f"%(deltaAngle2011[0]*100.,deltaAngle2012[0]*100.,fabs(MC.errcombo([err02,err12])/result[0])*100.))
        print("gain\t\t\t\t\t\t%0.2f"%(fabs(MC.errcombo([err07,err17])/result[0])*100.))
        print("MWPC eff\t%0.2f\t\t%0.2f\t\t%0.2f"%(deltaEff2011/(1+deltaEff2011)*100.,deltaEff2012/(1+deltaEff2012)*100.,fabs(MC.errcombo([err05,err15])/result[0])*100.))
        print("field\t\t\t\t\t\t%0.2f"%(fabs(MC.errcombo([err06,err16])/result[0])*100.))
        print("UCN Bg\t\t%0.2f\t\t%0.2f\t\t%0.2f"%(deltaNeutronBG*100.,deltaNeutronBG*100.,fabs(MC.errcombo([err09,err19])/result[0])*100.))
        print("muon\t\t\t\t\t\t%0.2f\n"%(fabs(MC.errcombo([err08,err18])/result[0])*100.))
        print("-----------------------------------------------")
        print("Stats\t\t\t\t\t\t%0.2f"%(fabs(MC.errcombo([stat0,stat1])/result[0])*100.))
        print("-----------------------------------------------")
        print("\tTheory Corrections")
        print("R.O.\t\t%0.2f\t\t%0.2f\t\t%0.2f"%(deltaRecoil2011[0]*100.,deltaRecoil2012[0]*100.,fabs(MC.errcombo([err010,err110])/result[0])*100.))#deltaRecoil[1]*100.))
        print("Rad\t\t%0.2f\t\t%0.2f\t\t%0.2f"%(deltaRadiative2011[0]*100.,deltaRadiative2012[0]*100.,fabs(MC.errcombo([err011,err111])/result[0])*100.))#deltaRadiative[1]*100.))
