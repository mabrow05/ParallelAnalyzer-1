#!/usr/bin/python

import sys
import csv
import os
from math import *
from PyxUtils import *
from DecayPhysics import *
#from review.ErrorEnvelope import *
from CorrectionsPlotter import *

# These are the error envelopes
limdat = {2008:[(0,5.0),(250,5.0),(500,500*0.013),(900,900*0.025),(1000,1000*0.025),(1200,1200*0.025)],
          2010:[(0,2.5),(200,200*0.0125),(500,500*0.0125),(1000,500*0.0125)],
          2011:[(0,0.017*130.3),(130.3,130.3*0.017),(368.49,0.010*368.49),(993.789,993.789*0.0065),(1200,993.789*0.0065)], 
          2012:[(0,0.017*130.3),(130.3,130.3*0.017),(368.49,0.010*368.49),(993.789,993.789*0.007),(1200,993.789*0.007)],
          #2012:[(0,0.017*130.3+10.),(130.3,130.3*0.017+10.),(368.49,0.010*368.49+10.),(993.789,993.789*0.0065+10.),(1200,993.789*0.0065+10.)]
}


def readCalEnvelope(year=2011,upper=True):
        envLower = []
        envUpper = []
        env = []
        with open('envolopeValues-%i-deg2-cal4-curves1000.tsv'%year,'rb') as tsvin:
                tsvin = csv.reader(tsvin, delimiter='\t')
                for row in tsvin:
                        envLower.append( ( float(row[0]),fabs(float(row[1])) ) )
                        envUpper.append( ( float(row[0]),float(row[2]) ) )
                        env.append( ( float(row[0]),max([fabs(float(row[1])),fabs(float(row[2]))]) ) )
                        print("%s\t%s\t%s"%(row[0],row[1],row[2]))
        limdat[year*10+1] = envLower
        limdat[year*10+2] = envUpper
        limdat[year] = env
        #if upper:
        #        limdat[year] = envUpper
        #else:
        #        limdat[year] = envLower
                        
def calEnvelope(E,year=2011):	
	i = 0
	while E > limdat[year][i+1][0]:
		i+=1
	l = (E-limdat[year][i][0])/(limdat[year][i+1][0]-limdat[year][i][0])
	return (1-l)*limdat[year][i][1]+l*limdat[year][i+1][1]

# uniformly spaced points
def unifrange(xmin,xmax,npts):
	return [ xmin + float(i)/float(npts-1)*(xmax-xmin) for i in range(npts)]


def bin_edges(w=10,n=100):
	return [w*i for i in range(n+1)]

baseOutPath = os.environ["ANALYSIS_RESULTS"]
os.system("mkdir %s/Corrections"%baseOutPath)

def writeUncertaintyTable(fout,dat):
	"""Write uncertainties table to file in correct format"""
	fout.write("#E_lo\tE_hi\tcorrection\tuncertainty\n")
	for d in dat:
		fout.write("%i\t%i\t%f\t%f\n"%tuple(d))

def ObsAsymApprox(KE,year):
	"""Phenomenological fit to MC observed asymmetry"""
        A0=0.
        p2=0.
        p3=0.
        p4=0.
        p5=0.
        p6=0.
        
        #MPM Values
        if year==2010 or year/10==2010:
                A0 = .1172
                p2 = .966648
                p3 = 0.0001174
                p4 = 56.5529
                p5 = 12.8861
                p6 = 2.18673

        #MB Values
        
        elif year==2011 or year/10==2011:
                A0 = .1184
                p2 = 0.967004   #[0]
                p3 = 0.000125724 #[1]
                p4 = 69.4554    #[3]
                p5 = 40.01147   #[4]
                p6 = 0.158159  #[2]

        elif year==2012 or year/10==2012:
                A0 = .1184
                p2 = 0.966143   #[0]
                p3 = 0.000130679 #[1]
                p4 = -230.862    #[3]
                p5 = 47.5465   #[4]
                p6 = 65.2350  #[2]

        # Functional Form is A0*beta(KE)*0.5*p2*(1+p3*KE)*(1+p6/(1+exp((KE-p4)/p5)))

        try: 
                exp_check = exp((KE-p4)/p5) 
        except OverflowError:
                exp_check = 10000000.

	return A0*beta(KE)*0.5*p2*(1+p3*KE)*(1+p6/(1+exp_check))

def simpleAsym(KE,year):
	return A0_PDG*beta(KE)*0.5

def energyErrorA(E,year):
	Eprim = E+calEnvelope(E,year) # the minus side gives larger uncert
	return ObsAsymApprox(Eprim,year)/ObsAsymApprox(E,year)-1.#simpleAsym(Eprim,year)/simpleAsym(E,year)-1.


def energyErrorSimple(E,year):
	Eprim = E+calEnvelope(E,year)
	return simpleAsym(Eprim,year)/simpleAsym(E,year)-1.

def energyErrorRC(E,year):
	Eprim = E+calEnvelope(E,year)
	return (WilkinsonRWM(E)-WilkinsonRWM(Eprim))/(1+WilkinsonRWM(E))

def energyErrorRad(E,year):
        Eprim = E+calEnvelope(E,year)
        return (shann_h_minus_g_a2pi(E)-shann_h_minus_g_a2pi(Eprim))/(1+shann_h_minus_g_a2pi(E)) 
	#return simpleAsym(Eprim,year)*(1+shann_h_minus_g_a2pi(Eprim))/(simpleAsym(E,year)*(1+shann_h_minus_g_a2pi(E)))-1.

def theoryUncertaintyTable(year=2011):
	"""Uncertainty due to energy calibration errors; see EnergyErrorsRevis.pdf"""
        
        year2 = year*10+1

	edges = bin_edges()
	dat = []
	errmax = 0
	for i in range(len(edges)-1)[-1::-1]:
		c = 0.5*(edges[i]+edges[i+1])
                Eprim = c+calEnvelope(c,year2)
		err = energyErrorRC(c,year2)
		if err > errmax:
			errmax = err
		#dat.append((edges[i],edges[i+1],0.0,errmax))
		dat.append((edges[i],edges[i+1],0.0,err))
		print c,Eprim,ObsAsymApprox(c,year2),err,errmax
                #print c,Eprim,simpleAsym(c,year),err,errmax
	dat = dat[::-1]
	fout = open(baseOutPath+"/Corrections/RecoilOrderUncertainty_%i.txt"%year,"w")
	fout.write("# Uncertainty from R.O. for %i data\n"%year)
	writeUncertaintyTable(fout,dat)

        dat = []
	errmax = 0
        for i in range(len(edges)-1)[-1::-1]:
		c = 0.5*(edges[i]+edges[i+1])
                Eprim = c+calEnvelope(c,year2)
		err = energyErrorRad(c,year2)
		if err > errmax:
			errmax = err
		#dat.append((edges[i],edges[i+1],0.0,errmax))
		dat.append((edges[i],edges[i+1],0.0,err))
		print c,Eprim,ObsAsymApprox(c,year2),err,errmax
                #print c,Eprim,simpleAsym(c,year),err,errmax
	dat = dat[::-1]
	fout = open(baseOutPath+"/Corrections/RadiativeUncertainty_%i.txt"%year,"w")
	fout.write("# Uncertainty from Radiative Corr for %i data\n"%year)
	writeUncertaintyTable(fout,dat)

def plotTheoryErrors(year=2011):

	gCx=graph.graphxy(width=15,height=8,
					  x=graph.axis.lin(title="Energy [keV]",min=0,max=780),
					  y=graph.axis.lin(title="uncertainty $\\Delta A/A$ [\\%]",min=0,max=3.),
					  key = graph.key.key(pos="tl"))
	setTexrunner(gCx)
			 
        gdat = [ [x,100*WilkinsonRWM(x)/(1+WilkinsonRWM(x)),100*shann_h_minus_g_a2pi(x)/(1+shann_h_minus_g_a2pi(x)) ] for x in unifrange(5.,795.,79)  ]


        gCx.plot(graph.data.points(gdat,x=1,y=2,title="Recoil Order"),
                 [ graph.style.line([style.linewidth.THick,color.rgb.red]),])
        gCx.plot(graph.data.points(gdat,x=1,y=3,title="Radiative"),
                 [ graph.style.line([style.linewidth.THick,style.linestyle.dotted,color.rgb.blue]),])
       
	
	print "%s Eavg Radiative = "%year,weightStats(gdat,220,670)
        print "%s Eavg Recoil = "%year,weightStats( [ [ x[0],x[2] ] for x in gdat ] ,220,670)

	gCx.writetofile("%s/Corrections/TheoryUncert%i.pdf"%(baseOutPath,year))

def linearityUncertaintyTable(year=2011):
	"""Uncertainty due to energy calibration errors; see EnergyErrorsRevis.pdf"""
        
        year2 = year*10+1

	edges = bin_edges()
	dat = []
	errmax = 0
	for i in range(len(edges)-1)[-1::-1]:
		c = 0.5*(edges[i]+edges[i+1])
		Eprim = c+calEnvelope(c,year)
		err = ObsAsymApprox(Eprim,year)/ObsAsymApprox(c,year)-1
		if err > errmax:
			errmax = err
		#dat.append((edges[i],edges[i+1],0.0,errmax))
		dat.append((edges[i],edges[i+1],0.0,err))
		print c,Eprim,ObsAsymApprox(c,year2),err,errmax
                #print c,Eprim,simpleAsym(c,year),err,errmax
	dat = dat[::-1]
	fout = open(baseOutPath+"/Corrections/EnergyLinearityUncertainty_%i.txt"%year,"w")
	fout.write("# Uncertainty from energy calibration linearity envelope for %i data\n"%year)
	writeUncertaintyTable(fout,dat)

def gainUncertaintyTable(year=2011,gainErr=0.005):
	"""Uncertainty due to gain calibration errors; see EnergyErrorsRevis.pdf"""
	edges = bin_edges()
	dat = []
	errmax = 0
	for i in range(len(edges)-1)[-1::-1]:
		c = 0.5*(edges[i]+edges[i+1])
		Eprim = c+gainErr*c
		err = ObsAsymApprox(Eprim,year)/ObsAsymApprox(c,year)-1
		if err > errmax:
			errmax = err
		#dat.append((edges[i],edges[i+1],0.0,errmax))
		dat.append((edges[i],edges[i+1],0.0,err))
		print c,Eprim,ObsAsymApprox(c,year),err,errmax
                #print c,Eprim,simpleAsym(c,year),err,errmax
	dat = dat[::-1]
	fout = open(baseOutPath+"/Corrections/gainUncertainty_%i.txt"%year,"w")
	fout.write("# Uncertainty from energy calibration linearity envelope for %i data\n"%year)
	writeUncertaintyTable(fout,dat)


def weightStats(xydat,e0,e1):
	s0sum = sum([beta(x[0])**2 * S0(x[0]) for x in xydat if e0<=x[0]<=e1])
	xsum = sum([x[1] * beta(x[0])**2 * S0(x[0]) for x in xydat if e0<=x[0]<=e1])
	return xsum/s0sum
	
def plotEnergyErrors(year=2011):

	gCx=graph.graphxy(width=15,height=8,
					  x=graph.axis.lin(title="Energy [keV]",min=0,max=800),
					  y=graph.axis.lin(title="uncertainty $\\Delta A/A$ [\\%]",min=0,max=1.),
					  key = graph.key.key(pos="tr"))
	setTexrunner(gCx)
			 
        #gdat = [ [x,100*fabs(energyErrorA(x,year)),100*energyErrorSimple(x,year),100*energyErrorRC(x,year)] for x in unifrange(50,850.,800) ]

	#gdatLower = [ [x,100*fabs(energyErrorA(x,year*10+1)),100*energyErrorSimple(x,year*10+1),100*energyErrorRC(x,year*10+1)] for x in unifrange(50,850.,800) ]
	#gdatUpper = [ [x,100*fabs(energyErrorA(x,year*10+2)),100*energyErrorSimple(x,year*10+2),100*energyErrorRC(x,year*10+2)] for x in unifrange(50,850.,800) ]
	gdat2010 = [ [x,100*energyErrorA(x,2010),100*energyErrorSimple(x,2010),100*energyErrorRC(x,2010)] for x in unifrange(50,850.,800) ]
	gdat2011 = [ [x,100*energyErrorA(x,2011),100*energyErrorSimple(x,2011),100*energyErrorRC(x,2011)] for x in unifrange(50,850.,800) ]
	gdat2012 = [ [x,100*energyErrorA(x,2012),100*energyErrorSimple(x,2012),100*energyErrorRC(x,2012)] for x in unifrange(50,850.,800) ]

	#gCx.plot(graph.data.points(gdat[::8],x=1,y=3,title="$A={\\beta \\over 2}A_0$"),
        #         [ graph.style.line([style.linewidth.THick,style.linestyle.dotted]),])
	#gCx.plot(graph.data.points(gdat,x=1,y=4,title="$A={\\beta \\over 2}(1+$R.C.$)A_0$"),
        #         [ graph.style.line([style.linewidth.THick,style.linestyle.dashed]),] )
        #gCx.plot(graph.data.points(gdat,x=1,y=2,title="%i Old Method"%year),
        #         [ graph.style.line([style.linewidth.THick]),])
	#gCx.plot(graph.data.points(gdatUpper,x=1,y=2,title="%i Upper Edge"%year),
        #         [ graph.style.line([style.linewidth.THick,color.rgb.red]),])
        #gCx.plot(graph.data.points(gdatLower,x=1,y=2,title="%i New Method"%year),
        #         [ graph.style.line([style.linewidth.THick,color.rgb.blue]),])
        #gCx.plot(graph.data.points(gdat2010,x=1,y=2,title="2010"),
        #         [ graph.style.line([style.linewidth.THick,color.rgb.red]),])
	gCx.plot(graph.data.points(gdat2011,x=1,y=2,title="2011-2012"),
                 [ graph.style.line([style.linewidth.THick,color.rgb.blue]),])
	gCx.plot(graph.data.points(gdat2012,x=1,y=2,title="2012-2013"),
                 [ graph.style.line([style.linewidth.THick,color.rgb.green]),])
			 
	#print "Eavg MC MPM 2010=",weightStats(gdat2010,220,670)
        print "Eavg MC MB 2011=",weightStats(gdat2011,190,740)
        print "Eavg MC MB 2012=",weightStats(gdat2012,190,740)
	
	#print "Eavg MC %i = "%year,weightStats(gdat,220,670)
	#print "Eavg MC Lower %i = "%year,weightStats(gdatLower,220,670)
        #print "Eavg MC Upper %i = "%year,weightStats(gdatUpper,220,670)
	gCx.writetofile("%s/Corrections/EnergyUncert%s.pdf"%(baseOutPath,"ALL" ))


def plotGainfluctErrors():

	gCx=graph.graphxy(width=15,height=8,
					  x=graph.axis.lin(title="Energy [keV]",min=0,max=800),
					  y=graph.axis.lin(title="uncertainty $\\Delta A/A$ [\\%]",min=-1,max=1),
					  key = graph.key.key(pos="tr"))
	setTexrunner(gCx)
	
	
	gfR = (lambda KE,delta:  SPol(KE,1)*SPol(KE*(1+delta),1)/SPol(KE,-1)/SPol(KE*(1-delta),-1))
	Asr = (lambda R: (1-sqrt(R))/(1+sqrt(R)))
	gfErr = (lambda KE,delta: Asr(gfR(KE,delta))/Asr(gfR(KE,0))-1. )
	
	gdat = [ [x,100*gfErr(x,-0.000125)] for x in unifrange(1,781.,100) ]
	gCx.plot(graph.data.points(gdat,x=1,y=2,title="spectra from theory"),
				 [ graph.style.line([style.linewidth.THick,style.linestyle.dotted]),])
	
	gfl = CorrFile(baseCorrPath+"GainFlucts.txt")
	gdat = [ [0.5*(d[0]+d[1]),100*d[3]] for d in gfl.dat if 20<d[0]<700]
	gCx.plot(graph.data.points(gdat,x=1,y=2,title="spectra from data"),
			 [ graph.style.line([style.linewidth.THick]),])
			 
	gCx.plot(graph.data.function("y(x)=0",title=None), [graph.style.line(lineattrs=[style.linestyle.dashed,]),])

 	print "Eavg MC =",weightStats(gdat,230,750)
	
	gCx.writetofile("/Users/michael/Desktop/GainfluctsUncert.pdf")



if __name__=="__main__":
	year = 2011
        #readCalEnvelope(year)
        #theoryUncertaintyTable(year)
        plotTheoryErrors(year)
        #gainUncertaintyTable(year,0.0064)
	#linearityUncertaintyTable(year)
	#gainFluctsUncertaintyTable()
        #readCalEnvelope(2011)
        #readCalEnvelope(2012)
	#plotEnergyErrors(20112)
	#plotGainfluctErrors()
