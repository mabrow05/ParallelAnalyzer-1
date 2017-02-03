#!/usr/bin/python

import sys
import os
sys.path.append("..")
from math import *
from ucnacore.PyxUtils import *
from ucnacore.DecayPhysics import *
#from review.ErrorEnvelope import *
from CorrectionsPlotter import *

# These are the error envelopes
limdat = {2008:[(0,5.0),(250,5.0),(500,500*0.013),(900,900*0.025),(1000,1000*0.025),(1200,1200*0.025)],
			2010:[(0,2.5),(200,200*0.0125),(500,500*0.0125),(1000,500*0.0125)],
			2011:[(0,0.017*130.3),(130.3,130.3*0.017),(368.49,0.014*368.49),(993.789,993.789*0.007),(1200,993.789*0.007)], #Parallel Analyzer
	                20112:[(0,6.),(200,6.),(500,10.),(1000,10.)], #Conservative envelope
	                #2011:[(0,200.*0.02),(200,200*0.02),(500,500*0.0150),(1000,500*0.0150)],#liberal envelope
			2012:[(0,0.02*180.3),(180.3,180.3*0.02),(418.4938,0.018*418.4938),(1200,1200.*0.011)]} #Parallel Analyzer

def bin_edges(w=10,n=100):
	return [w*i for i in range(n+1)]

#baseOutPath = os.environ["UCNA_AUX"]+"/Corrections/"
baseOutPath = os.environ["UCNA_ANALYSIS_OUTPUT_DIR"]
os.system("mkdir %s/Corrections"%baseOutPath)

def writeUncertaintyTable(fout,dat):
	"""Write uncertainties table to file in correct format"""
	fout.write("#E_lo\tE_hi\tcorrection\tuncertainty\n")
	for d in dat:
		fout.write("%i\t%i\t%f\t%f\n"%tuple(d))

def ObsAsymApprox(KE,year):
	"""Phenomenological fit to MC observed asymmetry"""
	A0 = .1172
	p2 = .966648
	p3 = 0.0001174
	p4 = 56.5529
	p5 = 12.8861
	p6 = 2.18673
	return A0*beta(KE)*0.5*p2*(1+p3*KE)*(1+p6/(1+exp((KE-p4)/p5)))

def simpleAsym(KE,year):
	return A0_PDG*beta(KE)*0.5

def energyErrorA(E,year):
	Eprim = E+calEnvelope(E,year)
	return ObsAsymApprox(Eprim,year)/ObsAsymApprox(E,year)-1.

def energyErrorSimple(E,year):
	Eprim = E+calEnvelope(E,year)
	return simpleAsym(Eprim,year)/simpleAsym(E,year)-1.

def energyErrorRC(E,year):
	Eprim = E+calEnvelope(E,year)
	return simpleAsym(Eprim,year)*(1+WilkinsonRWM(Eprim))/(simpleAsym(E,year)*(1+WilkinsonRWM(E)))-1.


def linearityUncertaintyTable(year=2011):
	"""Uncertainty due to energy calibration errors; see EnergyErrorsRevis.pdf"""
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
		print c,Eprim,ObsAsymApprox(c,year),err,errmax
	dat = dat[::-1]
	fout = open(baseOutPath+"/Corrections/EnergyLinearityUncertainty_%i.txt"%year,"w")
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
			 
	gdat = [ [x,100*energyErrorA(x,year),100*energyErrorSimple(x,year),100*energyErrorRC(x,year)] for x in unifrange(50,850.,800) ]
	gdat2010 = [ [x,100*energyErrorA(x,2010),100*energyErrorSimple(x,2010),100*energyErrorRC(x,2010)] for x in unifrange(50,850.,800) ]
	gdatComp = [ [x,100*energyErrorA(x,2011),100*energyErrorSimple(x,2011),100*energyErrorRC(x,2011)] for x in unifrange(50,850.,800) ]
	#gCx.plot(graph.data.points(gdat[::8],x=1,y=3,title="$A={\\beta \\over 2}A_0$"),
	#			 [ graph.style.line([style.linewidth.THick,style.linestyle.dotted]),])
	#gCx.plot(graph.data.points(gdat,x=1,y=4,title="$A={\\beta \\over 2}(1+$R.C.$)A_0$"),
	#			 [ graph.style.line([style.linewidth.THick,style.linestyle.dashed]),])
	#gCx.plot(graph.data.points(gdat,x=1,y=2,title="MC asymmetry"),
	#		 [ graph.style.line([style.linewidth.THick]),])
	gCx.plot(graph.data.points(gdat2010,x=1,y=2,title="2010"),
			 [ graph.style.line([style.linewidth.THick]),])
	gCx.plot(graph.data.points(gdat,x=1,y=2,title="2011-2012 Parallel Analyzer"),
			 [ graph.style.line([style.linewidth.THick,color.rgb.red]),])
	#gCx.plot(graph.data.points(gdat,x=1,y=2,title="2012-2013 Parallel Analyzer"),
	#		 [ graph.style.line([style.linewidth.THick,color.rgb.blue]),])
			 
	print "Eavg MC MPM 2010=",weightStats(gdat2010,220,670)
        print "Eavg MC MB %i="%year,weightStats(gdat,220,670)
	
	#print "Eavg plain = ",weightStats([(x[0],x[2]) for x in gdat],220,670)
				 			 
	gCx.writetofile("%s/Corrections/EnergyUncert%i_new.pdf"%(baseOutPath,year))


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

 	print "Eavg MC =",weightStats(gdat,220,670)
	
	gCx.writetofile("/Users/michael/Desktop/GainfluctsUncert.pdf")



if __name__=="__main__":
	year = 2011
	linearityUncertaintyTable(year)
	#gainFluctsUncertaintyTable()
	plotEnergyErrors(year)
	#plotGainfluctErrors()
