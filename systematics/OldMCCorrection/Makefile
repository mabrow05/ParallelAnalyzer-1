ObjSuf		= o
SrcSuf		= C
DllSuf		= so
OutPutOpt	= -o
ROOTLIBS	= $(shell root-config --libs)
ROOTGLIBS	= $(shell root-config --glibs)
ROOTCFLAGS	= $(shell root-config --cflags)

VPATH = ./:../../include

CXX		= g++
CXXFLAGS	= -O2 -Wall -fPIC -std=c++11
LD		= g++
SOFLAGS		= -shared
LIBS		= $(ROOTLIBS) $(ROOTGLIBS)
CXXFLAGS       += $(ROOTCFLAGS) -I../../include
LIBS	       += -lSpectrum -lMinuit

objects = Asymmetries.o MBUtils.o EvtRateHandler.o calibrationTools.o


.PHONY: all
all: MBAnalyzer.exe UK2SIMcomparison.exe

MBAnalyzer.exe: MBAnalyzer.o $(objects) BetaSpectrum.o
	$(CXX) $(CXXFLAGS) -o MBAnalyzer.exe MBAnalyzer.o $(objects) BetaSpectrum.o $(LIBS) 
	@echo "Done"


calibrationTools.o: calibrationTools.hh

Asymmetries.o: Asymmetries.hh MBUtils.o

BetaSpectrum.o: BetaSpectrum.hh

EvtRateHandler.o: EvtRateHandler.hh  MBUtils.o

MBUtils.o: MBUtils.hh

clean:
		@rm -f *.o *~  core
