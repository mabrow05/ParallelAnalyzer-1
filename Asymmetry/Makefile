ObjSuf		= o
SrcSuf		= C
DllSuf		= so
OutPutOpt	= -o
ROOTLIBS	= $(shell root-config --libs)
ROOTGLIBS	= $(shell root-config --glibs)
ROOTCFLAGS	= $(shell root-config --cflags)


CXX		= g++
CXXFLAGS	= -O2 -Wall -fPIC -std=c++11 -g
LD		= g++
SOFLAGS		= -shared
LIBS		= $(ROOTLIBS) $(ROOTGLIBS)
CXXFLAGS       += $(ROOTCFLAGS) -Iinclude
LIBS	       += -lSpectrum -lMinuit

objects = MBUtils.o Asymmetries.o SQLinterface.o SystematicCorrections.o EvtRateHandler.o 

.PHONY: all
all: MBAnalyzer.exe UK2SIMcomparison.exe

MBAnalyzer.exe: MBAnalyzer.o $(objects) 
	$(CXX) $(CXXFLAGS) -o MBAnalyzer.exe MBAnalyzer.o $(objects) $(LIBS) 
	@echo "Done"

UK2MPMcomparison.exe: UK2MPMcomparison.o $(objects) 
	$(CXX) $(CXXFLAGS) -o UK2MPMcomparison.exe UK2MPMcomparison.o $(objects) $(LIBS) 
	@echo "Done"

UK2SIMcomparison.exe: UK2SIMcomparison.o $(objects) 
	$(CXX) $(CXXFLAGS) -o UK2SIMcomparison.exe UK2SIMcomparison.o $(objects) $(LIBS) 
	@echo "Done"

Asymmetries.o: Asymmetries.hh

EvtRateHandler.o: EvtRateHandler.hh SQLinterface.hh

SQLinterface.o: SQLinterface.hh

SystematicCorrections.o: SystematicCorrections.hh

MBUtils.o: MBUtils.hh

clean:
		@rm -f *.o *~  core
