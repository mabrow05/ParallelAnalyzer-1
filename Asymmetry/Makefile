ObjSuf		= o
SrcSuf		= C
DllSuf		= so
OutPutOpt	= -o
ROOTLIBS	= $(shell root-config --libs)
ROOTGLIBS	= $(shell root-config --glibs)
ROOTCFLAGS	= $(shell root-config --cflags)


CXX		= g++
CXXFLAGS	= -O2 -Wall -fPIC
LD		= g++
SOFLAGS		= -shared
LIBS		= $(ROOTLIBS) $(ROOTGLIBS)
CXXFLAGS       += $(ROOTCFLAGS) -Iinclude
LIBS	       += -lSpectrum -lMinuit

objects = Asymmetries.o SQLinterface.o SystematicCorrections.o EvtRateHandler.o

.PHONY: all
all: MBAnalyzer.exe 

MBAnalyzer.exe: MBAnalyzer.o $(objects) 
	$(CXX) $(CXXFLAGS) -o MBAnalyzer.exe MBAnalyzer.o $(objects) $(LIBS) 
	@echo "Done"

Asymmetries.o: Asymmetries.hh

EvtRateHandler.o: EvtRateHandler.hh SQLinterface.hh

SQLintercace.o: SQLinterface.hh

SystematicCorrections.o: SystematicCorrections.hh

clean:
		@rm -f *.o *~  core
