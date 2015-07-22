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
CXXFLAGS       += $(ROOTCFLAGS) -I../include
LIBS	       += -lSpectrum -lMinuit

.PHONY: all
all: MWPC_fitSourcePeaksSIM.exe MWPC_fitSourcePeaks.exe

MWPC_fitSourcePeaksSIM.exe: MWPC_fitSourcePeaksSIM.o
	$(CXX) $(CXXFLAGS) -o MWPC_fitSourcePeaksSIM.exe MWPC_fitSourcePeaksSIM.o $(LIBS) 
	@echo "Done"

MWPC_fitSourcePeaks.exe: MWPC_fitSourcePeaks.o
	$(CXX) $(CXXFLAGS) -o MWPC_fitSourcePeaks.exe MWPC_fitSourcePeaks.o $(LIBS) 
	@echo "Done"

clean:
		@rm -f *.o *~  core
