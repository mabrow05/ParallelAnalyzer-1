# Makefile from Jianglai

ObjSuf        = o
SrcSuf        = C
DllSuf        = so
OutPutOpt     = -o  
ROOTLIBS      = $(shell root-config --libs)
ROOTGLIBS     = $(shell root-config --glibs)
ROOTCFLAGS    = $(shell root-config --cflags)


# Work with Linux with egcs	
CXX           = g++ 
CXXFLAGS      = -O2 -Wall -fPIC
LD            = g++
SOFLAGS       = -shared
LIBS          = $(ROOTLIBS) $(ROOTGLIBS)
CXXFLAGS     += $(ROOTCFLAGS) -I../include
LIBS         += -lSpectrum -lMinuit

.PHONY: all
all: source_peaks.exe source_peaks_EvisPMTbyPMT.exe sim_source_peaks_EvisPMTbyPMT.exe sim_source_peaks_EnergyPeak.exe source_peaks_EnergyPeak.exe 
source_peaks.exe: source_peaks.o
	$(CXX) $(CXXFLAGS) -o source_peaks.exe source_peaks.o $(LIBS)
	@echo "Done"

source_peaks_EvisPMTbyPMT.exe: source_peaks_EvisPMTbyPMT.o
	$(CXX) $(CXXFLAGS) -o source_peaks_EvisPMTbyPMT.exe source_peaks_EvisPMTbyPMT.o $(LIBS)
	@echo "Done"

sim_source_peaks_EvisPMTbyPMT.exe: sim_source_peaks_EvisPMTbyPMT.o
	$(CXX) $(CXXFLAGS) -o sim_source_peaks_EvisPMTbyPMT.exe sim_source_peaks_EvisPMTbyPMT.o $(LIBS)
	@echo "Done"

sim_source_peaks_EnergyPeak.exe: sim_source_peaks_EnergyPeak.o
	$(CXX) $(CXXFLAGS) -o sim_source_peaks_EnergyPeak.exe sim_source_peaks_EnergyPeak.o $(LIBS)
	@echo "Done"

source_peaks_EnergyPeak.exe: source_peaks_EnergyPeak.o
	$(CXX) $(CXXFLAGS) -o source_peaks_EnergyPeak.exe source_peaks_EnergyPeak.o $(LIBS)
	@echo "Done"

# -------------------------------------------------------------------------------
#  Generic compilation and linking step to make an executable from
#  a single *.cc file
#
#%: %.$(SrcSuf)
#	$(CXX) $(CXXFLAGS) -o $@ $< $(LIBS)
#	@echo "$@ done"

clean:
		@rm -f *.o *~  core
