
ObjSuf        = o
SrcSuf        = C
DllSuf        = so
OutPutOpt     = -o  
ROOTLIBS      = $(shell root-config --libs)
ROOTGLIBS     = $(shell root-config --glibs)
ROOTCFLAGS    = $(shell root-config --cflags)

VPATH = ./:../include:../Asymmetry

# Work with Linux with egcs	
CXX           = g++
CXXFLAGS      = -O2 -Wall -fPIC -std=c++11
LD            = g++
SOFLAGS       = -shared
LIBS          = $(ROOTLIBS) $(ROOTGLIBS)
CXXFLAGS     +=  $(ROOTCFLAGS) -I../include 
LIBS         += -lSpectrum -lMinuit

objects = SimulationProcessor.o

SimulationProcessor: $(objects)
	$(CXX) $(CXXFLAGS) -o $@ $(objects) $(LIBS)
	@echo "Done"

SimulationProcessor.o: SimulationProcessor.hh posMapReader.h


# -------------------------------------------------------------------------------
#  Generic compilation and linking step to make an executable from
#  a single *.cc file
#
#%: %.$(SrcSuf)
#	$(CXX) $(CXXFLAGS) -o $@ $< $(LIBS)
#	@echo "$@ done"

.PHONY : clean
clean:
		rm -f *.o *~ SimulationProcessor
