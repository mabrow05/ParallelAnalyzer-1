# Makefile from Jianglai

ObjSuf        = o
SrcSuf        = C
DllSuf        = so
OutPutOpt     = -o  
ROOTLIBS      = $(shell root-config --libs)
ROOTGLIBS     = $(shell root-config --glibs)
ROOTCFLAGS    = $(shell root-config --cflags)

VPATH = ./:../include/:../Asymmetry/
#:../Asymmetry/

# Work with Linux with egcs	
CXX           = g++ 
CXXFLAGS      = -O2 -Wall -fPIC -std=c++11
LD            = g++
SOFLAGS       = -shared
LIBS          = $(ROOTLIBS) $(ROOTGLIBS)
CXXFLAGS     += $(ROOTCFLAGS) -I../include
LIBS         += -lSpectrum -lMinuit

revCalSim.exe: revCalSim.o positionMapHandler.o SQLinterface.o calibrationTools.o MBUtils.o TriggerMap.o
	$(CXX) $(CXXFLAGS) -o revCalSim.exe revCalSim.o MBUtils.o calibrationTools.o positionMapHandler.o TriggerMap.o SQLinterface.o $(LIBS)
	@echo "Done"

MBUtils.o: MBUtils.hh

positionMapHandler.o: positionMapHandler.hh MBUtils.o

TriggerMap.o: TriggerMap.hh MBUtils.o

SQLinterface.o: SQLinterface.hh

calibrationTools.o: calibrationTools.hh

#posMapReader.o: posMapReader.h



# -------------------------------------------------------------------------------
#  Generic compilation and linking step to make an executable from
#  a single *.cc file
#
#%: %.$(SrcSuf)
#	$(CXX) $(CXXFLAGS) -o $@ $< $(LIBS)
#	@echo "$@ done"

clean:
		@rm -f *.o *~  core
