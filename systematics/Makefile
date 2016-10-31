# Makefile from Jianglai

ObjSuf        = o
SrcSuf        = C
DllSuf        = so
OutPutOpt     = -o  
ROOTLIBS      = $(shell root-config --libs)
ROOTGLIBS     = $(shell root-config --glibs)
ROOTCFLAGS    = $(shell root-config --cflags)

VPATH = ./:../include/

# Work with Linux with egcs	
CXX           = g++ 
CXXFLAGS      = -O2 -Wall -fPIC -std=c++11 -g
LD            = g++
SOFLAGS       = -shared
LIBS          = $(ROOTLIBS) $(ROOTGLIBS)
CXXFLAGS     += $(ROOTCFLAGS) -I../include
LIBS         += -lSpectrum -lMinuit

DeltaExpProcessor.exe: DeltaExpProcessor.o positionMapHandler.o SQLinterface.o calibrationTools.o MBUtils.o TriggerMap.o BetaSpectrum.o
	$(CXX) $(CXXFLAGS) -o DeltaExpProcessor.exe DeltaExpProcessor.o MBUtils.o calibrationTools.o positionMapHandler.o TriggerMap.o BetaSpectrum.o ../Asymmetry/SQLinterface.o $(LIBS)
	@echo "Done"

XuanStyle_DeltaExpProcessor.exe: XuanStyle_DeltaExpProcessor.o positionMapHandler.o SQLinterface.o calibrationTools.o MBUtils.o TriggerMap.o BetaSpectrum.o
	$(CXX) $(CXXFLAGS) -o XuanStyle_DeltaExpProcessor.exe XuanStyle_DeltaExpProcessor.o MBUtils.o calibrationTools.o positionMapHandler.o TriggerMap.o BetaSpectrum.o ../Asymmetry/SQLinterface.o $(LIBS)
	@echo "Done"

MBUtils.o: MBUtils.hh

BetaSpectrum.o: BetaSpectrum.hh

positionMapHandler.o: positionMapHandler.hh MBUtils.o

TriggerMap.o: TriggerMap.hh MBUtils.o

SQLinterface.o: ../Asymmetry/SQLinterface.hh

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
