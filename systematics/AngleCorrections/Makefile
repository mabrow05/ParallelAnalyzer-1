# Makefile from Jianglai

ObjSuf        = o
SrcSuf        = C
DllSuf        = so
OutPutOpt     = -o  
ROOTLIBS      = $(shell root-config --libs)
ROOTGLIBS     = $(shell root-config --glibs)
ROOTCFLAGS    = $(shell root-config --cflags)

VPATH = ./:../../include/

# Work with Linux with egcs	
CXX           = g++ 
CXXFLAGS      = -O2 -Wall -fPIC -g -std=c++11
LD            = g++
SOFLAGS       = -shared
LIBS          = $(ROOTLIBS) $(ROOTGLIBS)
CXXFLAGS     += $(ROOTCFLAGS) -I../../include 
LIBS         += -lSpectrum -lMinuit

objects = MBUtils.o calibrationTools.o BetaSpectrum.o

AngleCorrections.exe: AngleCorrections.o $(objects)
	$(CXX) $(CXXFLAGS) -o AngleCorrections.exe AngleCorrections.o $(objects) $(LIBS)
	@echo "Done"

processedDeltaAngle.exe: processedDeltaAngle.o $(objects)
	$(CXX) $(CXXFLAGS) -o processedDeltaAngle.exe processedDeltaAngle.o $(objects) $(LIBS)
	@echo "Done"

AsymmCalculator.exe: AsymmCalculator.o $(objects)
	$(CXX) $(CXXFLAGS) -o AsymmCalculator.exe AsymmCalculator.o $(objects) $(LIBS)
	@echo "Done"


MBUtils.o: MBUtils.hh

BetaSpectrum.o: BetaSpectrum.hh

positionMapHandler.o: positionMapHandler.hh MBUtils.o

TriggerMap.o: TriggerMap.hh MBUtils.o

calibrationTools.o: calibrationTools.hh

MWPCPositionResponse.o: MWPCPositionResponse.hh



# -------------------------------------------------------------------------------
#  Generic compilation and linking step to make an executable from
#  a single *.cc file
#
#%: %.$(SrcSuf)
#	$(CXX) $(CXXFLAGS) -o $@ $< $(LIBS)
#	@echo "$@ done"

clean:
		@rm -f *.o *~  core
