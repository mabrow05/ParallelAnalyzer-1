# Makefile from Michael Brown

ObjSuf        = o
SrcSuf        = C
DllSuf        = so
OutPutOpt     = -o  
ROOTLIBS      = $(shell root-config --libs)
ROOTGLIBS     = $(shell root-config --glibs)
ROOTCFLAGS    = $(shell root-config --cflags)


# Work with Linux with egcs	
CXX           = g++ 
CXXFLAGS      = -O2 -Wall -fPIC -std=c++11 -Wl,--no-as-needed
LD            = g++
SOFLAGS       = -shared
LIBS          = $(ROOTLIBS) $(ROOTGLIBS)
CXXFLAGS     += $(ROOTCFLAGS) -I../include
LIBS         += -lSpectrum -lMinuit


.PHONY: all
all: Etrue_Edep

Etrue_Edep: Etrue_Edep.o
	$(CXX) $(CXXFLAGS) -o Etrue_Edep Etrue_Edep.o $(LIBS)
	@echo "Done"

clean: 
	@rm -f *.o *~  core Etrue_Edep
