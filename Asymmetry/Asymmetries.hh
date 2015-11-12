#ifndef ASYMMETRIES_HH
#define ASYMMETRIES_HH

#include <map>

#include "SQLinterface.hh"
#inlcude "EvtRateHandler.hh"

// Class for determining asymmetry from one run
class AsymmetryBase {
public:
  AsymmetryBase(int oct) : octet(oct) {}
  ~AsymmetryBase() {}

  void readOctetFile(); //populates the map runType
protected: 
  int octet;
  std::map<string,int> runType; // the key is the run type, mapped val is run number
}
