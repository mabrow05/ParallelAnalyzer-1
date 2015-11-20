#ifndef ASYMMETRIES_HH
#define ASYMMETRIES_HH

#include <map>
#include <string>

#include "SQLinterface.hh"
#include "EvtRateHandler.hh"

// Class for determining asymmetry from one run
class AsymmetryBase {
public:
  AsymmetryBase(int oct) : UKdata(true), Simulation(false), octet(oct) {readOctetFile();}
  ~AsymmetryBase() {}

  void readOctetFile(); //populates the map runType
  int getRunsInOctet() {return runsInOctet;}
  bool UKdata;
  bool Simulation;
protected: 
  int octet;
  int runsInOctet;
  std::map <std::string,int> runType; // the key is the run type, mapped val is run number
};


class OctetAsymmetry : public AsymmetryBase {
public:
  OctetAsymmetry(int oct, double enBinWidth=10., double fidCut=45.);
  ~OctetAsymmetry() {std::cout << "\n\n\n";}

  void makePlots();
  void calcBGsubtractedEvts();
  double getNumBGsubtrEvts(double enWinLow, double enWinHigh, int evtType);
  bool checkIfBetaRun(std::string type);
  void writeRatesToFile(); //For now this only uses type0 evts

private:
  //double numBGsubtractedEvts; //double because this is calculated using the rate*runLength
  double energyBinWidth;
  double fiducialCut;
  std::vector< std::vector<double> > numEvtsByTypeByBin; //number of events for each evt type, each bin summed 
                                                         // over both sides
  std::vector<double> binLowerEdge;
  std::vector<double> binUpperEdge; //Hold the Energy of the upper and lower bin edges
};
#endif
