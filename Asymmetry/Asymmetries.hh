#ifndef ASYMMETRIES_HH
#define ASYMMETRIES_HH

#include <map>
#include <string>

#include "SQLinterface.hh"
#include "EvtRateHandler.hh"

// base class for loading octet information
class AsymmetryBase {
public:
  AsymmetryBase(int oct, double enBinWidth=10., double fidCut=45., bool ukdata=true, bool simulation=false);
  ~AsymmetryBase() {}

  void readOctetFile(); //populates the map runType
  int getRunsInOctet() {return runsInOctet;}
  bool isFullOctet(); //Checks whether the octet has all the necessary beta runs
  void loadRates(); //Loads the BG subtracted rates for each run and side
  bool checkIfBetaRun(std::string type); //Checks if run of type is a beta run or not
  void calcBGsubtractedEvts(); // Simply calculates and fills numEvtsByTypeByBin vector
  std::vector < double > getNumBGsubtrEvts(double enWinLow, double enWinHigh, int evtType); // returns the number of events summed over bin window
  void writeRatesToFile(); //For now this only uses type0 evts
  bool isAnaChoiceRateVectors() {return boolAnaChRtVecs;}
  std::vector < std::vector < std::vector<double> > > returnBGsubtractedRate(std::string runType); // returns A2, A5, etc below
  std::vector < std::vector < std::vector<double> > > returnBGsubtractedRateError(std::string runType);
  std::vector < double > returnRunLength(std::string runType); //Returns both beta and BG run length

  

protected:
  void makeAnalysisChoiceRateVectors(int anaChoice); //Makes new vectors based on the analysis choice (1-8)
  std::vector< std::vector<double> > numEvtsEastByTypeByBin; //number of events for each evt type, each bin summed 
  std::vector< std::vector<double> > numEvtsWestByTypeByBin; //number of events for each evt type, each bin summed 
  
  bool UKdata; //Boolean which is true for UK data and false if mpm
  bool Simulation; //Boolean to use simulated data
  int octet; // Holds the octet being analyzed
  
  double energyBinWidth;
  double fiducialCut;
  bool boolAnaChRtVecs;
  int runsInOctet;
  std::map <std::string,int> runType; // the key is the run type, mapped val is run number

  // The following vectors are bin by bin rates for each evt type on each side
  // for ecample, A2[type][side][bin] (A2[0-3][0-1][0-nBins])
  std::vector < std::vector < std::vector<double> > > A2, A5, A7, A10, B2, B5, B7, B10; //BG subtr rates for each beta run 
  std::vector < std::vector < std::vector<double> > > A2err, A5err, A7err, A10err, B2err, B5err, B7err, B10err; //BG subtr rate errors for each beta run

  //The following vectors store the length of the runs where A2len[0][0]=A2East length and A2[1][0] = A2East BG run length(A1)
  std::vector < std::vector < double > > A2len, A5len, A7len, A10len, B2len, B5len, B7len, B10len;

  // The following vectors sum over the appropriate event types for the analysisChoice [side][bin]
  std::vector < std::vector <double> > anaChoice_A2, anaChoice_A5, anaChoice_A7, anaChoice_A10, anaChoice_B2, anaChoice_B5, anaChoice_B7, anaChoice_B10; //BG subtr rates for each beta run 
  std::vector < std::vector <double> > anaChoice_A2err, anaChoice_A5err, anaChoice_A7err, anaChoice_A10err, anaChoice_B2err, anaChoice_B5err, anaChoice_B7err, anaChoice_B10err; //BG subtr rates for each beta run

  std::vector<double> binLowerEdge;
  std::vector<double> binUpperEdge; //Hold the Energy of the upper and lower bin edges
};


class OctetAsymmetry : public AsymmetryBase {
public:
  OctetAsymmetry(int oct, double enBinWidth=10., double fidCut=45., bool ukdata=true, bool simulation=false);
  ~OctetAsymmetry() {std::cout << "\n\n\n";}

  void makePlots();
  void calcAsymmetryBinByBin(int anaChoice=1); //Calculates the raw asymmetry bin by bin to be written to file and plotted
  void calcTotalAsymmetry(double enWinLow, double enWinHigh, int anaChoice=1); //Returns total raw asymmetry over energy window
  void calcSuperSum(int anaChoice=1); //Calculates the super sum over the entire octet for spectral comparisons
  void writeAsymToFile(int anaChoice);
  double returnTotalAsymmetry() {return totalAsymmetry;}
  double returnTotalAsymmetryError() {return totalAsymmetryError;}

  std::vector <double> returnSuperSum() {return superSum;}
  std::vector <double> returnSuperSumError() {return superSumError;}

private:
  
  std::vector <double> asymmetry; //Raw asymmetry in bins
  std::vector <double> asymmetryError; //Raw Asymmetry error in bins
  std::vector <double> superSum; //Raw asymmetry in bins
  std::vector <double> superSumError; //Raw Asymmetry error in bins
  double totalAsymmetry; //Bin summed asymmetry
  double totalAsymmetryError;
};
#endif
