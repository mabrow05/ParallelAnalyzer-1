#ifndef EVTRATEHANDLER_HH
#define EVTRATEHANDLER_HH

#include <TTree.h>
#include <TBranch.h>
#include <TLeaf.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH1D.h>
#include <TF1.h>
#include <TH1.h>
#include <vector>
#include <cstdlib>
#include <cmath>
#include "SQLinterface.hh"



class EvtRateHandler {
public:
  EvtRateHandler(int run, const std::string& inDir, double enBinWidth=10., double fidCut=100., bool ukdata=true, bool unblind=false, bool sep23=false);
  virtual ~EvtRateHandler();
  int runNumber; //Run Number being read in
  std::string inputDir; //input data directory
  double fiducialCut; //definition of a fiducial volume
  bool UKdata; //This is true if the tree format is UK style, False if it's official analyzer style
  double pol;  //Polarization of run as determined from number of events on each side
               
 

  int polarization(int run); // determines the polarization. Assigns value to pol and returns polarization
                             // This comes from the database. Flipper on -> 1, Flipper Off -> -1

  double returnRunLength(int side) {return runLength[side];} // Return the length of the run (s)

  void CalcRates(); //evtType (0,1,2,3) for now. This returns a histogram of
                                                                            // rates
  std::vector< std::vector<double> > getRateVectors(int side) { 
    return side==0 ? rateEvec : ( side==1 ? rateWvec : throw "BAD SIDE GIVEN TO getRateVectors"); 
  } ;

  std::vector< std::vector<double> > getRateErrors(int side) {
    return side==0 ? rateEerr : ( side==1 ? rateWerr : throw "BAD SIDE GIVEN TO getRateErrors"); 
  } ;

  TH1D getCountHist(int side, int evtType); //Return the raw counts histogram

protected:
  virtual void dataReader(); //Read in data and fill histograms
  
  bool unblinded; //Holds whether the result should be unblinded
  bool separate23; //Whether or not to separate the type 2/3 events
  unsigned int numEnergyBins;
  double runLength[2]; // E/W
  double UCNMonIntegral;
  std::vector <TH1D*> hisE;
  std::vector <TH1D*> hisW; // Total Counts histograms 
  std::vector< std::vector<double> > rateEvec;
  std::vector< std::vector<double> > rateWvec;
  std::vector< std::vector<double> > rateEerr; //Stores the statistical error for each Bin
  std::vector< std::vector<double> > rateWerr;
  
};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 
 
class SimEvtRateHandler: public EvtRateHandler {
public:
  SimEvtRateHandler(int run, const std::string& inDir, double enBinWidth=10., double fidCut=100., bool applyAsymm=true, bool unblind=false): EvtRateHandler(run, inDir, enBinWidth, fidCut, true, unblind),applyAsymmetry(applyAsymm) {}

protected:
  void dataReader(); //Different set of variables for reverse calibrated simulated data
  bool applyAsymmetry;
  
};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
class BGSubtractedRate {
public:
  BGSubtractedRate(int run, double enBin, double fidCut, bool ukdata=true, bool sim=false, bool applyAsym=true, bool unblind=false);
  ~BGSubtractedRate() {}

  std::vector<double> returnRunLengths(bool beta=true); // for BG do beta=false
  void calcBGSubtRates(); //Loads the rates and calculates BG subtr rates and errors
  std::vector<double> returnBGSubtRate(int side, int etype); //returns a vector holding all the BG subtracted rates.
  std::vector<double> returnBGSubtRateError(int side, int etype); //Returns background subtracted rate statistical error

  int getBackgroundRun(int run); //Returns the background run number for the run specified
  void CreateRateHistograms(); //Create, fill, and save rate histograms to file.

private:
 
  int runNumber;
  double EnergyBinWidth; //Width of energy bins used
  double fiducialCut;
  //int evtType; //either 0, 1, or 23
  bool UKdata; //Whether the replay was done using UK code or MPM code
  bool Simulation; //Whether the rate is from simulation or not, in which case there is no background run
  bool applyAsymmetry; //if rate is from simulation, whether or not to apply the event weight from A **LEGACY! SHOULD ALWAYS BE FALSE
  bool UNBLIND; //SHOULD BE FALSE UNTIL UNBLINDING

  std::vector< std::vector<double> > BetaRateE; //Save the rates here for beta run
  std::vector< std::vector<double> > BetaRateErrorE; //Save the error here for beta run
  std::vector< std::vector<double> > BGRateE;   // Save the background rates here
  std::vector< std::vector<double> > BGRateErrorE;   // Save the background rates here
  std::vector< std::vector<double> > FinalRateE; //This is the difference in the rates
  std::vector< std::vector<double> > FinalRateErrorE; //This is the statistical error in the difference in the rates
  

  std::vector< std::vector<double> > BetaRateW; //Save the rates here for beta run
  std::vector< std::vector<double> > BetaRateErrorW; //Save the error here for beta run
  std::vector< std::vector<double> > BGRateW;   // Save the background rates here
  std::vector< std::vector<double> > BGRateErrorW;   // Save the background rates here
  std::vector< std::vector<double> > FinalRateW; //This is the difference in the rates
  std::vector< std::vector<double> > FinalRateErrorW; //This is the statistical error in the difference in the rates

  std::vector <double> runLengthBG; // E/W
  std::vector <double> runLengthBeta;
  
  void LoadRatesByBin(); // Load rates and save them to vectors
  void CalcFinalRate();
  
  
};
#endif
