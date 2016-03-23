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
  EvtRateHandler(int run, const std::string& inDir, bool ukdata=true) : runNumber(run),inputDir(inDir),UKdata(ukdata) {}
  virtual ~EvtRateHandler();
  int runNumber; //Run Number being read in
  std::string inputDir; //input data directory
  bool UKdata; //This is true if the tree format is UK style, False if it's official analyzer style
  double pol;  //Polarization of run as determined from number of events on each side
               
  double fiducialCut; //definition of a fiducial volume

  int polarization(int run); // determines the polarization. Assigns value to pol and returns polarization
                             // This comes from the database. Flipper on -> 1, Flipper Off -> -1

  double returnRunLength(int side) {return runLength[side];} // Return the length of the run (s)

  void CalcRates(double enBinWidth=10., double fidCut=45.); //evtType (0,1,2,3) for now. This returns a histogram of
                                                                            // rates
  std::vector< std::vector<double> > getRateVectors(int side);
  std::vector< std::vector<double> > getRateErrors(int side);
  TH1D getRateHist(int side, int evtType);

protected:
  virtual void dataReader(); //Read in data and fill histograms
 
  unsigned int numEnergyBins;
  double runLength[2]; // E/W
  std::vector <TH1D*> rateE;
  std::vector <TH1D*> rateW; // Rate histograms 
  std::vector< std::vector<double> > rateEvec;
  std::vector< std::vector<double> > rateWvec;
  std::vector< std::vector<double> > rateEerr; //Stores the statistical error for each Bin
  std::vector< std::vector<double> > rateWerr;
  double EmwpcX, EmwpcY, WmwpcX, WmwpcY, TimeE, TimeW, Erecon; //Branch Variables being read in
  float EmwpcX_f, EmwpcY_f, WmwpcX_f, WmwpcY_f, TimeE_f, TimeW_f, Erecon_f; // For reading in data from MPM replays
  int PID, Side, Type;
};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 
 
class SimEvtRateHandler: public EvtRateHandler {
public:
  SimEvtRateHandler(int run, const std::string& inDir): EvtRateHandler(run, inDir, true) {}

protected:
  void dataReader(); //Different set of variables for reverse calibrated simulated data
  double AsymWeight;
  std::vector < std::vector < double > >  mwpcPos; //holds the position of the event in the MWPC for simulated data
};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
class BGSubtractedRate {
public:
  BGSubtractedRate(int run, double enBin, double fidCut, bool ukdata=true, bool sim=false);
  ~BGSubtractedRate() {}

  int runNumber;
  double EnergyBinWidth; //Width of energy bins used
  double fiducialCut;
  //int evtType; //either 0, 1, or 23
  bool UKdata; //Whether the replay was done using UK code or MPM code
  bool Simulation; //Whether the rate is from simulation or not, in which case there is no background run

  std::vector<double> returnRunLengths(bool beta=true); // for BG do beta=false
  void calcBGSubtRates(); //Loads the rates and calculates BG subtr rates and errors
  std::vector<double> returnBGSubtRate(int side, int etype); //returns a vector holding all the BG subtracted rates.
  std::vector<double> returnBGSubtRateError(int side, int etype); //Returns background subtracted rate statistical error

  int getBackgroundRun(int run); //Returns the background run number for the run specified
  void CreateRateHistograms(); //Create, fill, and save rate histograms to file.
private:
 
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

  double runLengthBG[2]; // E/W
  double runLengthBeta[2];
  
  void LoadRatesByBin(); // Load rates and save them to vectors
  void CalcFinalRate();
  
  
};
#endif
