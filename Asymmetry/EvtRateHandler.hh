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
#include "SQLinterface.hh"



class EvtRateHandler {
public:
  EvtRateHandler(int run, bool ukdata , const std::string& inDir) : runNumber(run),inputDir(inDir),UKdata(ukdata) {}
  ~EvtRateHandler() {if (rateE) delete rateE; if (rateW) delete rateW; } //May not want to delete these pointers...
  int runNumber; //Run Number being read in
  std::string inputDir; //input data directory
  bool UKdata; //This is true if the tree format is UK style, False if it's official analyzer style
  double pol;  //Polarization of run as determined from number of events on each side
               
  double fiducialCut; //definition of a fiducial volume

  int polarization(int run); // determines the polarization. Assigns value to pol and returns polarization
                             // This comes from the database. Flipper on -> 1, Flipper Off -> -1

  void CalcRates(int evtType, double enBins, double fidCut=45.); //evtType (0,1,23) for now. This returns a histogram of
                                                                            // rates
  std::vector<double> getRateVector(int side);
  TH1D getRateHist(int side);

protected:
  virtual void dataReader(int evtType); //Read in data and fill histogram
 
  int numEnergyBins;
  TH1D *rateE, *rateW; // Rate histogram 
  std::vector<double> rateEvec;
  std::vector<double> rateWvec;
  double EmwpcX, EmwpcY, WmwpcX, WmwpcY, TimeE, TimeW, Erecon; //Branch Variables being read in
  int PID, Side, Type;
};
  
class SimEvtRateHandler: public EvtRateHandler {
public:
  SimEvtRateHandler(int run, const std::string& inDir): EvtRateHandler(run, true, inDir) {}
  ~SimEvtRateHandler() {}
protected:
  virtual void dataReader(int evtType); //Different set of variables for reverse calibrated simulated data
};

class BGSubtractedRate {
public:
  BGSubtractedRate(int run, double enBin, double fidCut, int etype=0, bool sim=false): runNumber(run), EnergyBinWidth(enBin), fiducialCut(fidCut), evtType(etype), Simulation(sim) {}
  ~BGSubtractedRate();

  int runNumber;
  double EnergyBinWidth; //Width of energy bins used
  double fiducialCut;
  int evtType; //either 0, 1, or 23
  bool Simulation; //Whether the rate is from simulation or not, in which case there is no background run
  

  std::vector<double> ReturnBGSubtRate(int side); //This handles Loading rates by bin and creating the rate histograms via
                                          // the private methods below. It then returns a vector holding all the BG subtracted rates.
  int getBackgroundRun(int run); //Returns the background run number for the run specified
  void CreateRateHistograms(); //Create, fill, and save rate histograms to file.
private:
 
  std::vector<double> BetaRate; //Save the rates here for beta run
  std::vector<double> BGRate;   // Save the background rates here
  std::vector<double> FinalRate; //This is the difference in the rates
  
  void LoadRatesByBin(int side); // Load rates and save them to vectors
  void CalcFinalRate();
  
  
};
#endif
