// Utility code holding all pertinent functions, variables, and structs for SimulationAnalyzer


#ifndef SIMPROCUTILS_HH
#define SIMPROCUTILS_HH

#include <vector>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <utility>

#include <TRandom3.h>
#include <TTree.h>
#include <TFile.h>
#include <TH1F.h>
#include <TF1.h>
#include <TH1D.h>
#include <TChain.h>
#include <TMath.h>


// This is the default error envelope as given by Michael M. (2008,2010) and Michael B. (2011,2012)

const std::vector < std::pair <Double_t,Double_t> > envelope2010 = {std::make_pair(0,2.5),
								    std::make_pair(200,200*0.0125),
								    std::make_pair(500,500*0.0125),
								    std::make_pair(1000,500*0.0125)};

const std::vector < std::pair <Double_t,Double_t> > envelope2011 = {std::make_pair(0,0.025*130.3),
								    std::make_pair(130.3,130.3*0.025),
								    std::make_pair(368.4938,0.018*368.4938),
								    std::make_pair(993.789,993.789*0.013),
								    std::make_pair(1000,1000.*0.013)};

const std::vector < std::pair <Double_t,Double_t> > envelope2012 = {std::make_pair(0,0.025*130.3),
								    std::make_pair(130.3,130.3*0.025),
								    std::make_pair(368.4938,0.018*368.4938),
								    std::make_pair(993.789,993.789*0.013),
								    std::make_pair(1000,1000.*0.013)};

const std::map <std::string, std::vector < std::pair<Double_t, Double_t> > > envelopes = {{"2010",envelope2010},
											  {"2011/2012",envelope2011},
											  {"2012/2013",envelope2012}};


// These are the nominal peak values to compare against for each source
const std::map <std::string,std::pair<Double_t,Double_t> > peaks2010 = {{"Sn",std::make_pair(365.629,365.394)},
									{"Ce",std::make_pair(125.173,123.539)},
									{"Bi1",std::make_pair(988.766,988.476)},
									{"Bi2",std::make_pair(492.094,497.781)}};

const std::map <std::string,std::pair<Double_t,Double_t> > peaks2011 = {{"Sn",std::make_pair(365.629,365.394)},
									{"Ce",std::make_pair(125.173,123.539)},
									{"Bi1",std::make_pair(0.,0.)},
									{"Bi2",std::make_pair(0.,0.)}};

const std::map <std::string,std::pair<Double_t,Double_t> > peaks2012 = {{"Sn",std::make_pair(365.629,365.394)},
									{"Ce",std::make_pair(0.,0.)},
									{"Bi1",std::make_pair(0.,0.)},
									{"Bi2",std::make_pair(0.,0.)}};



Int_t PID, type, side, primaryID; // basic analysis tags

Double_t initialMomentum[3]; //For reconstructing which detector type 1 events initially struck

Double_t Eprim, AsymWeight; // initial energy from simulation, and weight of event based on A*beta*cos(theta)

Double_t Erecon; // smeared, weighted, and trigger func corrected energy with the EQ2Etrue conversion applied



///////// Function declarations /////////////////

// This function holds all the meat
void revCalSimulation (std::string source, std::string geometry, UInt_t numEvents, bool linCorr, std::vector < std::vector <Double_t> > params, Int_t index);

//Function to return the x value of the max bin in a histogram 
Double_t GetXatMax(TH1D* hist, Double_t xmin=-1., Double_t xmax=-1.);

//Function to make initial check on state of parameter set.. i.e. directly comparing against error envelope
bool checkParamSetStatus(std::vector <Double_t> params,std::string source, std::string geometry); 

bool checkPeakStatus(std::vector < std::vector < Double_t> > meanAndSig, std::string source, std::string geometry); 

//Function to return a 2 element vector holding the fit mean and sigma of a peak
std::vector <Double_t> FitGaus(TH1D* histToFit, Double_t gausMean, Double_t min, Double_t max);

//Function to return the trigger function for each side in a std::vector in the form vec[side][param]
// where side==0 is East and side==1 is West
std::vector < std::vector < Double_t > > getTriggerFunctionParams(Int_t XeRunPeriod, Int_t nParams);

//Return the (nPE/keV) values for each PMT
std::vector < Double_t > getAlphaValues(Int_t fileNum);

//Function to return the probability of a trigger at a given energy given the parameters
Double_t triggerProbability(std::vector <Double_t> params, Double_t En);

//Get the conversion from EQ2Etrue
std::vector < std::vector < std::vector <double> > > getEQ2EtrueParams(std::string geometry);
  
//Function to return the correction to linearity offset
Double_t applyLinearityTwiddle (std::vector <Double_t> &params, Double_t EQ);

// Sets up the output tree for the simulated data
void SetUpOutputTree(TTree& tree);
  

//////////////////////// Data Structures /////////////////////////

struct Evis {
  double EvisE;
  double EvisW;
} evis;

struct Edep {
  double EdepE;
  double EdepW;
} edep;

struct EdepQ {
  double EdepQE;
  double EdepQW;
} edepQ;

struct time {
  double timeE;
  double timeW;
} Time;

struct MWPCEnergy {
  double MWPCEnergyE; 
  double MWPCEnergyW;
} mwpcE;

struct MWPCPos {
  double MWPCPosE[3];
  double MWPCPosW[3];
} mwpc_pos;

struct ScintPos {
  double ScintPosE[3];
  double ScintPosW[3];
} scint_pos;

struct ScintPosAdjusted {
  double ScintPosAdjE[3];
  double ScintPosAdjW[3];
} scint_pos_adj;

struct PMT_Evis {
  double Evis[8];
  double weight[8]; 
} pmt_Evis;


#endif
