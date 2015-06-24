#include <TTree.h>
#include <TBranch.h>
#include <TLeaf.h>
#include <TClonesArray.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH1D.h>
#include <TF1.h>
#include <TH1.h>
#include <TApplication.h>
#include <TChain.h>
#include <vector>
#include <cstdlib>

using namespace std;

class Sim2Data
{
public:
  Sim2Data(int run, double fidCut);
  ~Sim2Data();
  int run;
  double pol; //polarization of run, determined from data
  float fiducialCut; // Fiducial cut used for analysis
  int type0evts; // used for normalization
  double type0integral; // used for normalization
  TFile *outfile;

  vector <TH1D*> eastData; //Types refer to index (2->2/3)
  vector <TH1D*> westData;
  vector <TH1D*> eastSim;
  vector <TH1D*> westSim;

  void histoCreator(); //Create histograms for both data and sim
  void dataReader(); //Read in data and fill data histograms
  void reverseCalSims(); //Read in reverse calibrated data, apply physics weight and normalize
  
};
