//////////////////////////////////////////////////////////
// Takes replay_pass4_#.root files and creates a spec_#.root
// file to be used in other analyzers. The spec file doesn't
// include all the values as of yet... 
//////////////////////////////////////////////////////////

#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>

// ROOT libraries
#include "TRandom3.h"
#include <TTree.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TF1.h>
#include <TH1D.h>

#include "fullTreeVariables.h"
#include "MWPCGeometry.h"
#include "pedestals.h"
#include "cuts.h"
#include "basic_reconstruction.h"

#include "replay_pass2.h"
#include "replay_pass3.h"
#include "replay_pass4.h"

using namespace std;

struct mwpc {
  float center;
  float width;
  float cathSum;
  float maxValue;
  int maxWire;
  int mult;
  int nClipped;
  int err;
  float rawCenter;
  float height;
};

int main(int argc, char *argv[])
{
  cout.setf(ios::fixed, ios::floatfield);
  cout.precision(12);

  cout << "Run " << argv[1] << " ..." << endl;
  cout << "... Creating spec_" << argv[1] << ".root from replay_pass4 data ..." << endl;


  // Open output ntuple
  char tempOut[500];
  sprintf(tempOut, "%s/spec_%s.root",getenv("UK_SPEC_REPLAY"), argv[1]);
  TFile *fileOut = new TFile(tempOut,"RECREATE");
  TTree *Tout = new TTree("phys", "phys");

  // Variables
  Int_t EvtN=0, TriggerNum=0;
  mwpc xEmpm={}, yEmpm={}, xWmpm={}, yWmpm={}; //if i ever get the chance to add in other pieces of michaels position stuff these are to be used
  Float_t UBtime_BB_float, AnodeE_float, AnodeW_float, EvisE_float, EvisW_float, EvisTot_float, timeE_float, timeW_float;

  Tout->Branch("EvtN",&EvtN,"EvtN/I");
  Tout->Branch("TriggerNum",&TriggerNum,"TriggerNum/I");
  Tout->Branch("Tof",&UBtime_BB_float,"Tof/F");
  
  Tout->Branch("xEmpm",&xEmpm,"center/F:width:cathSum:maxValue:maxWire/I:mult:nClipped:err:rawCenter/I:height");
  Tout->Branch("yEmpm",&yEmpm,"center/F:width:cathSum:maxValue:maxWire/I:mult:nClipped:err:rawCenter/I:height");
  Tout->Branch("xWmpm",&xWmpm,"center/F:width:cathSum:maxValue:maxWire/I:mult:nClipped:err:rawCenter/I:height");
  Tout->Branch("yWmpm",&yWmpm,"center/F:width:cathSum:maxValue:maxWire/I:mult:nClipped:err:rawCenter/I:height");

  Tout->Branch("AnodeE", &AnodeE_float, "AnodeE/F");
  Tout->Branch("EvisE", &EvisE_float, "EvisE/F");
  Tout->Branch("AnodeW", &AnodeW_float, "AnodeW/F"); 
  Tout->Branch("EvisW", &EvisW_float, "EvisW/F");
  Tout->Branch("Erecon", &EvisTot_float, "Erecon/F");

  Tout->Branch("TimeE", &timeE_float, "TimeE/F");
  Tout->Branch("TimeW", &timeW_float, "TimeW/F");

  Tout->Branch("PID",  &PID_pass4,  "PID/I");
  Tout->Branch("Type", &type_pass4, "Type/I");
  Tout->Branch("Side", &side_pass4, "Side/I");
  //Tout->Branch("posError_pass2", &posError_pass2, "posError_pass2/I");
  

  // Open input ntuple
  char tempIn[500];
  sprintf(tempIn, "%s/replay_pass4_%s.root", getenv("REPLAY_PASS4"),argv[1]);

  TFile *fileIn = new TFile(tempIn, "READ");
  TTree *Tin = (TTree*)(fileIn->Get("pass4"));

  Tin->SetBranchAddress("AnodeE", &AnodeE);
  Tin->SetBranchAddress("AnodeW", &AnodeW);

  Tin->SetBranchAddress("timeE", &timeE);
  Tin->SetBranchAddress("timeW", &timeW);
  Tin->SetBranchAddress("timeE_BB", &timeE_BB);
  Tin->SetBranchAddress("timeW_BB", &timeW_BB);
  Tin->SetBranchAddress("UBtime", &UBtime);
  Tin->SetBranchAddress("UBtime_BB", &UBtime_BB);
  Tin->SetBranchAddress("twoFoldE", &twoFoldE);
  Tin->SetBranchAddress("twoFoldW", &twoFoldW);

  //SetBranchAddress with leaves to store the Evis and weights for all 8 PMTs
  Tin->SetBranchAddress("pmt_Evis", &pmt_Evis);

  Tin->SetBranchAddress("xE_pass4", &xE_pass4);
  Tin->SetBranchAddress("yE_pass4", &yE_pass4);
  Tin->SetBranchAddress("xW_pass4", &xW_pass4);
  Tin->SetBranchAddress("yW_pass4", &yW_pass4);

  Tin->SetBranchAddress("EvisTot", &EvisTot);
  Tin->SetBranchAddress("EvisE", &EvisE);
  Tin->SetBranchAddress("EvisW", &EvisW);
  
  //Tin->SetBranchAddress("EreconTot", &EreconTot, "EreconTot/D");
  //Tin->SetBranchAddress("EreconE", &EreconE, "EreconE/D");
  //Tin->SetBranchAddress("EreconW", &EreconW, "EreconW/D");

  Tin->SetBranchAddress("PID_pass4",  &PID_pass4);
  Tin->SetBranchAddress("type_pass4", &type_pass4);
  Tin->SetBranchAddress("side_pass4", &side_pass4);
  Tin->SetBranchAddress("posError_pass4", &posError_pass4);

  
  int nEvents = Tin->GetEntries();
  cout << "... Processing nEvents = " << nEvents << endl;

  // Loop over events
  for (int i=0; i<nEvents; i++) {
    Tin->GetEvent(i);

    if (type_pass4==-1) {
      type_pass4=4;
	}

    //Position Variables
    xEmpm.center = (float)xE_pass4;
    yEmpm.center = (float)yE_pass4;
    xWmpm.center = (float)xW_pass4;
    yWmpm.center = (float)yW_pass4;
    
    timeE_float = (float)timeE;
    timeW_float = (float)timeW;
    AnodeE_float = (float)AnodeE;
    AnodeW_float = (float)AnodeW;
    EvisE_float = (float)EvisE;
    EvisW_float = (float)EvisW;
    EvisTot_float = (float)EvisTot;
    UBtime_BB_float = (float)UBtime_BB;
    
    

    Tout->Fill();
  }

  // Write output ntuple
  fileOut->Write();
  fileOut->Close();

  return 0;
}
