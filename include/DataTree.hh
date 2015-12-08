#ifndef DATATREE_H
#define DATATREE_H

#include <TTree.h>
#include <string>

struct MWPC {
  Double_t center;
  Double_t width;
  Double_t cathSum;
  Double_t maxValue;
  Int_t maxWire;
  Int_t mult;
  Int_t nClipped;
  Int_t err;
  Double_t rawCenter;
  Double_t height;
};

struct Scint {
  Double_t q1, q2, q3, q4;
  Double_t e1, de1, e2, de2, e3, de3, e4, de4;
  Double_t energy, denergy; 
  Double_t nPE1, nPE2, nPE3, nPE4;
};

class DataTree {
public:
  DataTree(); //Constructor
  ~DataTree();

  void makeOutputTree(std::string outputFile, std::string outputTree); 
  void fillOutputTree() {outputTree->Fill();}
  void writeOutputFile() {outputFile->Write();)
  void setupInputTree(std::string inputFile, std::string inputTree);
  void getEvent(UInt_t N) {inputTree->GetEvent(N);}
  Int_t getEntries() {return inputTree->GetEntriesFast();}

  //List of final variables
  Int_t TriggerNum; //trigger number from full.root file
  Int_t EvtN;       // Event number when read in 
  Int_t Sis00;      // Trigger source
  Double_t DeltaT;  // Time since last event
  Double_t Tof;     // time since last beam burst
  Double_t TimeE;   // blinded time East
  Double_t TimeW;   // blinded time West
  Double_t TDCE, TDCW; // Trigger TDC by side
  Double_t TDCE1, TDCE2, TDCE3, TDCE4; //East ind PMT trigger TDC
  Double_t TDCW1, TDCW2, TDCW3, TDCW4; //West ind PMT trigger TDC
  MWPC xE, xW, yE, yW; //MWPC structures
  Double_t Cathodes_Ex[16], Cathodes_Ey[16]; //array of cathode values
  Double_t Cathodes_Wx[16], Cathodes_Wy[16]; //array of cathode values
  Scint ScintE, ScintW; //Scintillator variable information
  Double_t EvisE, EvisW; // Visible (quenched) energy on each side
  Double_t CathSumE, CathSumW; //sum of both x-y cathode planes
  Double_t CathMaxE, CathMaxW; //smaller max-value of each plane
  Double_t EMWPC_E, EMWPC_W; // estimated wirechamber energy deposition
  Double_t AnodeE, AnodeW; // Wirechamber anode signals
  bool PassedAnoE, PassedAnoW; //Whether or not passed anode threshold
  bool PassedCathE, PassedCathW; //Whether or not passed cathode threshold
  bool TaggedBackE, TaggedBackW; // Whether or not tagged by backing veto
  bool TaggedTopE, TaggedTopW;  // Whether tagged as muon by top veto
                                // which was only on East side so west is 0
  bool TaggedDriftE, TaggedDriftW; // Whether tagged by muon veto drift tubes
  Double_t EastBackADC, WestBackADC; // muon backing veto ADC
  Double_t EastBackTDC, WestBackTDC; // muon backing veto TDC
  Double_t EastDriftVetoADC, WestDriftVetoADC; // muon drift veto ADC
  Double_t EastTopVetoADC; //East Top veto ADC (No west)
  Double_t EastTopVetoTDC; // East top veto TDC (No west)
  bool EvnbGood, BkhfGood; // DAQ header and footer quality
  
  Int_t PID; //Particle ID
  Int_t Type; //Event type
  Int_t Side; //Earlier Trigger side
  Double_t ProbIII; //Probability of type 3 event
  Double_t Erecon; //Final reconstructed energy of an event

private:
  TFile *inputFile, *outputFile; 
  TTree *inputTree, *outputTree; 
  //std::string inputTreeName, outputTreeName; 
}
  



#endif


