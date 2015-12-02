#ifndef MAKEOUTPUTTREE_H
#define MAKEOUTPUTTREE_H

#include <TTree.h>

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

};

class OutputTree {
public:
  OutputTree(std::string treeName); //Constructor
  ~OutputTree();

  TTree makeOutputTree(); 

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
  

}
  



#endif


