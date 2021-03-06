#include <TChain.h>
#include <TH1D.h>
#include <TFile.h>
#include <TTree.h>
#include <TMath.h>

#include <vector>
#include <iostream>

int main() {

  /*East trigg = 25.7736 +/- 0.0356907
    East no trigg = 0.0427383 +/- 0.00702683
    East efficiency = 0.998345 +/- 0.000271735
    West trigg = 25.2393 +/- 0.0353165
    West no trigg = 0.0360786 +/- 0.00640295
    West efficiency = 0.998573 +/- 0.000252966
  */
  Double_t eastEff = 0.99912;
  Double_t westEff = 0.99974;

  std::vector <int> runs {19899,19900,19902,19904,19924,19925,19927,19929};//{18144,18147,18149,18152,18156,18159,18161,18164};
  //{18144,18147,18149,18152,18156,18159,18161,18164};
  //{19899,19900,19902,19904,19924,19925,19927,19929};//

  TChain *chain = new TChain("revCalSim");

  for (auto run:runs) {
    chain->Add(TString::Format("%s/sources/revCalSim_%i_Sn113.root",getenv("REVCALSIM"),run));
    //chain->Add(TString::Format("%s/beta_highStatistics/revCalSim_%i_Beta.root",getenv("REVCALSIM"),run));
  }
  
  Double_t eastThresh = 0.5;
  Double_t westThresh = 0.5;  
  Double_t eff = 1.;

  Double_t trigg, total;

  while (eastEff<eff) {
    eastThresh+=0.1;
    total = chain->GetEntries("side==0 && type==0 && PID==1");
    trigg = chain->GetEntries(TString::Format("side==0 && type==0 && PID==1 && MWPCEnergyE>%f",eastThresh));
    eff = trigg/total;
    std::cout << eastThresh << "\tEff: " << eff << std::endl;
  }
  eastThresh-=0.1; eff=1.; 
  while (eastEff<eff) {
    eastThresh+=0.01;
    total = chain->GetEntries("side==0 && type==0 && PID==1");
    trigg = chain->GetEntries(TString::Format("side==0 && type==0 && PID==1 && MWPCEnergyE>%f",eastThresh));
    eff = trigg/total;  
    std::cout << eastThresh << "\tEff: " << eff << std::endl;
  }
  eastThresh-=0.01; eff=1.;
  while (eastEff<eff) {
    eastThresh+=0.001;
    total = chain->GetEntries("side==0 && type==0 && PID==1");
    trigg = chain->GetEntries(TString::Format("side==0 && type==0 && PID==1 && MWPCEnergyE>%f",eastThresh));
    eff = trigg/total;  
    std::cout << eastThresh << "\tEff: " << eff << std::endl;
  }

  //////////////////WEST///////////////////////
  eff=1.;

  while (westEff<eff) {
    westThresh+=0.1;
    total = chain->GetEntries("side==0 && type==0 && PID==1");
    trigg = chain->GetEntries(TString::Format("side==0 && type==0 && PID==1 && MWPCEnergyE>%f",westThresh));
    eff = trigg/total;   
    std::cout << westThresh << "\tEff: " << eff << std::endl;
  }
  westThresh-=0.1; eff=1.;  
  while (westEff<eff) {
    westThresh+=0.01;
    total = chain->GetEntries("side==0 && type==0 && PID==1");
    trigg = chain->GetEntries(TString::Format("side==0 && type==0 && PID==1 && MWPCEnergyE>%f",westThresh));
    eff = trigg/total;  
    std::cout << westThresh << "\tEff: " << eff << std::endl;
  }
  westThresh-=0.01; eff=1.;
  while (westEff<eff) {
    westThresh+=0.001;
    total = chain->GetEntries("side==0 && type==0 && PID==1");
    trigg = chain->GetEntries(TString::Format("side==0 && type==0 && PID==1 && MWPCEnergyE>%f",westThresh));
    eff = trigg/total;  
    std::cout << westThresh << "\tEff: " << eff << std::endl;
  }

  std::cout << "East Threshold = " << eastThresh << "\n";
  std::cout << "West Threshold = " << westThresh << "\n";
}
