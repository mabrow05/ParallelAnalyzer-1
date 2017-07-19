#include <TChain.h>
#include <TH1D.h>
#include <TFile.h>
#include <TTree.h>
#include <TMath.h>

#include <vector>
#include <iostream>

int main() {

  Double_t eastEff = 0.9988;
  Double_t westEff = 1.;

  std::vector <int> runs {19899,19900,19902,19904,19924,19925,19927,19929};

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
    std::cout << eastThresh << "\tEff: " << eff << std::endl;
  }
  westThresh-=0.1; eff=1.;  
  while (westEff<eff) {
    westThresh+=0.01;
    total = chain->GetEntries("side==0 && type==0 && PID==1");
    trigg = chain->GetEntries(TString::Format("side==0 && type==0 && PID==1 && MWPCEnergyE>%f",westThresh));
    eff = trigg/total;  
    std::cout << eastThresh << "\tEff: " << eff << std::endl;
  }
  westThresh-=0.01; eff=1.;
  while (westEff<eff) {
    westThresh+=0.001;
    total = chain->GetEntries("side==0 && type==0 && PID==1");
    trigg = chain->GetEntries(TString::Format("side==0 && type==0 && PID==1 && MWPCEnergyE>%f",westThresh));
    eff = trigg/total;  
    std::cout << eastThresh << "\tEff: " << eff << std::endl;
  }

  std::cout << "East Threshold = " << eastThresh << "\n";
  std::cout << "West Threshold = " << westThresh << "\n";
}
