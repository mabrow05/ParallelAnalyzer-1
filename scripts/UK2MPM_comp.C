/*
  Script to compare UK and MPM data which has been run through the calibrations
  from each respective analyzer

  run using... root -l 'UK2MPM_comp.C(int runNumber)'
 */

#include <string>
#include "../include/MButils.hh"
  
void UK2MPM_comp(int runNumber) {

  gStyle->SetOptStat(0);

  //Filepaths
  std::string mpmData = std::string(getenv("UCNAOUTPUTDIR"))+"/hists/spec_" + itos(runNumber) + ".root";
  std::string ukData = std::string(getenv("REPLAY_PASS4"))+"/replay_pass4_" + itos(runNumber) + ".root";

  //Open files
  TFile *mpmFile = new TFile(mpmData.c_str(),"READ");
  TFile *ukFile = new TFile(ukData.c_str(),"READ");

  //Get Trees
  TTree *mpmTree = (TTree*)mpmFile->Get("phys");
  TTree *ukTree = (TTree*)ukFile->Get("pass4");

  //Making histograms for comparison  
  TH1F *EvisE0_mpm = new TH1F("EvisE0_mpm","EvisE Type 0", 100, 0., 800.);
  TH1F *EvisW0_mpm = new TH1F("EvisW0_mpm","EvisW Type 0", 100, 0., 800.);
  TH1F *EvisE0_uk = new TH1F("EvisE0_uk","EvisE Type 0", 100, 0., 800.);
  TH1F *EvisW0_uk = new TH1F("EvisW0_uk","EvisW Type 0", 100, 0., 800.);

  TH1F *EvisE1_mpm = new TH1F("EvisE1_mpm","EvisE Type 1", 200, 0., 800.);
  TH1F *EvisW1_mpm = new TH1F("EvisW1_mpm","EvisW Type 1", 200, 0., 800.);
  TH1F *EvisE1_uk = new TH1F("EvisE1_uk","EvisE Type 1", 200, 0., 800.);
  TH1F *EvisW1_uk = new TH1F("EvisW1_uk","EvisW Type 1", 200, 0., 800.);

  TH1F *EvisE23_mpm = new TH1F("EvisE23_mpm","EvisE Type 2/3", 200, 0., 800.);
  TH1F *EvisW23_mpm = new TH1F("EvisW23_mpm","EvisW Type 2/3", 200, 0., 800.);
  TH1F *EvisE23_uk = new TH1F("EvisE23_uk","EvisE Type 2/3", 200, 0., 800.);
  TH1F *EvisW23_uk = new TH1F("EvisW23_uk","EvisW Type 2/3", 200, 0., 800.); 

  TH1I *Type_mpm = new TH1I("Type_mpm","Event Types", 4, 0, 4);
  TH1I *Type_uk = new TH1I("Type_uk","Event Types", 4, 0, 4);
  
  EvisE0_mpm->SetLineColor(kRed);
  EvisW0_mpm->SetLineColor(kRed);
  EvisE0_uk->SetLineColor(kBlue);
  EvisW0_uk->SetLineColor(kBlue);

  EvisE1_mpm->SetLineColor(kRed);
  EvisW1_mpm->SetLineColor(kRed);
  EvisE1_uk->SetLineColor(kBlue);
  EvisW1_uk->SetLineColor(kBlue);

  EvisE23_mpm->SetLineColor(kRed);
  EvisW23_mpm->SetLineColor(kRed);
  EvisE23_uk->SetLineColor(kBlue);
  EvisW23_uk->SetLineColor(kBlue);

  Type_mpm->SetLineColor(kRed);
  Type_uk->SetLineColor(kBlue);



  

  //Plot what is of interest

  TCanvas *cEvis0 = new TCanvas("cEvis"," ", 1400., 600.);
  cEvis0->Divide(2,1);
  cEvis0->cd(1);
  
  mpmTree->Draw("EvisE>>EvisE0_mpm", "Type==0 && Side==0 && PID==1");
  ukTree->Draw("EvisE>>EvisE0_uk", "Type==0 && Side==0 && PID==1 && EvisE>0.1", "SAME");

  cEvis0->cd(2);
  
  mpmTree->Draw("EvisW>>EvisW0_mpm", "Type==0 && Side==1 && PID==1");
  ukTree->Draw("EvisW>>EvisW0_uk", "Type==0 && Side==1 && PID==1 && EvisW>0.1", "SAME");

  TCanvas *cEvis1 = new TCanvas("cEvis1"," ", 1400., 600.);
  cEvis1->Divide(2,1);
  cEvis1->cd(1);
  
  mpmTree->Draw("EvisE>>EvisE1_mpm", "Type==1 && Side==0 && PID==1");
  ukTree->Draw("EvisE>>EvisE1_uk", "Type==1 && Side==0 && PID==1 && EvisE>0.1", "SAME");

  cEvis1->cd(2);
  
  mpmTree->Draw("EvisW>>EvisW1_mpm", "Type==1 && Side==1 && PID==1");
  ukTree->Draw("EvisW>>EvisW1_uk", "Type==1 && Side==1 && PID==1 && EvisW>0.1", "SAME");

  TCanvas *cEvis23 = new TCanvas("cEvis23"," ", 1400., 600.);
  cEvis23->Divide(2,1);
  cEvis23->cd(1);
  
  mpmTree->Draw("EvisE>>EvisE23_mpm", "(Type==2 || Type==3)&& Side==0 && PID==1");
  ukTree->Draw("EvisE>>EvisE23_uk", "(Type==2 || Type==3) && Side==0 && PID==1 && EvisE>0.1", "SAME");

  cEvis23->cd(2);
  
  mpmTree->Draw("EvisW>>EvisW23_mpm", "(Type==2 || Type==3) && Side==1 && PID==1");
  ukTree->Draw("EvisW>>EvisW23_uk", "(Type==2 || Type==3) && Side==1 && PID==1 && EvisW>0.1", "SAME");


  TCanvas *cTypes = new TCanvas("cTypes"," ", 700., 600.);
  
  mpmTree->Draw("Type>>Type_mpm", "PID==1");
  ukTree->Draw("Type>>Type_uk", "PID==1", "SAME");
}
