#include <iostream>
#include <vector>
#include <cstdlib>

void trigger_threshold(Int_t XeRunPeriod) {

  Char_t temp[500];
  //Create file to write parameters of fit.
  sprintf(temp,"%s/trigger_functions_XePeriod_%i.dat",getenv("TRIGGER_FUNC"),XeRunPeriod);
  ofstream triggFunc(temp);

  // Read in the Xe runs in this runPeriod
  vector <Int_t> XeRuns;
  sprintf(temp,"../run_lists/Xenon_Calibration_Run_Period_%i.dat",XeRunPeriod);
  ifstream XeRunList(temp);

  Int_t run, numRuns=0;
  while (XeRunList >> run) {
    XeRuns.push_back(run);
    numRuns++;
  }

  if (numRuns!=XeRuns.size()) {
    cout << "Runs!=XeRuns.size()\n";
    exit(0);
  }
  XeRunList.close();
  cout << "There are " << numRuns << " in this Xe run period\n";

  TChain *chain = new TChain("pass4");

  for (Int_t i=0; i<numRuns; i++) {
    sprintf(temp,"/extern/UCNA/replay_pass4_MB/replay_pass4_%i.root",XeRuns[i]);
    chain->AddFile(temp);
  }
 
  //sprintf(temp,"/extern/UCNA/replay_pass4_MB/replay_pass4_%i.root",runNumber);
  //TFile *file = new TFile(temp,"READ");

  //TTree *t = (TTree*)file->Get("pass4");

  Int_t nbins = 75;
  Double_t East_upper_limit = 100.;
  Double_t West_upper_limit = 60.;
  TH1F *Etype1 = new TH1F("Etype1","West Type 1: EvisE",nbins,0.,East_upper_limit);
  TH1F *Etype23 = new TH1F("Etype23","West Type 2/3: EvisE",nbins,0.,East_upper_limit);
  TH1F *Etotal = new TH1F("Etotal","West Type 1,2/3: EvisE",nbins,0.,East_upper_limit);
  TH1F *Etrigg = new TH1F("Etrigg","East Trigger Probability",nbins,0.,East_upper_limit);
  Etrigg->SetMarkerStyle(20);

  TH1F *Wtype1 = new TH1F("Wtype1","East Type 1: EvisW",nbins,0.,West_upper_limit);
  TH1F *Wtype23 = new TH1F("Wtype23","East Type 2/3: EvisW",nbins,0.,West_upper_limit);
  TH1F *Wtotal = new TH1F("Wtotal","East Type 1,2/3: EvisW",nbins,0.,West_upper_limit);
  TH1F *Wtrigg = new TH1F("Wtrigg","West Trigger Probability",nbins,0.,West_upper_limit);
  Wtrigg->SetMarkerStyle(20);

  //TF1 *erf = new TF1("erf","([5]*TMath::Erf((x-[0])/[1])+0.5)+[2]*TMath::Gaus(x,[3],[4])",0.,150.);
  TF1 *erf = new TF1("erf","([0]+[1]*TMath::Erf((x-[2])/[3]))+[4]*TMath::Gaus(x,[5],[6])",0.,150.);
  
  erf->SetParameter(0,0.5); //Constant offset of erf
  erf->SetParameter(1,0.5); //Scaling of erf
  erf->SetParameter(2,25.); //Mean of gaussian integrated for erf
  erf->SetParameter(3,13.); //std. dev. of gaussian integrated for erf
  erf->SetParameter(4,-0.1); //Scaling of additional gaussian for cancellation of erf overshooting distribution
  erf->SetParameter(5,33.); //Mean of additional gaussian
  erf->SetParameter(6,19.); //Std. dev of additional gaussian
  //erf->SetParLimits(6,1.,25.);
  //erf->SetParLimits(4,-5.,0.);

  TCanvas *c1 = new TCanvas("c1"," ",1200.,1600.);
  c1->Divide(2,2);
  c1->cd(1);
  chain->Draw("EvisE>>Etype1","type_pass4==1 && side_pass4==1 && PID_pass4==1 && EvisE>0.");
  c1->cd(2);  
  chain->Draw("EvisE>>Etype23","type_pass4==2 && side_pass4==1 && PID_pass4==1 && EvisE>0.");
  c1->cd(3);
  Etotal->Add(Etype1,Etype23);
  Etotal->Draw();
  c1->cd(4);
  Etrigg->Divide(Etype1,Etotal);
  Etrigg->SetStats(0);
  //Etrigg->Draw("P");

  Etrigg->Fit("erf","","",0.,East_upper_limit);
  
  Etrigg->Draw("P");
  triggFunc << erf->GetParameter(0) << " " << erf->GetParameter(1) << " " << erf->GetParameter(2) << " " << erf->GetParameter(3) << " " << erf->GetParameter(4) << endl;


  TCanvas *c2 = new TCanvas("c2"," ",1200.,1600.);
  c2->Divide(2,2);
  c2->cd(1);
  chain->Draw("EvisW>>Wtype1","type_pass4==1 && side_pass4==0 && PID_pass4==1 && EvisW>0.");
  c2->cd(2);  
  chain->Draw("EvisW>>Wtype23","type_pass4==2 && side_pass4==0 && PID_pass4==1 && EvisW>0.");
  c2->cd(3);
  Wtotal->Add(Wtype1,Wtype23);
  Wtotal->Draw();
  c2->cd(4);
  Wtrigg->Divide(Wtype1,Wtotal);
  Wtrigg->SetStats(0);
  //Wtrigg->Draw("P");

  erf->SetParameter(0,0.5); //Constant offset of erf
  erf->SetParameter(1,0.5); //Scaling of erf
  erf->SetParameter(2,10.); //Mean of gaussian integrated for erf
  erf->SetParameter(3,5.); //std. dev. of gaussian integrated for erf
  erf->SetParameter(4,-0.1); //Scaling of additional gaussian for cancellation of erf overshooting distribution
  erf->SetParameter(5,14.); //Mean of additional gaussian
  erf->SetParameter(6,6.); //Std. dev of additional gaussian
  
  Wtrigg->Fit("erf","","",0.,West_upper_limit);
  //Wtrigg->Fit("erf","R");

  Wtrigg->Draw("P");

  triggFunc << erf->GetParameter(0) << " " << erf->GetParameter(1) << " " << erf->GetParameter(2) << " " << erf->GetParameter(3) << " " << erf->GetParameter(4);

  triggFunc.close();
  
  //Write the canvases to file
  sprintf(temp,"%i",XeRunPeriod);
  TString canvas1 = TString(getenv("TRIGGER_FUNC")) + TString("trigger_functions_XePeriod_")+TString(temp)+TString(".pdf(");
  TString canvas2 = TString(getenv("TRIGGER_FUNC")) + TString("trigger_functions_XePeriod_")+TString(temp)+TString(".pdf)");
  c1->Print(canvas1);
  c1->Print(canvas2);
  
}
