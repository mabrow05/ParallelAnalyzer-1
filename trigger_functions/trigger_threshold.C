#include <iostream>
#include <vector>
#include <cstdlib>

void trigger_threshold(Int_t XeRunPeriod) {

  bool mpmData=true;
  Char_t temp[500];
  //Create file to write parameters of fit.
  if (mpmData) sprintf(temp,"%s/trigger_functions_XePeriod_MPM_%i.dat",getenv("TRIGGER_FUNC"),XeRunPeriod);
  else sprintf(temp,"%s/trigger_functions_XePeriod_%i.dat",getenv("TRIGGER_FUNC"),XeRunPeriod);
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

  TChain *chain;
  if (mpmData) chain = new TChain("phys");
  else chain = new TChain("pass4");

  if (mpmData) {
    for (Int_t i=0; i<numRuns; i++) {
      sprintf(temp,"%s/hists/spec_%i.root",getenv("UCNAOUTPUTDIR"),XeRuns[i]);
      chain->AddFile(temp);
    }
  }
  else {
    for (Int_t i=0; i<numRuns; i++) {
      sprintf(temp,"%s/replay_pass4_%i.root",getenv("REPLAY_PASS4"),XeRuns[i]);
      chain->AddFile(temp);
    }
  }
 
  //sprintf(temp,"/extern/UCNA/replay_pass4_MB/replay_pass4_%i.root",runNumber);
  //TFile *file = new TFile(temp,"READ");

  //TTree *t = (TTree*)file->Get("pass4");

  Int_t nbins = 75;
  Double_t East_upper_limit = 100.;
  Double_t West_upper_limit = 100.;
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
  //TF1 *erf = new TF1("erf","([0]+[1]*TMath::Erf((x-[2])/[3]))+[4]*TMath::Gaus(x,[5],[6])",0.,150.);
  TF1 *erf = new TF1("erf","([0]+[1]*TMath::Erf((x-[2])/[3])*TMath::TanH([7]*x)+[4]*TMath::Gaus(x,[5],[6]))",0.,150.);
  //TF1 *erf = new TF1("erf","([0]+[1]*TMath::TanH([2]*x-[3])+[4]*TMath::Gaus(x,[5],[6]))",0.,150.);
  erf->SetParameter(0,0.5); //Constant offset of erf
  erf->SetParameter(1,0.5); //Scaling of erf
  erf->SetParameter(2,30.); //Mean of gaussian integrated for erf
  erf->SetParameter(3,12.); //std. dev. of gaussian integrated for erf
  erf->SetParameter(4,-0.1); //Scaling of additional gaussian for cancellation of erf overshooting distribution
  erf->SetParameter(5,40.); //Mean of additional gaussian
  //erf->SetParLimits(5,35.,40.);
  erf->SetParameter(6,12.); //Std. dev of additional gaussian
  //erf->SetParLimits(6,9.,13.);
  erf->SetParameter(7,35.);
  //erf->SetParLimits(7,30.,45.);

  //Working Parameters of the east trigger function for erf*tanh+gaus
  /*erf->SetParameter(0,0.5); //Constant offset of erf
  erf->SetParameter(1,0.5); //Scaling of erf
  erf->SetParameter(2,30.); //Mean of gaussian integrated for erf
  erf->SetParameter(3,12.); //std. dev. of gaussian integrated for erf
  erf->SetParameter(4,-0.1); //Scaling of additional gaussian for cancellation of erf overshooting distribution
  erf->SetParameter(5,40.); //Mean of additional gaussian
  erf->SetParameter(6,12.); //Std. dev of additional gaussian
  erf->SetParameter(7,35.);*/
  //erf->SetParameter(8,7.);
//erf->SetParLimits(6,1.,25.);
  //erf->SetParLimits(4,-5.,0.);

  TCanvas *c1 = new TCanvas("c1"," ",1200.,1600.);
  c1->Divide(2,2);
  c1->cd(1);
  chain->Draw("EvisE>>Etype1","Type==1 && Side==1 && PID==1 && EvisE>0.");
  c1->cd(2);  
  chain->Draw("EvisE>>Etype23","Type==2 && Side==1 && PID==1 && EvisE>0.");
  c1->cd(3);
  Etotal->Add(Etype1,Etype23);
  Etotal->Draw();
  c1->cd(4);
  Etrigg->Divide(Etype1,Etotal);
  Etrigg->SetStats(0);
  //Etrigg->Draw("P");

  Etrigg->Fit("erf","","",0.,East_upper_limit);
  
  Etrigg->Draw("P");
  triggFunc << erf->GetParameter(0) << " " << erf->GetParameter(1) << " " << erf->GetParameter(2) << " " << erf->GetParameter(3) << " " << erf->GetParameter(4) << " " << erf->GetParameter(5) << " " <<  endl; //erf->GetParameter(6) << endl;


  TCanvas *c2 = new TCanvas("c2"," ",1200.,1600.);
  c2->Divide(2,2);
  c2->cd(1);
  chain->Draw("EvisW>>Wtype1","Type==1 && Side==0 && PID==1 && EvisW>0.");
  c2->cd(2);  
  chain->Draw("EvisW>>Wtype23","Type==2 && Side==0 && PID==1 && EvisW>0.");
  c2->cd(3);
  Wtotal->Add(Wtype1,Wtype23);
  Wtotal->Draw();
  c2->cd(4);
  Wtrigg->Divide(Wtype1,Wtotal);
  Wtrigg->SetStats(0);
  //Wtrigg->Draw("P");

  erf->SetParameter(0,0.5); //Constant offset of erf
  erf->SetParameter(1,0.5); //Scaling of erf
  erf->SetParameter(2,25.); //Mean of gaussian integrated for erf
  erf->SetParameter(3,13.); //std. dev. of gaussian integrated for erf
  erf->SetParameter(4,-0.1); //Scaling of additional gaussian for cancellation of erf overshooting distribution
  erf->SetParameter(5,33.); //Mean of additional gaussian
  erf->SetParameter(6,19.); //Std. dev of additional gaussian
  erf->SetParameter(7,25.);
  //erf->SetParameter(7,0.036);
  //erf->SetParameter(8,25.);

  //erf->SetParameter(0,0.5); //Constant offset of erf
  //erf->SetParameter(1,0.5); //Scaling of erf
  //erf->SetParameter(2,10.); //Mean of gaussian integrated for erf
  //erf->SetParameter(3,5.); //std. dev. of gaussian integrated for erf
  //erf->SetParameter(4,-0.1); //Scaling of additional gaussian for cancellation of erf overshooting distribution
  //erf->SetParameter(5,14.); //Mean of additional gaussian
  //erf->SetParameter(6,6.); //Std. dev of additional gaussian
  //erf->SetParameter(7,25.);
  //erf->SetParameter(7,0.036);
  //erf->SetParameter(8,25.);


  Wtrigg->Fit("erf","","",0.,West_upper_limit);
  //Wtrigg->Fit("erf","R");

  Wtrigg->Draw("P");

  triggFunc << erf->GetParameter(0) << " " << erf->GetParameter(1) << " " << erf->GetParameter(2) << " " << erf->GetParameter(3) << " " << erf->GetParameter(4) << " " << erf->GetParameter(5) << " " << endl;// erf->GetParameter(6);

  triggFunc.close();
  
  //Write the canvases to file
  sprintf(temp,"%i",XeRunPeriod);
  TString canvas1 = TString(getenv("TRIGGER_FUNC")) + TString("trigger_functions_XePeriod_")+TString(temp)+TString(".pdf(");
  TString canvas2 = TString(getenv("TRIGGER_FUNC")) + TString("trigger_functions_XePeriod_")+TString(temp)+TString(".pdf)");
  c1->Print(canvas1);
  c2->Print(canvas2);
  
}
