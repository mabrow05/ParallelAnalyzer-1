#include <TChain.h>
#include <TH1D.h>
#include <TFile.h>
#include <TTree.h>
#include <TMath.h>

#include <vector>
#include <iostream>

int main() {

  Double_t enlow = 323.;
  Double_t enhigh = 403.;

  std::vector <int> runs {19924,19925,19927,19929};
  //{19899,19900,19902,19904,19924,19925,19927,19929}; // Sn
  //{18144,18147,18149,18152,18156,18159,18161,18164};// Octet 16
  std::vector <int> bgruns {19926,19928,19930};
  //{19901,19903,19926,19928,19930};
  //{18143,18146,18151,18154,18155,18158,18163,18166};

  Double_t bgtime=0.,fgtime=0.;
  Double_t bgTriggE=0,bgNoTriggE=0,fgTriggE=0,fgNoTriggE=0;
  Double_t bgTriggW=0,bgNoTriggW=0,fgTriggW=0,fgNoTriggW=0;
  Double_t bgTriggEerr=0.,bgNoTriggEerr=0.,fgTriggEerr=0.,fgNoTriggEerr=0.;
  Double_t bgTriggWerr=0.,bgNoTriggWerr=0.,fgTriggWerr=0.,fgNoTriggWerr=0.;

  //  TString eastTriggCuts = "Type==0 && PID==1 && Side==0 && Erecon>320. && Erecon<400.";
  //TString westTriggCuts = TString::Format("(ScintW.q1+ScintW.q2+ScintW.q3+ScintW.q4)>100. && PassedCathW && Erecon>%f && Erecon<%f",enlow,enhigh);

  TString eastTriggCuts = TString::Format("(Sis00==1 || Sis00==33) && PassedCathE && Erecon>%f && Erecon<%f",enlow,enhigh);
  TString westTriggCuts = TString::Format("(Sis00==2 || Sis00==34) && PassedCathW && Erecon>%f && Erecon<%f",enlow,enhigh);

  TString eastNoTriggCuts = TString::Format("(Sis00==1) && !PassedCathE && Erecon>%f && Erecon<%f",enlow,enhigh);
  TString westNoTriggCuts = TString::Format("(Sis00==2) && !PassedCathW && Erecon>%f && Erecon<%f",enlow,enhigh);
  
  //TString eastNoTriggCuts = TString::Format("(ScintE.q1+ScintE.q2+ScintE.q3+ScintE.q4)>50. && !PassedCathE && Sis00<4 && Erecon>%f && Erecon<%f",enlow,enhigh);
  //TString westNoTriggCuts = TString::Format("(ScintW.q1+ScintW.q2+ScintW.q3+ScintW.q4)>100. && !PassedCathW && Sis00<4 && Erecon>%f && Erecon<%f",enlow,enhigh);
    

  for (auto run:runs) {
    Double_t time;
    TFile *f = new TFile(TString::Format("%s/replay_pass3_%i.root",getenv("REPLAY_PASS3"),run),"READ");
    TTree *t = (TTree*)f->Get("pass3");
    t->SetBranchAddress("Time",&time);
    t->GetEvent(t->GetEntriesFast()-1);
    fgtime+=time;

    fgTriggE += (Double_t)t->GetEntries(eastTriggCuts);
    fgTriggW += (Double_t)t->GetEntries(westTriggCuts);
    fgNoTriggE += (Double_t)t->GetEntries(eastNoTriggCuts);
    fgNoTriggW += (Double_t)t->GetEntries(westNoTriggCuts);

    delete f;
  }
  
  for (auto run:bgruns) {
    Double_t time;
    TFile *f = new TFile(TString::Format("%s/replay_pass3_%i.root",getenv("REPLAY_PASS3"),run),"READ");
    TTree *t = (TTree*)f->Get("pass3");
    t->SetBranchAddress("Time",&time);
    t->GetEvent(t->GetEntriesFast()-1);
    bgtime+=time;

    bgTriggE += (Double_t)t->GetEntries(eastTriggCuts);
    bgTriggW += (Double_t)t->GetEntries(westTriggCuts);
    bgNoTriggE += (Double_t)t->GetEntries(eastNoTriggCuts);
    bgNoTriggW += (Double_t)t->GetEntries(westNoTriggCuts);

    delete f;
  }

  fgTriggEerr = TMath::Sqrt(fgTriggE);
  fgTriggWerr = TMath::Sqrt(fgTriggW);
  fgNoTriggEerr = TMath::Sqrt(fgNoTriggE);
  fgNoTriggWerr = TMath::Sqrt(fgNoTriggW);
  bgTriggEerr = TMath::Sqrt(bgTriggE);
  bgTriggWerr = TMath::Sqrt(bgTriggW);
  bgNoTriggEerr = TMath::Sqrt(bgNoTriggE);
  bgNoTriggWerr = TMath::Sqrt(bgNoTriggW);

  fgTriggE/=fgtime;
  fgTriggW/=fgtime;
  fgNoTriggE/=fgtime;
  fgNoTriggW/=fgtime;

  bgTriggE/=bgtime;
  bgTriggW/=bgtime;
  bgNoTriggE/=bgtime;
  bgNoTriggW/=bgtime;

  fgTriggEerr/=fgtime;
  fgTriggWerr/=fgtime;
  fgNoTriggEerr/=fgtime;
  fgNoTriggWerr/=fgtime;

  bgTriggEerr/=bgtime;
  bgTriggWerr/=bgtime;
  bgNoTriggEerr/=bgtime;
  bgNoTriggWerr/=bgtime;
  
  fgTriggE-=bgTriggE;
  fgTriggW-=bgTriggW;
  fgNoTriggE-=bgNoTriggE;
  fgNoTriggW-=bgNoTriggW;

  fgTriggEerr = TMath::Sqrt(TMath::Power(fgTriggEerr,2.)+TMath::Power(bgTriggEerr,2.));
  fgTriggWerr = TMath::Sqrt(TMath::Power(fgTriggWerr,2.)+TMath::Power(bgTriggWerr,2.));
  fgNoTriggEerr = TMath::Sqrt(TMath::Power(fgNoTriggEerr,2.)+TMath::Power(bgNoTriggEerr,2.));
  fgNoTriggWerr = TMath::Sqrt(TMath::Power(fgNoTriggWerr,2.)+TMath::Power(bgNoTriggWerr,2.));

  Double_t EastTotal = fgTriggE+fgNoTriggE;
  Double_t EastTotalErr = TMath::Sqrt(TMath::Power(fgTriggEerr,2.)+TMath::Power(fgNoTriggEerr,2.));
  Double_t EastEff = fgTriggE/EastTotal;
  Double_t EastEffErr = TMath::Abs(fgNoTriggEerr*EastEff/EastTotal);

  Double_t WestTotal = fgTriggW+fgNoTriggW;
  Double_t WestTotalErr = TMath::Sqrt(TMath::Power(fgTriggWerr,2.)+TMath::Power(fgNoTriggWerr,2.));
  Double_t WestEff = fgTriggW/WestTotal;
  Double_t WestEffErr = TMath::Abs(fgNoTriggWerr*WestEff/WestTotal);

  std::cout << "East trigg = " << fgTriggE << " +/- " << fgTriggEerr << std::endl;
  std::cout << "East no trigg = " << fgNoTriggE << " +/- " << fgNoTriggEerr << std::endl;
  std::cout << "East efficiency = " << EastEff << " +/- " << EastEffErr << std::endl;
  std::cout << "West trigg = " << fgTriggW << " +/- " << fgTriggWerr << std::endl;
  std::cout << "West no trigg = " << fgNoTriggW << " +/- " << fgNoTriggWerr << std::endl;
  std::cout << "West efficiency = " << WestEff << " +/- " << WestEffErr << std::endl;

  
 
}
