#include <TChain.h>
#include <TH1D.h>
#include <TFile.h>
#include <TTree.h>
#include <TMath.h>

#include <vector>
#include <iostream>

int main() {


  std::vector <int> runs {19899,19900,19902,19904,19924,19925,19927,19929};
  std::vector <int> bgruns {19901,19903,19926,19928,19930};

  Double_t bgtime=0.,fgtime=0.;
  Double_t bgTriggE=0,bgNoTriggE=0,fgTriggE=0,fgNoTriggE=0;
  Double_t bgTriggW=0,bgNoTriggW=0,fgTriggW=0,fgNoTriggW=0;
  Double_t bgTriggEerr=0.,bgNoTriggEerr=0.,fgTriggEerr=0.,fgNoTriggEerr=0.;
  Double_t bgTriggWerr=0.,bgNoTriggWerr=0.,fgTriggWerr=0.,fgNoTriggWerr=0.;

  for (auto run:runs) {
    Double_t time;
    TFile *f = new TFile(TString::Format("%s/replay_pass3_%i.root",getenv("REPLAY_PASS3"),run),"READ");
    TTree *t = (TTree*)f->Get("pass3");
    t->SetBranchAddress("Time",&time);
    t->GetEvent(t->GetEntriesFast()-1);
    fgtime+=time;

    fgTriggE += (Double_t)t->GetEntries("Type==0 && PID==1 && Side==0 && Erecon>320. && Erecon<400.");
    fgTriggW += (Double_t)t->GetEntries("Type==0 && PID==1 && Side==1 && Erecon>320. && Erecon<400.");
    fgNoTriggE += (Double_t)t->GetEntries("PID==0 && Side==0 && Erecon>320. && Erecon<400. && Sis00==1");
    fgNoTriggW += (Double_t)t->GetEntries("PID==0 && Side==1 && Erecon>320. && Erecon<400. && Sis00==2");
  }
  
  for (auto run:bgruns) {
    Double_t time;
    TFile *f = new TFile(TString::Format("%s/replay_pass3_%i.root",getenv("REPLAY_PASS3"),run),"READ");
    TTree *t = (TTree*)f->Get("pass3");
    t->SetBranchAddress("Time",&time);
    t->GetEvent(t->GetEntriesFast()-1);
    bgtime+=time;

    bgTriggE += (Double_t)t->GetEntries("Type==0 && PID==1 && Side==0 && Erecon>320. && Erecon<400.");
    bgTriggW += (Double_t)t->GetEntries("Type==0 && PID==1 && Side==1 && Erecon>320. && Erecon<400.");
    bgNoTriggE += (Double_t)t->GetEntries("PID==0 && Side==0 && Erecon>320. && Erecon<400. && Sis00==1");
    bgNoTriggW += (Double_t)t->GetEntries("PID==0 && Side==1 && Erecon>320. && Erecon<400. && Sis00==2");
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

  
  /*
      TChain *fgchain = new TChain("pass3");
  for (auto run : runs) {    
    fgchain->Add(TString::Format("%s/replay_pass3_%i.root",getenv("REPLAY_PASS3"),run));
  }

  TChain *bgchain = new TChain("pass3");
  for (auto run : bgruns) {
    bgchain->Add(TString::Format("%s/replay_pass3_%i.root",getenv("REPLAY_PASS3"),run));
  }
  

  TH1D *triggE = new TH1D("triggE","East wirechamber trigger",100,280.,440.);
  TH1D *notriggE = new TH1D("notriggE","East no wirechamber trigger",100,280.,440.);

  TH1D *triggW = new TH1D("triggW","West wirechamber trigger",100,280.,440.);
  TH1D *notriggW = new TH1D("notriggW","West no wirechamber trigger",100,280.,440.);

  TCanvas *c1 = new TCanvas("c1","c1",1000,800);
  c1->Divide(2,2);
  
  c1->cd(1);
  chain->Draw("Erecon>>triggE","Sis00==1 && PID==1 && Side==0 && Erecon>320. && Erecon<400.");
  c1->cd(2);
  chain->Draw("Erecon>>notriggE","Sis00==1 && PID==0 && Side==0 && Erecon>320. && Erecon<400.");
  c1->cd(3);
  chain->Draw("Erecon>>triggW","Sis00==2 && Type==0 && PID==1 && Side==1 && Erecon>320. && Erecon<400.");
  c1->cd(4);
  chain->Draw("Erecon>>notriggW","Sis00==2 && PID==0 && Side==1 && Erecon>320. && Erecon<400.");
  */

}
