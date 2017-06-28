
#include <TString.h>
#include <vector>
#include <TMath.h>
#include <TFile.h>
#include <TH1D.h>
#include <TTree.h>

#include <iostream>
#include <cstdlib>
#include <fstream>
#include <vector>
#include <iomanip>

// Function to return beta when given the kinetic energy of an electron
Double_t returnBeta(Double_t En) { 
  Double_t me = 510.998928; //rest mass energy of electron in keV
  return sqrt(En*En+2.*En*me)/(En+me);
};


void FieldAsymmetry(TString field,TString year) {


  std::vector <Double_t> rE_polE (100,0.);
  std::vector <Double_t> rE_polE_err (100,0.);  
  std::vector <Double_t> rW_polE (100,0.);
  std::vector <Double_t> rW_polE_err (100,0.);  
  std::vector <Double_t> rE_polW (100,0.);
  std::vector <Double_t> rE_polW_err (100,0.);  
  std::vector <Double_t> rW_polW (100,0.);
  std::vector <Double_t> rW_polW_err (100,0.);  

  TString basepath = TString::Format("/extern/mabrow05/ucna/geant4work/output/bigSim/%sField_%s/",field.Data(),year.Data());

  TFile *f;
  
  TH1D *hisE_polE = new TH1D("hisE_polE","hisE",100,0.,1000.);
  TH1D *hisW_polE = new TH1D("hisW_polE","hisW",100,0.,1000.);

  TH1D *hisE_polW = new TH1D("hisE_polE","hisE",100,0.,1000.);
  TH1D *hisW_polW = new TH1D("hisW_polE","hisW",100,0.,1000.);


  Double_t EdepQ[2], MWPCEnergy[2],primKE,primTheta;

  //Start with East polarization
  for (int i=0; i<5; ++i) {
    std::cout << "In Folder " << i << std::endl;
    for (int j=0; j<1000; ++j) {
      f = new TFile(TString::Format("%s/%s%i/analyzed_%i.root",
					   basepath.Data(),
					   (field==TString("good")?"G":
					    (field==TString("bad")?"B":"F")),i,j),"READ");
      TTree *t = (TTree*)f->Get("anaTree");
      t->SetBranchAddress("EdepQ",EdepQ);
      t->SetBranchAddress("MWPCEnergy",MWPCEnergy);
      t->SetBranchAddress("primKE",&primKE);
      t->SetBranchAddress("primTheta",&primTheta);
      
  
      for (UInt_t i=0; i<t->GetEntriesFast() ; ++i) {

	t->GetEntry(i);

	if ( EdepQ[0]>0. && EdepQ[1]==0. && MWPCEnergy[0]>0. && MWPCEnergy[1]==0. )
	  hisE_polE->Fill(primKE,1+0.12*sqrt(primKE*primKE+2.*primKE*510.998928)/(primKE+510.998928)*TMath::Cos(primTheta));

	if ( EdepQ[1]>0. && EdepQ[0]==0. && MWPCEnergy[1]>0. && MWPCEnergy[0]==0. )
	  hisW_polE->Fill(primKE,1+0.12*sqrt(primKE*primKE+2.*primKE*510.998928)/(primKE+510.998928)*TMath::Cos(primTheta));

      }   
      delete f;
      
    }
  }


  //Now West polarization
  for (int i=5; i<10; ++i) {
    std::cout << "In Folder " << i << std::endl;
    for (int j=0; j<1000; ++j) {
      f = new TFile(TString::Format("%s/%s%i/analyzed_%i.root",
					   basepath.Data(),
					   (field==TString("good")?"G":
					    (field==TString("bad")?"B":"F")),i,j),"READ");
 
      TTree *t = (TTree*)f->Get("anaTree");
 
      t->SetBranchAddress("EdepQ",EdepQ);
      t->SetBranchAddress("MWPCEnergy",MWPCEnergy);
      t->SetBranchAddress("primKE",&primKE);
      t->SetBranchAddress("primTheta",&primTheta);
      
  
      for (UInt_t i=0; i<t->GetEntriesFast() ; ++i) {

	t->GetEntry(i);

	if ( EdepQ[0]>0. && EdepQ[1]==0. && MWPCEnergy[0]>0. && MWPCEnergy[1]==0. )
	  hisE_polW->Fill(primKE,1-0.12*sqrt(primKE*primKE+2.*primKE*510.998928)/(primKE+510.998928)*TMath::Cos(primTheta));

	if ( EdepQ[1]>0. && EdepQ[0]==0. && MWPCEnergy[1]>0. && MWPCEnergy[0]==0. )
	  hisW_polW->Fill(primKE,1-0.12*sqrt(primKE*primKE+2.*primKE*510.998928)/(primKE+510.998928)*TMath::Cos(primTheta));

      }   
      delete f;
      
    }
  }

  for (UInt_t b=0; b<100; ++b) {
    rE_polE[b] = hisE_polE->GetBinContent(b+1); rE_polE_err[b] = TMath::Sqrt(rE_polE[b]);
    rW_polE[b] = hisW_polE->GetBinContent(b+1); rW_polE_err[b] = TMath::Sqrt(rW_polE[b]);
    rE_polW[b] = hisE_polW->GetBinContent(b+1); rE_polW_err[b] = TMath::Sqrt(rE_polW[b]);
    rW_polW[b] = hisW_polW->GetBinContent(b+1); rW_polW_err[b] = TMath::Sqrt(rW_polW[b]);
  }

  std::vector <Double_t> asymm(100,0.);
  std::vector <Double_t> asymmErr(100,0.);
  
  for (UInt_t b=0;b<100;++b) {
    Double_t R = (rE_polE[b]*rW_polW[b]) / (rE_polW[b]*rW_polE[b]);
    Double_t deltaR = TMath::Sqrt( TMath::Power( rE_polE_err[b]*rW_polW[b]/(rE_polW[b]*rW_polE[b]), 2 ) + 
				   TMath::Power( rW_polW_err[b]*rE_polE[b]/(rE_polW[b]*rW_polE[b]), 2 ) +
				   TMath::Power( rE_polW_err[b]*(rE_polE[b]*rW_polW[b]) / (rE_polW[b]*rE_polW[b]*rW_polE[b]), 2 ) +
				   TMath::Power( rW_polE_err[b]*(rE_polE[b]*rW_polW[b]) / (rE_polW[b]*rW_polE[b]*rW_polE[b]), 2 ) );
    
    
    asymm[b] = (1.-TMath::Sqrt(R))/(1.+TMath::Sqrt(R));
    asymmErr[b] = (deltaR)/(TMath::Sqrt(R)*TMath::Power((TMath::Sqrt(R)+1.),2));
    
  }
  std::cout << " A = " << asymm[20] << " +/- " << asymmErr[20] << std::endl;
  
  ofstream ofile(TString::Format("asymm_%sField.txt",field.Data()).Data());
  ofile << std::setprecision(10);

  for (int b=0;b<100;++b) {
    ofile << 10.*b+5. << "\t" << asymm[b] << "\t" << asymmErr[b] << std::endl;
  }
  ofile.close();
    


};

int main(int argc, char *argv[]) {

  if (argc!=2) {
    std::cout << "usage: ./FieldDipSystematic.exe [field=good,flat] [year]\n";
    exit(0);
  }

  FieldAsymmetry(TString(argv[1]),TString(argv[2]));

}

