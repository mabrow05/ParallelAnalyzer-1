
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

#include "calibrationTools.hh"

// Function to return beta when given the kinetic energy of an electron
Double_t returnBeta(Double_t En) { 
  Double_t me = 510.998928; //rest mass energy of electron in keV
  return sqrt(En*En+2.*En*me)/(En+me);
};


void FieldAsymmetry(TString field,TString year, int startFileNum, int endFileNum) {

  int badFiles[21] = {926,944,1605,2327,2614,2703,3109,3449,4494,4968,
		      5531,5869,6216,6687,6863,7766,9225,9718,9798,9821,9959};
  

  //Load the simulated relationship between EQ and Etrue
  int fakeRunNumber = ( year==TString("2011-2012")?17126:
			( year==TString("2012-2013")?22024:21537 ));
  EreconParameterization eRecon(fakeRunNumber);

  //Initializing the separator
  BackscatterSeparator sep;
  sep.LoadCutCurve(fakeRunNumber);

  TString field_prefix = ( field==TString("flat")?TString("FLAT_FIELD"):
			   (field==TString("dip")?TString("FIELD_DIP"):
			    (field==TString("symm")?TString("SYMMETRIC_FIELD_DIP"):"") ) );

  if ( field_prefix==TString("") ) {
    std::cout << "bad field type\n"; exit(0);
  }

  std::vector <Double_t> rE_polE (100,0.);
  std::vector <Double_t> rE_polE_err (100,0.);  
  std::vector <Double_t> rW_polE (100,0.);
  std::vector <Double_t> rW_polE_err (100,0.);  
  std::vector <Double_t> rE_polW (100,0.);
  std::vector <Double_t> rE_polW_err (100,0.);  
  std::vector <Double_t> rW_polW (100,0.);
  std::vector <Double_t> rW_polW_err (100,0.);  

  TString basepath = TString::Format("/scratch/mabr239/output/%s_1e-3_%s/",field_prefix.Data(),year.Data());

  TFile *f;
  
  TH1D *hisE_polE = new TH1D("hisE_polE","hisE",100,0.,1000.);
  TH1D *hisW_polE = new TH1D("hisW_polE","hisW",100,0.,1000.);

  TH1D *hisE_polW = new TH1D("hisE_polW","hisE",100,0.,1000.);
  TH1D *hisW_polW = new TH1D("hisW_polW","hisW",100,0.,1000.);


  Double_t EdepQ[2],MWPCEnergy[2],time[2],primKE,primTheta;

  //Start with East polarization
  for (int j=startFileNum; j<endFileNum; j+=2) {
    
    bool skip = false;
    for (int n=0;n<21;++n) {
      if (j==badFiles[n]) skip=true;
    }

    if (skip) continue;
    std::cout << "in file " << j << std::endl;

    f = new TFile(TString::Format("%s/analyzed_%i.root",basepath.Data(),j),"READ");
    TTree *t = (TTree*)f->Get("anaTree");
    t->SetBranchAddress("EdepQ",EdepQ);
    t->SetBranchAddress("MWPCEnergy",MWPCEnergy);
    t->SetBranchAddress("primKE",&primKE);
    t->SetBranchAddress("primTheta",&primTheta);
    t->SetBranchAddress("time",time);
    
    for (UInt_t i=0; i<t->GetEntriesFast() ; ++i) {
      
      t->GetEntry(i);
      double erecon = 0.;
      double weight = 1.+0.12*sqrt(primKE*primKE+2.*primKE*510.998928)/(primKE+510.998928)*TMath::Cos(primTheta);
      
      // Type 0
      if ( EdepQ[0]>0. && EdepQ[1]==0. && MWPCEnergy[0]>0. && MWPCEnergy[1]==0. ) {
	erecon = eRecon.getErecon(0,0,EdepQ[0]);
	hisE_polE->Fill(erecon,weight);
      }
      else if ( EdepQ[1]>0. && EdepQ[0]==0. && MWPCEnergy[1]>0. && MWPCEnergy[0]==0. ) {
	erecon = eRecon.getErecon(1,0,EdepQ[1]);
	hisW_polE->Fill(erecon,weight);
      }

      //Type 1
      else if ( EdepQ[0]>0. && EdepQ[1]>0. && MWPCEnergy[0]>0. && MWPCEnergy[1]>0. ) {
	if (time[0]<time[1]) {
	  erecon = eRecon.getErecon(0,1,EdepQ[0]+EdepQ[1]);
	  hisE_polE->Fill(erecon,weight);
	}
	else {
	  erecon = eRecon.getErecon(1,1,EdepQ[0]+EdepQ[1]);
	  hisW_polE->Fill(erecon,weight);
	}
      }
      
      //Type 2/3
      else if ( MWPCEnergy[0]>0. && MWPCEnergy[1]>0. ) {
	int side = -1;
	if (EdepQ[0]>0. && EdepQ[1]==0.) {
	  int type = sep.separate23(MWPCEnergy[0]);
	  erecon = eRecon.getErecon(0,2,EdepQ[0]);
	  side = type==2 ? 1 : 0;
	}
	else if (EdepQ[0]==0. && EdepQ[1]>0.){
	  int type = sep.separate23(MWPCEnergy[1]);
	  erecon = eRecon.getErecon(1,2,EdepQ[1]);
	  side = type==2 ? 0 : 1;
	}
	if ( side==0 ) hisE_polE->Fill(erecon,weight);
	else if ( side==1 ) hisW_polE->Fill(erecon,weight);
      }
      
    }   
    delete f;
    
  }



  //Now West polarization
  
  for (int j=startFileNum+1; j<endFileNum; j+=2) {
    std::cout << "in file " << j << std::endl;

    bool skip = false;
    for (int n=0;n<21;++n) {
      if (j==badFiles[n]) skip=true;
    }

    if (skip) continue;
    
    f = new TFile(TString::Format("%s/analyzed_%i.root",basepath.Data(),j),"READ");
    TTree *t = (TTree*)f->Get("anaTree");
    
    t->SetBranchAddress("EdepQ",EdepQ);
    t->SetBranchAddress("MWPCEnergy",MWPCEnergy);
    t->SetBranchAddress("primKE",&primKE);
    t->SetBranchAddress("primTheta",&primTheta);
    t->SetBranchAddress("time",time);
    
    for (UInt_t i=0; i<t->GetEntriesFast() ; ++i) {
      
      t->GetEntry(i);
      double erecon = 0.;
      double weight = 1-0.12*sqrt(primKE*primKE+2.*primKE*510.998928)/(primKE+510.998928)*TMath::Cos(primTheta);

      // Type 0
      if ( EdepQ[0]>0. && EdepQ[1]==0. && MWPCEnergy[0]>0. && MWPCEnergy[1]==0. ) {
	erecon = eRecon.getErecon(0,0,EdepQ[0]);
	hisE_polW->Fill(erecon,weight);
      }
      else if ( EdepQ[1]>0. && EdepQ[0]==0. && MWPCEnergy[1]>0. && MWPCEnergy[0]==0. ) {
	erecon = eRecon.getErecon(1,0,EdepQ[1]);
	hisW_polW->Fill(erecon,weight);
      }

      //Type 1
      else if ( EdepQ[0]>0. && EdepQ[1]>0. && MWPCEnergy[0]>0. && MWPCEnergy[1]>0. ) {
	if (time[0]<time[1]) {
	  erecon = eRecon.getErecon(0,1,EdepQ[0]+EdepQ[1]);
	  hisE_polW->Fill(erecon,weight);
	}
	else {
	  erecon = eRecon.getErecon(1,1,EdepQ[0]+EdepQ[1]);
	  hisW_polW->Fill(erecon,weight);
	}
      }
      
      //Type 2/3
      else if ( MWPCEnergy[0]>0. && MWPCEnergy[1]>0. ) {
	int side = -1;
	if (EdepQ[0]>0. && EdepQ[1]==0.) {
	  int type = sep.separate23(MWPCEnergy[0]);
	  erecon = eRecon.getErecon(0,2,EdepQ[0]);
	  side = type==2 ? 1 : 0;
	}
	else if (EdepQ[0]==0. && EdepQ[1]>0.){
	  int type = sep.separate23(MWPCEnergy[1]);
	  erecon = eRecon.getErecon(1,2,EdepQ[1]);
	  side = type==2 ? 0 : 1;
	}
	if ( side==0 ) hisE_polW->Fill(erecon,weight);
	else if ( side==1 ) hisW_polW->Fill(erecon,weight);
      } 
    }   
    delete f;
    
  }

  
  for (UInt_t b=0; b<100; ++b) {
    rE_polE[b] = hisE_polE->GetBinContent(b+1); rE_polE_err[b] = TMath::Sqrt(rE_polE[b]);
    rW_polE[b] = hisW_polE->GetBinContent(b+1); rW_polE_err[b] = TMath::Sqrt(rW_polE[b]);
    rE_polW[b] = hisE_polW->GetBinContent(b+1); rE_polW_err[b] = TMath::Sqrt(rE_polW[b]);
    rW_polW[b] = hisW_polW->GetBinContent(b+1); rW_polW_err[b] = TMath::Sqrt(rW_polW[b]);
  }
  
  delete hisE_polE;
  delete hisW_polE;
  delete hisE_polW;
  delete hisW_polW;
  
  std::vector <Double_t> asymm(100,0.);
  std::vector <Double_t> asymmErr(100,0.);
  
  for (UInt_t b=0;b<100;++b) {
    Double_t R = (rE_polW[b]*rW_polE[b]) / (rE_polE[b]*rW_polW[b]);
    Double_t deltaR = TMath::Sqrt( TMath::Power( rE_polW_err[b]*rW_polE[b]/(rE_polE[b]*rW_polW[b]), 2 ) + 
				   TMath::Power( rW_polE_err[b]*rE_polW[b]/(rE_polE[b]*rW_polW[b]), 2 ) +
				   TMath::Power( rE_polE_err[b]*(rE_polW[b]*rW_polE[b]) / (rE_polE[b]*rE_polE[b]*rW_polW[b]), 2 ) +
				   TMath::Power( rW_polW_err[b]*(rE_polW[b]*rW_polE[b]) / (rE_polE[b]*rW_polW[b]*rW_polW[b]), 2 ) );
    
    
    asymm[b] = (1.-TMath::Sqrt(R))/(1.+TMath::Sqrt(R));
    asymmErr[b] = (deltaR)/(TMath::Sqrt(R)*TMath::Power((TMath::Sqrt(R)+1.),2));
    
  }
  std::cout << " A = " << asymm[20] << " +/- " << asymmErr[20] << std::endl;
  
  ofstream ofile(TString::Format("%s_asymm_%sField_files_%i-%i.txt",year.Data(),field.Data(),startFileNum,endFileNum).Data());
  ofile << std::setprecision(10);

  for (int b=0;b<100;++b) {
    ofile << 10.*b+5. << "\t" << asymm[b] << "\t" << asymmErr[b] << std::endl;
  }
  ofile.close();
    


};

int main(int argc, char *argv[]) {

  if (argc!=5) {
    std::cout << "usage: ./FieldDipSystematic.exe [field=dip,flat,symm] [year] [minFileNum] [maxFileNum]\n";
    exit(0);
  }
  
  Int_t startFileNum = atoi(argv[3]);
  Int_t endFileNum = atoi(argv[4]);
  FieldAsymmetry(TString(argv[1]),TString(argv[2]),startFileNum,endFileNum);

}

