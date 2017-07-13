#include <TRandom3.h>
#include <TString.h>
#include <vector>
#include <TMath.h>
#include <TFile.h>
#include <TH1D.h>
#include <TTree.h>
#include <TBranch.h>
#include <TLeaf.h>
#include <TF1.h>

#include <iostream>
#include <cstdlib>
#include <fstream>
#include <vector>
#include <iomanip>

#include "calibrationTools.hh"


class trigg {

public: 
  trigg(TString year) ;
  ~trigg();
  bool triggerE(Double_t *eq,Double_t rand) {return (fE->EvalPar(eq,pE)>rand)?true:false;};
  bool triggerW(Double_t *eq,Double_t rand) {return (fW->EvalPar(eq,pW)>rand)?true:false;};
protected:
  Double_t pE[2];
  Double_t pW[2];
  TF1 *fE,*fW;
};

trigg::trigg(TString year) {
  
  if (year==TString("2011-2012")) {   
    pE[0] = -18.3615;
    pE[1] = 47.5979;
    pW[0] = -15.8085;
    pW[1] = 46.3180;
  }
  else {
    pE[0] = -13.6101;
    pE[1] = 45.7509;
    pW[0] = -16.9604;
    pW[1] = 46.5911;
  }
  
  fE = new TF1("fE","0.5+0.5*TMath::Erf((x-[0])/[1])",0.,1000.);
  fW = new TF1("fW","0.5+0.5*TMath::Erf((x-[0])/[1])",0.,1000.);
  
};

trigg::~trigg() {delete fE; delete fW;}


// Function to return beta when given the kinetic energy of an electron
Double_t returnBeta(Double_t En) { 
  Double_t me = 510.998928; //rest mass energy of electron in keV
  return sqrt(En*En+2.*En*me)/(En+me);
};


void EfficiencyCorr(Double_t eastEff, Double_t westEff, TString year, int startFileNum, int endFileNum) {

  bool usePrimValues=false;

  double fidCut = 50.;

  //Load the simulated relationship between EQ and Etrue
  int fakeRunNumber = ( year==TString("2011-2012")?17126:
  			( year==TString("2012-2013")?22024:21537 ));
  EreconParameterization eRecon(fakeRunNumber);

  //Initializing the separator
  BackscatterSeparator sep;
  sep.LoadCutCurve(fakeRunNumber);

  trigg scintTrigg(year);

  //For smearing
  double alpha = 0.4; // nPE/keV of roughly 400 PE per 1 GeV
  Double_t g_d = 4.*16.;
  Double_t g_rest = 4.*12500.;
  
  Double_t tot_rE_polE=0.;
  Double_t tot_rW_polE=0.;
  Double_t tot_rE_polW=0.;
  Double_t tot_rW_polW=0.;

  Double_t tot_rE_polE_err=0.;
  Double_t tot_rW_polE_err=0.;
  Double_t tot_rE_polW_err=0.;
  Double_t tot_rW_polW_err=0.;

  std::vector <Double_t> rE_polE (100,0.);
  std::vector <Double_t> rE_polE_err (100,0.);  
  std::vector <Double_t> rW_polE (100,0.);
  std::vector <Double_t> rW_polE_err (100,0.);  
  std::vector <Double_t> rE_polW (100,0.);
  std::vector <Double_t> rE_polW_err (100,0.);  
  std::vector <Double_t> rW_polW (100,0.);
  std::vector <Double_t> rW_polW_err (100,0.);  

  std::vector <UInt_t> eastType0(90,0);
  std::vector <UInt_t> eastType1(90,0);
  std::vector <UInt_t> eastType23(90,0);
  std::vector <UInt_t> eastType3(90,0);

  std::vector <UInt_t> westType0(90,0);
  std::vector <UInt_t> westType1(90,0);
  std::vector <UInt_t> westType23(90,0);
  std::vector <UInt_t> westType3(90,0);

  TString basepath = TString::Format("%s/",(year==TString("2011-2012")?getenv("SIM_2011_2012"):(year==TString("2012-2013")?getenv("SIM_2012_2013"):getenv("SIM_2012_2013_ISOBUTANE"))));

  TFile *f;
  
  TH1D *hisE_polE = new TH1D("hisE_polE","hisE",100,0.,1000.);
  TH1D *hisW_polE = new TH1D("hisW_polE","hisW",100,0.,1000.);

  TH1D *hisE_polW = new TH1D("hisE_polW","hisE",100,0.,1000.);
  TH1D *hisW_polW = new TH1D("hisW_polW","hisW",100,0.,1000.);

  TRandom3 rand;

  Double_t EdepQ[2],MWPCEnergy[2],time[2],primKE,primTheta,ScintPosE[3],ScintPosW[3],keInSD[24],primPos[4];

  //Start with East polarization
  for (int j=startFileNum; j<endFileNum; ++j) {
    std::cout << "East Pol file " << j << std::endl;

    f = new TFile(TString::Format("%s/Beta_polE/analyzed_%i.root",basepath.Data(),j),"READ");
    TTree *t = (TTree*)f->Get("anaTree");
    t->SetBranchAddress("EdepQ",EdepQ);
    t->SetBranchAddress("MWPCEnergy",MWPCEnergy);
    t->SetBranchAddress("primKE",&primKE);
    t->SetBranchAddress("primTheta",&primTheta);
    t->SetBranchAddress("time",time);
    t->GetBranch("ScintPos")->GetLeaf("ScintPosE")->SetAddress(ScintPosE);
    t->GetBranch("ScintPos")->GetLeaf("ScintPosW")->SetAddress(ScintPosW);
    t->SetBranchAddress("keInSD",keInSD);
    t->SetBranchAddress("primPos",primPos);
    
    for (UInt_t i=0; i<t->GetEntriesFast() ; ++i) {
      
      t->GetEntry(i);
      //double erecon = 0.;
      double rE = TMath::Sqrt(ScintPosE[0]*ScintPosE[0]+ScintPosE[1]*ScintPosE[1])*TMath::Sqrt(0.6)*10.;
      double rW = TMath::Sqrt(ScintPosW[0]*ScintPosW[0]+ScintPosW[1]*ScintPosW[1])*TMath::Sqrt(0.6)*10.;

      double rPrim = primPos[3]*100.;

      if (rE>fidCut || rW>fidCut || rPrim>fidCut) continue;

      
      double cosTheta = TMath::Cos(primTheta);
      int realside = cosTheta>0.?1:0;
      int windowSide = (keInSD[5]>keInSD[15])?0:1;
      int side = -1;
      int type=-1;
      bool trigger = false;

      EdepQ[0] = (1./(alpha*g_d*g_rest)) * (rand.Poisson(g_rest*rand.Poisson(g_d*rand.Poisson(alpha*EdepQ[0]))));
      EdepQ[1] = (1./(alpha*g_d*g_rest)) * (rand.Poisson(g_rest*rand.Poisson(g_d*rand.Poisson(alpha*EdepQ[1]))));
	  
      
      bool mwpcEastTrigg = (rand.Rndm()<eastEff && MWPCEnergy[0]>0.)?true:false;
      bool mwpcWestTrigg = (rand.Rndm()<westEff && MWPCEnergy[1]>0.)?true:false;
      bool scintEastTrigg = EdepQ[0]>0. && scintTrigg.triggerE(&EdepQ[0],rand.Rndm());
      bool scintWestTrigg = EdepQ[1]>0. && scintTrigg.triggerW(&EdepQ[1],rand.Rndm());

      // Type 0
      if ( scintEastTrigg && !scintWestTrigg && mwpcEastTrigg && !mwpcWestTrigg ) {
	trigger=true;
	side=0;
	type=0;
      }
      else if ( !scintEastTrigg && scintWestTrigg && !mwpcEastTrigg && mwpcWestTrigg ) {
	trigger=true;
	side=1;
	type=0;
      }
      //Type 1
      else if ( scintEastTrigg && scintWestTrigg && mwpcEastTrigg && mwpcWestTrigg ) {
	trigger=true;
	side=time[0]<time[1]?0:1;
	type=1;
      }
      //Type 2/3
      else if ( mwpcEastTrigg && mwpcWestTrigg && ( scintEastTrigg || scintWestTrigg ) ) {
	trigger=true;
	side=scintEastTrigg?0:1;
	type=2;
      }
      
      if (trigger) {
	double evis = type==1?(EdepQ[0]+EdepQ[1]):EdepQ[side];
	double erecon = eRecon.getErecon(side,type,evis);
	if (type==2) {
	  if (side==0) {
	      type = sep.separate23(MWPCEnergy[0]);
	      side = type==2 ? 1 : 0;
	    }
	    else if (side==1) {
	      type = sep.separate23(MWPCEnergy[1]);
	      side = type==2 ? 0 : 1;
	    }
	}
	if (!usePrimValues) {
	  if (side==0) { hisE_polE->Fill(erecon); tot_rE_polE+=(erecon>220. && erecon<670.)?1.:0.;}
	  else if (side==1) { hisW_polE->Fill(erecon); tot_rW_polE+=(erecon>220. && erecon<670.)?1.:0.;}
	}
	else {
	  if (realside==0) { hisE_polE->Fill(primKE); tot_rE_polE+=(primKE>220. && primKE<670.)?1.:0.;}
	  else if (realside==1) { hisW_polE->Fill(primKE); tot_rW_polE+=(primKE>220. && primKE<670.)?1.:0.;}
	}
      }
    }  
    delete f;
    
  }



  //Now West polarization
  
  for (int j=startFileNum; j<endFileNum; ++j) {
    std::cout << "West Pol file " << j << std::endl;
    
    f = new TFile(TString::Format("%s/Beta_polW/analyzed_%i.root",basepath.Data(),j),"READ");
    TTree *t = (TTree*)f->Get("anaTree");
    
    t->SetBranchAddress("EdepQ",EdepQ);
    t->SetBranchAddress("MWPCEnergy",MWPCEnergy);
    t->SetBranchAddress("primKE",&primKE);
    t->SetBranchAddress("primTheta",&primTheta);
    t->SetBranchAddress("time",time);
    t->GetBranch("ScintPos")->GetLeaf("ScintPosE")->SetAddress(ScintPosE);
    t->GetBranch("ScintPos")->GetLeaf("ScintPosW")->SetAddress(ScintPosW);
    t->SetBranchAddress("keInSD",keInSD);
    t->SetBranchAddress("primPos",primPos);
 
    for (UInt_t i=0; i<t->GetEntriesFast() ; ++i) {
      
      t->GetEntry(i);
      //double erecon = 0.;
      double rE = TMath::Sqrt(ScintPosE[0]*ScintPosE[0]+ScintPosE[1]*ScintPosE[1])*TMath::Sqrt(0.6)*10.;
      double rW = TMath::Sqrt(ScintPosW[0]*ScintPosW[0]+ScintPosW[1]*ScintPosW[1])*TMath::Sqrt(0.6)*10.;

      double rPrim = primPos[3]*100.;

      if (rE>fidCut || rW>fidCut || rPrim>fidCut) continue;

      double weight = 1-0.1184*sqrt(primKE*primKE+2.*primKE*510.998928)/(primKE+510.998928)*TMath::Cos(primTheta);
      double cosTheta = TMath::Cos(primTheta);
      int realside = cosTheta>0.?1:0;
      int windowSide = (keInSD[5]>keInSD[15])?0:1;
      int side = -1;
      int type=-1;
      bool trigger = false;
     
      EdepQ[0] = (1./(alpha*g_d*g_rest)) * (rand.Poisson(g_rest*rand.Poisson(g_d*rand.Poisson(alpha*EdepQ[0]))));
      EdepQ[1] = (1./(alpha*g_d*g_rest)) * (rand.Poisson(g_rest*rand.Poisson(g_d*rand.Poisson(alpha*EdepQ[1]))));

      
      bool mwpcEastTrigg = (rand.Rndm()<eastEff && MWPCEnergy[0]>0.)?true:false;
      bool mwpcWestTrigg = (rand.Rndm()<westEff && MWPCEnergy[1]>0.)?true:false;
      bool scintEastTrigg = EdepQ[0]>0. && scintTrigg.triggerE(&EdepQ[0],rand.Rndm());
      bool scintWestTrigg = EdepQ[1]>0. && scintTrigg.triggerW(&EdepQ[1],rand.Rndm());

      // Type 0
      if ( scintEastTrigg && !scintWestTrigg && mwpcEastTrigg && !mwpcWestTrigg ) {
	trigger=true;
	side=0;
	type=0;
      }
      else if ( !scintEastTrigg && scintWestTrigg && !mwpcEastTrigg && mwpcWestTrigg ) {
	trigger=true;
	side=1;
	type=0;
      }
      //Type 1
      else if ( scintEastTrigg && scintWestTrigg && mwpcEastTrigg && mwpcWestTrigg ) {
	trigger=true;
	side=time[0]<time[1]?0:1;
	type=1;
      }
      //Type 2/3
      else if ( mwpcEastTrigg && mwpcWestTrigg && ( scintEastTrigg || scintWestTrigg ) ) {
	trigger=true;
	side=scintEastTrigg?0:1;
	type=2;
      }
      
      if (trigger) {
	double evis = type==1?(EdepQ[0]+EdepQ[1]):EdepQ[side];
	double erecon = eRecon.getErecon(side,type,evis);
	if (type==2) {
	  if (side==0) {
	      type = sep.separate23(MWPCEnergy[0]);
	      side = type==2 ? 1 : 0;
	    }
	    else if (side==1) {
	      type = sep.separate23(MWPCEnergy[1]);
	      side = type==2 ? 0 : 1;
	    }
	}
	if (!usePrimValues) {
	  if (side==0) { hisE_polW->Fill(erecon); tot_rE_polW+=(erecon>220. && erecon<670.)?1.:0.;}
	  else if (side==1) { hisW_polW->Fill(erecon); tot_rW_polW+=(erecon>220. && erecon<670.)?1.:0.;}
	}
	else {
	  if (realside==0) { hisE_polW->Fill(primKE); tot_rE_polW+=(primKE>220. && primKE<670.)?1.:0.;}
	  else if (realside==1) { hisW_polW->Fill(primKE); tot_rW_polW+=(primKE>220. && primKE<670.)?1.:0.;}
	}
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

  tot_rE_polE_err = TMath::Sqrt(tot_rE_polE);
  tot_rW_polE_err = TMath::Sqrt(tot_rW_polE);
  tot_rE_polW_err = TMath::Sqrt(tot_rE_polW);
  tot_rW_polW_err = TMath::Sqrt(tot_rW_polW);
  
  delete hisE_polE;
  delete hisW_polE;
  delete hisE_polW;
  delete hisW_polW;
  
  std::vector <Double_t> asymm(100,0.);
  std::vector <Double_t> asymmErr(100,0.);

  Double_t tot_asymm=0., tot_asymmErr=0.;
  
  for (UInt_t b=0;b<100;++b) {
    Double_t R = (rE_polW[b]*rW_polE[b]) / (rE_polE[b]*rW_polW[b]);
    Double_t deltaR = TMath::Sqrt( TMath::Power( rE_polW_err[b]*rW_polE[b]/(rE_polE[b]*rW_polW[b]), 2 ) + 
				   TMath::Power( rW_polE_err[b]*rE_polW[b]/(rE_polE[b]*rW_polW[b]), 2 ) +
				   TMath::Power( rE_polE_err[b]*(rE_polW[b]*rW_polE[b]) / (rE_polE[b]*rE_polE[b]*rW_polW[b]), 2 ) +
				   TMath::Power( rW_polW_err[b]*(rE_polW[b]*rW_polE[b]) / (rE_polE[b]*rW_polW[b]*rW_polW[b]), 2 ) );
    
    
    asymm[b] = (1.-TMath::Sqrt(R))/(1.+TMath::Sqrt(R));
    asymmErr[b] = (deltaR)/(TMath::Sqrt(R)*TMath::Power((TMath::Sqrt(R)+1.),2));
    
  }

  Double_t R = (tot_rE_polW*tot_rW_polE) / (tot_rE_polE*tot_rW_polW);
  Double_t deltaR = TMath::Sqrt( TMath::Power( tot_rE_polW_err*tot_rW_polE/(tot_rE_polE*tot_rW_polW), 2 ) + 
				 TMath::Power( tot_rW_polE_err*tot_rE_polW/(tot_rE_polE*tot_rW_polW), 2 ) +
				 TMath::Power( tot_rE_polE_err*(tot_rE_polW*tot_rW_polE) / (tot_rE_polE*tot_rE_polE*tot_rW_polW), 2 ) +
				 TMath::Power( tot_rW_polW_err*(tot_rE_polW*tot_rW_polE) / (tot_rE_polE*tot_rW_polW*tot_rW_polW), 2 ) );
  
  
  tot_asymm = (1.-TMath::Sqrt(R))/(1.+TMath::Sqrt(R));
  tot_asymmErr = (deltaR)/(TMath::Sqrt(R)*TMath::Power((TMath::Sqrt(R)+1.),2));
  std::cout << " A = " << tot_asymm << " +/- " << tot_asymmErr << std::endl;
  
  ofstream ofile(TString::Format("%s_asymm%s_E%0.3f_W%0.3f.txt",year.Data(),
				 usePrimValues?"PrimKE":"Erecon",eastEff,westEff).Data());
  ofile << std::setprecision(10);

  for (int b=0;b<100;++b) {
    ofile << 10.*b+5. << "\t" << asymm[b] << "\t" << asymmErr[b] << std::endl;
  }
  ofile.close();

 


};

int main(int argc, char *argv[]) {

  if (argc!=6) {
    std::cout << "usage: ./EfficiencyCorr.exe [east eff] [west eff] [year] [minFileNum] [maxFileNum]\n";
    exit(0);
  }
  
  Int_t startFileNum = atoi(argv[4]);
  Int_t endFileNum = atoi(argv[5]);
  EfficiencyCorr(atof(argv[1]),atof(argv[2]),TString(argv[3]),startFileNum,endFileNum);

}

