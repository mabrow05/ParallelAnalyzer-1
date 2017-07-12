
#include <TString.h>
#include <vector>
#include <TMath.h>
#include <TFile.h>
#include <TH1D.h>
#include <TTree.h>
#include <TBranch.h>
#include <TLeaf.h>

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


void AngleCorr(TString year, int startFileNum, int endFileNum) {

  //Load the simulated relationship between EQ and Etrue
  int fakeRunNumber = ( year==TString("2011-2012")?17126:
			( year==TString("2012-2013")?22024:21537 ));
  EreconParameterization eRecon(fakeRunNumber);

  //Initializing the separator
  BackscatterSeparator sep;
  sep.LoadCutCurve(fakeRunNumber);

  std::vector <Double_t> aveThetaPureE(100,0.);
  std::vector <Double_t> PureEntriesE(100,0.);
  std::vector <Double_t> aveThetaTriggE(100,0.);
  std::vector <Double_t> TriggEntriesE(100,0.);

  std::vector <Double_t> aveThetaPureW(100,0.);
  std::vector <Double_t> PureEntriesW(100,0.);
  std::vector <Double_t> aveThetaTriggW(100,0.);
  std::vector <Double_t> TriggEntriesW(100,0.);
  
  TString basepath = TString::Format("%s/",year==TString("2011-2012")?getenv("SIM_2011_2012"):
				     (year==TString("2012-2013")?getenv("SIM_2012_2013"):getenv("SIM_2012_2013_ISOBUTANE")));

  TFile *f;

  Double_t EdepQ[2],MWPCEnergy[2],time[2],primKE,primTheta,ScintPosE[3],ScintPosW[3],primPos[4];

  //Start with East polarization
  for (int j=startFileNum; j<endFileNum; ++j) {
    
    f = new TFile(TString::Format("%s/Beta_polE/analyzed_%i.root",basepath.Data(),j),"READ");
    TTree *t = (TTree*)f->Get("anaTree");
    t->SetBranchAddress("EdepQ",EdepQ);
    t->SetBranchAddress("MWPCEnergy",MWPCEnergy);
    t->SetBranchAddress("primKE",&primKE);
    t->SetBranchAddress("primTheta",&primTheta);
    t->SetBranchAddress("time",time);
    t->GetBranch("ScintPos")->GetLeaf("ScintPosE")->SetAddress(ScintPosE);
    t->GetBranch("ScintPos")->GetLeaf("ScintPosW")->SetAddress(ScintPosW);
    t->SetBranchAddress("primPos",primPos);

    std::cout << "In East file " << j << std::endl;
    
    for (UInt_t i=0; i<t->GetEntriesFast() ; ++i) {
      
      t->GetEntry(i);
      
      Double_t rE = TMath::Sqrt(ScintPosE[0]*ScintPosE[0]+ScintPosE[1]*ScintPosE[1])*TMath::Sqrt(0.6)*10.;
      Double_t rW = TMath::Sqrt(ScintPosW[0]*ScintPosW[0]+ScintPosW[1]*ScintPosW[1])*TMath::Sqrt(0.6)*10.;
      Double_t rPrim = primPos[3]*100.;

      if (rE>50. || rW>50. || rPrim>50.) continue;

      Int_t side = TMath::Cos(primTheta)>0.?1:0;

      if (side==0) {
	aveThetaPureE[(int)(primKE/10.)]+=TMath::Abs(TMath::Cos(primTheta));
	PureEntriesE[(int)(primKE/10.)]+=1.;
      }
      else if (side==1) {
	aveThetaPureW[(int)(primKE/10.)]+=TMath::Abs(TMath::Cos(primTheta));
	PureEntriesW[(int)(primKE/10.)]+=1.;
      }

      bool trigger = false;
      // Type 0
      if ( EdepQ[0]>0. && EdepQ[1]==0. && MWPCEnergy[0]>0. && MWPCEnergy[1]==0. ) {
	trigger=true;
      }
      else if ( EdepQ[1]>0. && EdepQ[0]==0. && MWPCEnergy[1]>0. && MWPCEnergy[0]==0. ) {
	trigger=true;
      }

      //Type 1
      else if ( EdepQ[0]>0. && EdepQ[1]>0. && MWPCEnergy[0]>0. && MWPCEnergy[1]>0. ) {
	trigger=true;
      }
      
      //Type 2/3
      else if ( MWPCEnergy[0]>0. && MWPCEnergy[1]>0. && (EdepQ[0]>0. || EdepQ[1]>0.) ) {
	trigger=true;
      }

      if (trigger && side==0) {
	aveThetaTriggE[(int)(primKE/10.)]+=TMath::Abs(TMath::Cos(primTheta));
	TriggEntriesE[(int)(primKE/10.)]+=1.;
      }
      else if (trigger && side==1) {
	aveThetaTriggW[(int)(primKE/10.)]+=TMath::Abs(TMath::Cos(primTheta));
	TriggEntriesW[(int)(primKE/10.)]+=1.;
      }

    }   
    delete f;
    
  }


  
  //Then West polarization
  for (int j=startFileNum; j<endFileNum; ++j) {
    
    f = new TFile(TString::Format("%s/Beta_polW/analyzed_%i.root",basepath.Data(),j),"READ");
    TTree *t = (TTree*)f->Get("anaTree");
    t->SetBranchAddress("EdepQ",EdepQ);
    t->SetBranchAddress("MWPCEnergy",MWPCEnergy);
    t->SetBranchAddress("primKE",&primKE);
    t->SetBranchAddress("primTheta",&primTheta);
    t->SetBranchAddress("time",time);
    t->GetBranch("ScintPos")->GetLeaf("ScintPosE")->SetAddress(ScintPosE);
    t->GetBranch("ScintPos")->GetLeaf("ScintPosW")->SetAddress(ScintPosW);
    t->SetBranchAddress("primPos",primPos);
	
    std::cout << "In West file " << j << std::endl;
    
    for (UInt_t i=0; i<t->GetEntriesFast() ; ++i) {
      
      t->GetEntry(i);
      Double_t rE = TMath::Sqrt(ScintPosE[0]*ScintPosE[0]+ScintPosE[1]*ScintPosE[1])*TMath::Sqrt(0.6)*10.;
      Double_t rW = TMath::Sqrt(ScintPosW[0]*ScintPosW[0]+ScintPosW[1]*ScintPosW[1])*TMath::Sqrt(0.6)*10.;
      Double_t rPrim = primPos[3]*100.;

      if (rE>50. || rW>50. || rPrim>50.) continue;

      Int_t side = TMath::Cos(primTheta)>0.?1:0;

      if (side==0) {
	aveThetaPureE[(int)(primKE/10.)]+=TMath::Abs(TMath::Cos(primTheta));
	PureEntriesE[(int)(primKE/10.)]+=1.;
      }
      else if (side==1) {
	aveThetaPureW[(int)(primKE/10.)]+=TMath::Abs(TMath::Cos(primTheta));
	PureEntriesW[(int)(primKE/10.)]+=1.;
      }

      bool trigger = false;
      // Type 0
      if ( EdepQ[0]>0. && EdepQ[1]==0. && MWPCEnergy[0]>0. && MWPCEnergy[1]==0. ) {
	trigger=true;
      }
      else if ( EdepQ[1]>0. && EdepQ[0]==0. && MWPCEnergy[1]>0. && MWPCEnergy[0]==0. ) {
	trigger=true;
      }

      //Type 1
      else if ( EdepQ[0]>0. && EdepQ[1]>0. && MWPCEnergy[0]>0. && MWPCEnergy[1]>0. ) {
	trigger=true;
      }
      
      //Type 2/3
      else if ( MWPCEnergy[0]>0. && MWPCEnergy[1]>0. && (EdepQ[0]>0. || EdepQ[1]>0.) ) {
	trigger=true;
      }

      if (trigger && side==0) {
	aveThetaTriggE[(int)(primKE/10.)]+=TMath::Abs(TMath::Cos(primTheta));
	TriggEntriesE[(int)(primKE/10.)]+=1.;
      }
      else if (trigger && side==1) {
	aveThetaTriggW[(int)(primKE/10.)]+=TMath::Abs(TMath::Cos(primTheta));
	TriggEntriesW[(int)(primKE/10.)]+=1.;
      }
      
    }   
    delete f;
    
  }



  for (UInt_t b=0; b<100; ++b) {
    aveThetaPureE[b]= ( PureEntriesE[b]>0.?aveThetaPureE[b]/PureEntriesE[b]:0. );
    aveThetaTriggE[b]= ( TriggEntriesE[b]>0.?aveThetaTriggE[b]/TriggEntriesE[b]:0. );
    aveThetaPureW[b]= ( PureEntriesW[b]>0.?aveThetaPureW[b]/PureEntriesW[b]:0. );
    aveThetaTriggW[b]= ( TriggEntriesW[b]>0.?aveThetaTriggW[b]/TriggEntriesW[b]:0. );
  }
  
  
  ofstream ofile(TString::Format("angleCorr_%s.txt",year.Data()).Data());
  ofile << std::setprecision(10);

  for (int b=0;b<100;++b) {
    double deltaEast = aveThetaTriggE[b]/aveThetaPureE[b];
    double deltaWest = aveThetaTriggW[b]/aveThetaPureW[b];

    ofile << 10.*b+5. << "\t" 
	  << deltaEast << "\t" 
	  << deltaWest << "\t"
	  << 1./(7./8.*deltaWest+1./8.*deltaEast)-1. << std::endl;
  }
  ofile.close();
    


};

int main(int argc, char *argv[]) {

  if (argc!=4) {
    std::cout << "usage: ./AngleCorrections.exe [year] [minFileNum] [maxFileNum]\n";
    exit(0);
  }
  
  Int_t startFileNum = atoi(argv[2]);
  Int_t endFileNum = atoi(argv[3]);
  AngleCorr(TString(argv[1]),startFileNum,endFileNum);

}

