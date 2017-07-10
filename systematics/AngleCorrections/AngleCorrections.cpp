
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


void AngleCorr(TString year, int startFileNum, int endFileNum) {

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

  std::vector <Double_t> aveThetaPure(100,0.);
  std::vector <Double_t> PureEntries(100,0.);
  std::vector <Double_t> aveThetaTrigg(100,0.);
  std::vector <Double_t> TriggEntries(100,0.);
  
  TString basepath = TString::Format("%s/",year==TString("2011-2012")?getenv("SIM_2011_2012"):
				     (year==TString("2012-2013")?getenv("SIM_2012_2013"):getenv("SIM_2012_2013_ISOBUTANE")));

  TFile *f;

  Double_t EdepQ[2],MWPCEnergy[2],time[2],primKE,primTheta;

  //Start with East polarization
  for (int j=startFileNum; j<endFileNum; ++j) {
    
    f = new TFile(TString::Format("%s/Beta_polE/analyzed_%i.root",basepath.Data(),j),"READ");
    TTree *t = (TTree*)f->Get("anaTree");
    t->SetBranchAddress("EdepQ",EdepQ);
    t->SetBranchAddress("MWPCEnergy",MWPCEnergy);
    t->SetBranchAddress("primKE",&primKE);
    t->SetBranchAddress("primTheta",&primTheta);
    t->SetBranchAddress("time",time);
    
    for (UInt_t i=0; i<t->GetEntriesFast() ; ++i) {
      
      t->GetEntry(i);

      aveThetaPure[(int)(primKE/10.)]+=TMath::Abs(TMath::Cos(primTheta));
      PureEntries[(int)(primKE/10.)]+=1.;

      double erecon = 0.;
            
      // Type 0
      if ( EdepQ[0]>0. && EdepQ[1]==0. && MWPCEnergy[0]>0. && MWPCEnergy[1]==0. ) {
	erecon = eRecon.getErecon(0,0,EdepQ[0]);
	aveThetaTrigg[(int)(erecon/10.)]+=TMath::Abs(TMath::Cos(primTheta));
	TriggEntries[(int)(erecon/10.)]+=1.;

      }
      else if ( EdepQ[1]>0. && EdepQ[0]==0. && MWPCEnergy[1]>0. && MWPCEnergy[0]==0. ) {
	erecon = eRecon.getErecon(1,0,EdepQ[1]);
	aveThetaTrigg[(int)(erecon/10.)]+=TMath::Abs(TMath::Cos(primTheta));
	TriggEntries[(int)(erecon/10.)]+=1.;
      }

      //Type 1
      else if ( EdepQ[0]>0. && EdepQ[1]>0. && MWPCEnergy[0]>0. && MWPCEnergy[1]>0. ) {
	if (time[0]<time[1]) {
	  erecon = eRecon.getErecon(0,1,EdepQ[0]+EdepQ[1]);
	  aveThetaTrigg[(int)(erecon/10.)]+=TMath::Abs(TMath::Cos(primTheta));
	  TriggEntries[(int)(erecon/10.)]+=1.;
	}
	else {
	  erecon = eRecon.getErecon(1,1,EdepQ[0]+EdepQ[1]);
	  aveThetaTrigg[(int)(erecon/10.)]+=TMath::Abs(TMath::Cos(primTheta));
	  TriggEntries[(int)(erecon/10.)]+=1.;
	
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
	aveThetaTrigg[(int)(erecon/10.)]+=TMath::Abs(TMath::Cos(primTheta));
	TriggEntries[(int)(erecon/10.)]+=1.;
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
    
    for (UInt_t i=0; i<t->GetEntriesFast() ; ++i) {
      
      t->GetEntry(i);

      aveThetaPure[(int)(primKE/10.)]+=TMath::Abs(TMath::Cos(primTheta));
      PureEntries[(int)(primKE/10.)]+=1.;

      double erecon = 0.;
            
      // Type 0
      if ( EdepQ[0]>0. && EdepQ[1]==0. && MWPCEnergy[0]>0. && MWPCEnergy[1]==0. ) {
	erecon = eRecon.getErecon(0,0,EdepQ[0]);
	aveThetaTrigg[(int)(erecon/10.)]+=TMath::Abs(TMath::Cos(primTheta));
	TriggEntries[(int)(erecon/10.)]+=1.;

      }
      else if ( EdepQ[1]>0. && EdepQ[0]==0. && MWPCEnergy[1]>0. && MWPCEnergy[0]==0. ) {
	erecon = eRecon.getErecon(1,0,EdepQ[1]);
	aveThetaTrigg[(int)(erecon/10.)]+=TMath::Abs(TMath::Cos(primTheta));
	TriggEntries[(int)(erecon/10.)]+=1.;
      }

      //Type 1
      else if ( EdepQ[0]>0. && EdepQ[1]>0. && MWPCEnergy[0]>0. && MWPCEnergy[1]>0. ) {
	if (time[0]<time[1]) {
	  erecon = eRecon.getErecon(0,1,EdepQ[0]+EdepQ[1]);
	  aveThetaTrigg[(int)(erecon/10.)]+=TMath::Abs(TMath::Cos(primTheta));
	  TriggEntries[(int)(erecon/10.)]+=1.;
	}
	else {
	  erecon = eRecon.getErecon(1,1,EdepQ[0]+EdepQ[1]);
	  aveThetaTrigg[(int)(erecon/10.)]+=TMath::Abs(TMath::Cos(primTheta));
	  TriggEntries[(int)(erecon/10.)]+=1.;
	
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
	aveThetaTrigg[(int)(erecon/10.)]+=TMath::Abs(TMath::Cos(primTheta));
	TriggEntries[(int)(erecon/10.)]+=1.;
      }
      
    }   
    delete f;
    
  }



  for (UInt_t b=0; b<100; ++b) {
    aveThetaPure[b]= ( PureEntries[b]>0.?aveThetaPure[b]/PureEntries[b]:0. );
    aveThetaTrigg[b]= ( TriggEntries[b]>0.?aveThetaTrigg[b]/TriggEntries[b]:0. );
  }
  
  
  ofstream ofile(TString::Format("angleCorr_%s.txt",year.Data()).Data());
  ofile << std::setprecision(10);

  for (int b=0;b<100;++b) {
    ofile << 10.*b+5. << "\t" << aveThetaPure[b]/aveThetaTrigg[b] << std::endl;
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

