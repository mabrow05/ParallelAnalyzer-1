
#include <TString.h>
#include <vector>
#include <TMath.h>
#include <TFile.h>
#include <TH1D.h>
#include <TTree.h>
#include <TBranch.h>
#include <TLeaf.h>
#include <TRandom3.h>

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


void AngleCorr(TString year, int startFileNum, int endFileNum) {

  //Load the simulated relationship between EQ and Etrue
  int fakeRunNumber = ( year==TString("2011-2012")?17126:
			( year==TString("2012-2013")?22024:21537 ));
  EreconParameterization eRecon(fakeRunNumber);

  //Initializing the separator
  BackscatterSeparator sep;
  sep.LoadCutCurve(fakeRunNumber);

  trigg scintTrigg(year);

  //For smearing on a detector basis (roughly)                                                                        
  double alpha = 0.4; // nPE/keV of roughly 400 PE per 1 GeV                                                                                  
  Double_t g_d = 4.*16.;
  Double_t g_rest = 4.*12500.;

  std::vector <Double_t> aveThetaPure(100,0.);
  std::vector <Double_t> PureEntries(100,0.);
  std::vector <Double_t> aveThetaTrigg(100,0.);
  std::vector <Double_t> TriggEntries(100,0.);

  std::vector <Double_t> aveThetaPureE(100,0.);
  std::vector <Double_t> PureEntriesE(100,0.);
  std::vector <Double_t> aveThetaTriggE(100,0.);
  std::vector <Double_t> TriggEntriesE(100,0.);

  std::vector <Double_t> aveThetaPureW(100,0.);
  std::vector <Double_t> PureEntriesW(100,0.);
  std::vector <Double_t> aveThetaTriggW(100,0.);
  std::vector <Double_t> TriggEntriesW(100,0.);

  //These are the successive corrected rates
  std::vector <std::vector <std::vector <Double_t> > > rate_polE (4, std::vector<std::vector<Double_t>>(2, std::vector<Double_t>(100,0.)));
  std::vector <std::vector <std::vector <Double_t> > > rate_polE_err (4, std::vector<std::vector<Double_t>>(2, std::vector<Double_t>(100,0.)));
  std::vector <std::vector <std::vector <Double_t> > > rate_polW (4, std::vector<std::vector<Double_t>>(2, std::vector<Double_t>(100,0.)));
  std::vector <std::vector <std::vector <Double_t> > > rate_polW_err (4, std::vector<std::vector<Double_t>>(2, std::vector<Double_t>(100,0.)));

  std::vector <std::vector<Double_t>> UNCORR_rate_polE (2, std::vector<Double_t>(100,0.));  
  std::vector <std::vector<Double_t>> UNCORR_rate_polE_err (2, std::vector<Double_t>(100,0.));
  std::vector <std::vector<Double_t>> UNCORR_rate_polW (2, std::vector<Double_t>(100,0.));  
  std::vector <std::vector<Double_t>> UNCORR_rate_polW_err (2, std::vector<Double_t>(100,0.)); 
  
  TString basepath = TString::Format("%s/",year==TString("2011-2012")?getenv("SIM_2011_2012"):
				     (year==TString("2012-2013")?getenv("SIM_2012_2013"):getenv("SIM_2012_2013_ISOBUTANE")));

  TFile *f;

  Double_t EdepQ[2],MWPCEnergy[2],time[2],primKE,primTheta,ScintPosE[3],ScintPosW[3],primPos[4];

  TRandom3 rand;

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

      Int_t realSide = TMath::Cos(primTheta)>0.?1:0;
      Int_t side=-1;
      Int_t type=-1;
      if (realSide==0) {
	aveThetaPure[(int)(primKE/10.)]+=TMath::Abs(returnBeta(primKE)*TMath::Cos(primTheta));
	PureEntries[(int)(primKE/10.)]+=1.;
      }
      else if (realSide==1) {
	aveThetaPure[(int)(primKE/10.)]+=TMath::Abs(returnBeta(primKE)*TMath::Cos(primTheta));
	PureEntries[(int)(primKE/10.)]+=1.;
      }

      bool trigger = false;

      EdepQ[0] = (1./(alpha*g_d*g_rest)) * (rand.Poisson(g_rest*rand.Poisson(g_d*rand.Poisson(alpha*EdepQ[0]))));
      EdepQ[1] = (1./(alpha*g_d*g_rest)) * (rand.Poisson(g_rest*rand.Poisson(g_d*rand.Poisson(alpha*EdepQ[1]))));


      bool mwpcEastTrigg = (MWPCEnergy[0]>0.)?true:false;
      bool mwpcWestTrigg = (MWPCEnergy[1]>0.)?true:false;
      bool scintEastTrigg = EdepQ[0]>0. && scintTrigg.triggerE(&EdepQ[0],rand.Rndm());//EdepQ[0]>0.;// 
      bool scintWestTrigg = EdepQ[1]>0. && scintTrigg.triggerW(&EdepQ[1],rand.Rndm());//

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

	if (type==0) {
	  UNCORR_rate_polE[side][erecon/10.]+=1.;
	  rate_polE[0][realSide][erecon/10.]+=1.;
	  rate_polE[1][realSide][erecon/10.]+=1.;
	  rate_polE[2][realSide][erecon/10.]+=1.;
	  rate_polE[3][realSide][erecon/10.]+=1.;
	} else if (type==1) {
	  UNCORR_rate_polE[side][erecon/10.]+=1.;
	  rate_polE[0][side][erecon/10.]+=1.;
	  rate_polE[1][realSide][erecon/10.]+=1.;
	  rate_polE[2][realSide][erecon/10.]+=1.;
	  rate_polE[3][realSide][erecon/10.]+=1.;
	} else if (type==2) {
          if (side==0) {
	    type = sep.separate23(MWPCEnergy[0]);
	    side = type==2 ? 1 : 0;
	  }
	  else if (side==1) {
	    type = sep.separate23(MWPCEnergy[1]);
	    side = type==2 ? 0 : 1;
	  }
	  if (type==2) {
	    UNCORR_rate_polE[side][erecon/10.]+=1.;
	    rate_polE[0][side][erecon/10.]+=1.;
	    rate_polE[1][side][erecon/10.]+=1.;
	    rate_polE[2][realSide][erecon/10.]+=1.;
	    rate_polE[3][realSide][erecon/10.]+=1.;
	  } else if (type==3) {
	    UNCORR_rate_polE[side][erecon/10.]+=1.;
	    rate_polE[0][side][erecon/10.]+=1.;
	    rate_polE[1][side][erecon/10.]+=1.;
	    rate_polE[2][side][erecon/10.]+=1.;
	    rate_polE[3][realSide][erecon/10.]+=1.;
	  }
	}

	aveThetaTrigg[(int)(erecon/10.)]+=TMath::Abs(returnBeta(primKE)*TMath::Cos(primTheta));
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

      Int_t realSide = TMath::Cos(primTheta)>0.?1:0;
      Int_t side=-1;
      Int_t type=-1;
      if (realSide==0) {
	aveThetaPure[(int)(primKE/10.)]+=TMath::Abs(returnBeta(primKE)*TMath::Cos(primTheta));
	PureEntries[(int)(primKE/10.)]+=1.;
      }
      else if (realSide==1) {
	aveThetaPure[(int)(primKE/10.)]+=TMath::Abs(returnBeta(primKE)*TMath::Cos(primTheta));
	PureEntries[(int)(primKE/10.)]+=1.;
      }

      bool trigger = false;

      EdepQ[0] = (1./(alpha*g_d*g_rest)) * (rand.Poisson(g_rest*rand.Poisson(g_d*rand.Poisson(alpha*EdepQ[0]))));
      EdepQ[1] = (1./(alpha*g_d*g_rest)) * (rand.Poisson(g_rest*rand.Poisson(g_d*rand.Poisson(alpha*EdepQ[1]))));


      bool mwpcEastTrigg = (MWPCEnergy[0]>0.)?true:false;
      bool mwpcWestTrigg = (MWPCEnergy[1]>0.)?true:false;
      bool scintEastTrigg = EdepQ[0]>0. && scintTrigg.triggerE(&EdepQ[0],rand.Rndm());//EdepQ[0]>0.;// 
      bool scintWestTrigg = EdepQ[1]>0. && scintTrigg.triggerW(&EdepQ[1],rand.Rndm());//

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

	if (type==0) {
	  UNCORR_rate_polW[side][erecon/10.]+=1.;
	  rate_polW[0][realSide][erecon/10.]+=1.;
	  rate_polW[1][realSide][erecon/10.]+=1.;
	  rate_polW[2][realSide][erecon/10.]+=1.;
	  rate_polW[3][realSide][erecon/10.]+=1.;
	} else if (type==1) {
	  UNCORR_rate_polW[side][erecon/10.]+=1.;
	  rate_polW[0][side][erecon/10.]+=1.;
	  rate_polW[1][realSide][erecon/10.]+=1.;
	  rate_polW[2][realSide][erecon/10.]+=1.;
	  rate_polW[3][realSide][erecon/10.]+=1.;
	} else if (type==2) {
          if (side==0) {
	    type = sep.separate23(MWPCEnergy[0]);
	    side = type==2 ? 1 : 0;
	  }
	  else if (side==1) {
	    type = sep.separate23(MWPCEnergy[1]);
	    side = type==2 ? 0 : 1;
	  }
	  if (type==2) {
	    UNCORR_rate_polW[side][erecon/10.]+=1.;
	    rate_polW[0][side][erecon/10.]+=1.;
	    rate_polW[1][side][erecon/10.]+=1.;
	    rate_polW[2][realSide][erecon/10.]+=1.;
	    rate_polW[3][realSide][erecon/10.]+=1.;
	  } else if (type==3) {
	    UNCORR_rate_polW[side][erecon/10.]+=1.;
	    rate_polW[0][side][erecon/10.]+=1.;
	    rate_polW[1][side][erecon/10.]+=1.;
	    rate_polW[2][side][erecon/10.]+=1.;
	    rate_polW[3][realSide][erecon/10.]+=1.;
	  }
	}
      
	aveThetaTrigg[(int)(erecon/10.)]+=TMath::Abs(returnBeta(primKE)*TMath::Cos(primTheta));
	TriggEntries[(int)(erecon/10.)]+=1.;
      }
      
    }  
    delete f;
    
  }


  std::vector <Double_t> asymm(100,0.);
  std::vector <Double_t> asymmErr(100,0.);

  std::vector < std::vector<Double_t> > bsCorr(4,std::vector<Double_t>(100,1.));
  
  for (UInt_t b=0;b<100;++b) {

    Double_t uncorrAsymm = 0.;

    Double_t R = (UNCORR_rate_polW[0][b]*UNCORR_rate_polE[1][b]) / (UNCORR_rate_polE[0][b]*UNCORR_rate_polW[1][b]);
    /*Double_t deltaR = TMath::Sqrt( TMath::Power( UNCORR_rate_polW[0]_err[b]*UNCORR_rate_polE[1][b]/(UNCORR_rate_polE[0][b]*UNCORR_rate_polW[1][b]), 2 ) + 
				   TMath::Power( UNCORR_rate_polE[1]_err[b]*UNCORR_rate_polW[0][b]/(UNCORR_rate_polE[0][b]*UNCORR_rate_polW[1][b]), 2 ) +
				   TMath::Power( UNCORR_rate_polE[0]_err[b]*(UNCORR_rate_polW[0][b]*UNCORR_rate_polE[1][b]) / (UNCORR_rate_polE[0][b]*UNCORR_rate_polE[0][b]*UNCORR_rate_polW[1][b]), 2 ) +
				   TMath::Power( UNCORR_rate_polW[1]_err[b]*(UNCORR_rate_polW[0][b]*UNCORR_rate_polE[1][b]) / (UNCORR_rate_polE[0][b]*UNCORR_rate_polW[1][b]*UNCORR_rate_polW[1][b]), 2 ) );*/
    
    
    uncorrAsymm = (1.-TMath::Sqrt(R))/(1.+TMath::Sqrt(R));
    //asymmErr[j] = (deltaR)/(TMath::Sqrt(R)*TMath::Power((TMath::Sqrt(R)+1.),2));
    
    std::vector <Double_t> asymm(4,0.);

    R=0.;

    for (int j=0;j<4;++j) {

      Double_t R = (rate_polW[j][0][b]*rate_polE[j][1][b]) / (rate_polE[j][0][b]*rate_polW[j][1][b]);

    
    asymm[j] = (1.-TMath::Sqrt(R))/(1.+TMath::Sqrt(R));
    //asymmErr[j] = (deltaR)/(TMath::Sqrt(R)*TMath::Power((TMath::Sqrt(R)+1.),2));
    
    }

    bsCorr[0][b] = asymm[0]/uncorrAsymm-1.;
    bsCorr[1][b] = asymm[1]/asymm[0]-1.;
    bsCorr[2][b] = asymm[2]/asymm[1]-1.;
    bsCorr[3][b] = asymm[3]/asymm[2]-1.;
  }
  
 
  for (UInt_t b=0; b<100; ++b) {
    aveThetaPure[b]= ( PureEntries[b]>0.?aveThetaPure[b]/PureEntries[b]:0. );
    aveThetaTrigg[b]= ( TriggEntries[b]>0.?aveThetaTrigg[b]/TriggEntries[b]:0. );
    //aveThetaPure[b]= ( PureEntries[b]>0.?aveThetaPure[b]/PureEntries[b]:0. );
    //aveThetaTrigg[b]= ( TriggEntries[b]>0.?aveThetaTrigg[b]/TriggEntries[b]:0. );
  }
  
  
  std::ofstream ofile(TString::Format("angleCorr_%s.txt",year.Data()).Data());
  ofile << std::setprecision(10);

  for (int b=0;b<100;++b) {
    double ratio = (returnBeta(10.*b+5.)/2.) / aveThetaTrigg[b];
    //double deltaEast = aveThetaTriggE[b]/aveThetaPureE[b];
    //double deltaWest = aveThetaTriggW[b]/aveThetaPureW[b];

    ofile << 10.*b+5. << "\t" 
	  << ratio << "\t" 
	  << ratio-1. << std::endl;
  }
  ofile.close();

  for (int j=0;j<4;++j) {
    ofile.open(TString::Format("%s_deltaBS_%i.txt",year.Data(),j).Data());
    ofile << std::setprecision(10);

    for (int b=0;b<100;++b) {

      ofile << 10.*b+5. << "\t" 
	    << bsCorr[j][b] << std::endl;
    }
    
    ofile.close();
  }

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

