#include <TString.h>
#include <vector>
#include "calibrationTools.hh"
#include <fstream>
#include <iomanip>
#include <iostream>

#include <TTree.h>
#include <TBranch.h>
#include <TLeaf.h>
#include <TFile.h>

std::vector <Int_t> badOct {7,9,59,60,61,62,63,64,65,66,67,91,93,101,107,121}; 
std::vector <TString> anaChoices {"A","B","C","D","E","F","G","H","J","K"};

// Function to return beta when given the kinetic energy of an electron
Double_t returnBeta(Double_t En) { 
  Double_t me = 510.998928; //rest mass energy of electron in keV
  return sqrt(En*En+2.*En*me)/(En+me);
};


int main(int argc, char* argv[]) {
  
  TString year = TString(argv[1]);
  
  int octetNumStart = year==TString("2011-2012")?0:60;
  int octetNumEnd = year==TString("2011-2012")?59:121;

   BackscatterSeparator sep;
  

  std::vector <Int_t> runs;

  //Need to read in list of all beta decay runs
  for (int octetNum=octetNumStart; octetNum<octetNumEnd+1; octetNum++) {
    
    if (std::find(badOct.begin(), badOct.end(),octetNum) != badOct.end()) { continue; } //Checking if octet should be ignored for data quality reasons
    
    std::ifstream octFile(TString::Format("%s/All_Octets/octet_list_%i.dat",
					  getenv("OCTET_LIST"),octetNum));
    Int_t runNumber = 0;
    std::string runType = "";
    
    while ( octFile >> runType >> runNumber ) {
      if ( runType == "A2" || runType == "A10" || runType == "B5" || runType == "B7" 
	   || runType == "B2" || runType == "B10" || runType == "A5" || runType == "A7" ) {
	runs.push_back(runNumber);
      }
    }
  }
	

  // Vectors to hold the BetaCosTheta in every energy bin
  // and the total number of entries
  std::vector <std::vector<Double_t> > aveBetaCosTheta(anaChoices.size(),std::vector<Double_t>(1200,0.));
  std::vector <std::vector<Double_t> > counts(anaChoices.size(),std::vector<Double_t>(1200,0.));

  // loop over all files of these runs, fill vectors

  TString infile;
  TFile *input;
  TTree *Tin;
  
  double TimeE=0., TimeW=0., Erecon=0., MWPCEnergyE=0., MWPCEnergyW=0.; //Branch Variables being read in
  int PID, Side, Type;
  //double scintPosE[3]={0.};
  //double scintPosW[3]={0.};
  //double cathRespPosE[3]={0.}; //Position reconstructed the same way it is in data
  //double cathRespPosW[3]={0.};
  //double mwpcPosE[3]={0.};
  //double mwpcPosW[3]={0.}; //holds the position of the event in the MWPC for simulated data
  double EmwpcX=0., EmwpcY=0., WmwpcX=0., WmwpcY=0.;
  int xE_nClipped, yE_nClipped, xW_nClipped, yW_nClipped;
  int xE_mult=0, yE_mult=0, xW_mult=0, yW_mult=0;

  double primKE, primTheta;
  
  double r2E,r2W;

  UInt_t nevents;
  for (auto run : runs) {
    std::cout << "run number " << run << "\n";
    infile = TString::Format("%s/beta_highStatistics/revCalSim_%i_Beta.root",getenv("REVCALSIM"),run);
    input = new TFile(infile.Data(), "READ");
    Tin = (TTree*)input->Get("revCalSim");
    
    Tin->SetBranchAddress("PID", &PID);
    Tin->SetBranchAddress("type", &Type);
    Tin->SetBranchAddress("side", &Side); 
    Tin->SetBranchAddress("Erecon",&Erecon);
    //Tin->GetBranch("time")->GetLeaf("timeE")->SetAddress(&TimeE);
    //Tin->GetBranch("time")->GetLeaf("timeW")->SetAddress(&TimeW);
    //Tin->GetBranch("ScintPosAdjusted")->GetLeaf("ScintPosAdjE")->SetAddress(scintPosE);
    //Tin->GetBranch("ScintPosAdjusted")->GetLeaf("ScintPosAdjW")->SetAddress(scintPosW);
    //Tin->GetBranch("cathRespPos")->GetLeaf("cathRespPosE")->SetAddress(cathRespPosE);
    //Tin->GetBranch("cathRespPos")->GetLeaf("cathRespPosW")->SetAddress(cathRespPosW);
    //Tin->GetBranch("MWPCPos")->GetLeaf("MWPCPosE")->SetAddress(mwpcPosE);
    //Tin->GetBranch("MWPCPos")->GetLeaf("MWPCPosW")->SetAddress(mwpcPosW);
    Tin->GetBranch("xE")->GetLeaf("center")->SetAddress(&EmwpcX);
    Tin->GetBranch("yE")->GetLeaf("center")->SetAddress(&EmwpcY);
    Tin->GetBranch("xW")->GetLeaf("center")->SetAddress(&WmwpcX);
    Tin->GetBranch("yW")->GetLeaf("center")->SetAddress(&WmwpcY);
    Tin->GetBranch("xE")->GetLeaf("nClipped")->SetAddress(&xE_nClipped);
    Tin->GetBranch("yE")->GetLeaf("nClipped")->SetAddress(&yE_nClipped);
    Tin->GetBranch("xW")->GetLeaf("nClipped")->SetAddress(&xW_nClipped);
    Tin->GetBranch("yW")->GetLeaf("nClipped")->SetAddress(&yW_nClipped);
    Tin->GetBranch("xE")->GetLeaf("mult")->SetAddress(&xE_mult);
    Tin->GetBranch("yE")->GetLeaf("mult")->SetAddress(&yE_mult);
    Tin->GetBranch("xW")->GetLeaf("mult")->SetAddress(&xW_mult);
    Tin->GetBranch("yW")->GetLeaf("mult")->SetAddress(&yW_mult);
    Tin->GetBranch("MWPCEnergy")->GetLeaf("MWPCEnergyE")->SetAddress(&MWPCEnergyE);
    Tin->GetBranch("MWPCEnergy")->GetLeaf("MWPCEnergyW")->SetAddress(&MWPCEnergyW);

    //Tin->SetBranchAddress("primPos",primPos);
    Tin->SetBranchAddress("primKE",&primKE);
    Tin->SetBranchAddress("primTheta",&primTheta);
    
    sep.LoadCutCurve(run);

    nevents = Tin->GetEntriesFast();
    
    //Run over all events in file, fill output histogram with rate
    
    r2E = 0.; //position of event squared
    r2W = 0.;

    
    for (unsigned int i=0; i<nevents; ++i) {
      
      Tin->GetEvent(i);

      if ( Type>3  || Side>1 || Erecon<=0. ) continue;
      
      if (PID==1) {

	//Cut out clipped events
	if ( Type!=0 ) {
	  //if (xE_nClipped>0 || yE_nClipped>0 || xW_nClipped>0 || yW_nClipped>0) continue;
	  if (xE_mult<1 || yE_mult<1 || xW_mult<1 || yW_mult<1) continue;
	}
	else {
	  if ( Side==0 ) {
	    //if ( xE_nClipped>0 || yE_nClipped>0 )  continue;
	    if ( xE_mult<1 || yE_mult<1 )  continue;
	  }
	  else if ( Side==1 ) {
	    //if ( xW_nClipped>0 || yW_nClipped>0 )  continue;
	    if ( xW_mult<1 || yW_mult<1 )  continue;
	  }
	}


	r2E = ( EmwpcX*EmwpcX + EmwpcY*EmwpcY ) ; //Transforming to decay trap coords
	r2W = ( WmwpcX*WmwpcX + WmwpcY*WmwpcY );
	  
	  
	if ( r2E>2500. || r2W>2500. ) continue;


	if (Type==2) {
	  
	  if (Side==0) {
	    Type = sep.separate23(MWPCEnergyE);
	    Side = Type==2 ? 1 : 0;
	  }
	  else if (Side==1) {
	    Type = sep.separate23(MWPCEnergyW);
	    Side = Type==2 ? 0 : 1;
	  }
	}
	

	// Fill histograms according to event types in this analysis choice
	
	//Type0
	if (Type==0 ) { 
	  aveBetaCosTheta[0][(int)(Erecon/10.)]+=TMath::Abs(returnBeta(primKE)*TMath::Cos(primTheta));
	  aveBetaCosTheta[1][(int)(Erecon/10.)]+=TMath::Abs(returnBeta(primKE)*TMath::Cos(primTheta));
	  aveBetaCosTheta[2][(int)(Erecon/10.)]+=TMath::Abs(returnBeta(primKE)*TMath::Cos(primTheta));
	  aveBetaCosTheta[3][(int)(Erecon/10.)]+=TMath::Abs(returnBeta(primKE)*TMath::Cos(primTheta));
	  counts[0][(int)(Erecon/10.)]+=1.;
	  counts[1][(int)(Erecon/10.)]+=1.;
	  counts[2][(int)(Erecon/10.)]+=1.;
	  counts[3][(int)(Erecon/10.)]+=1.;
	}
	//Type1
	else if ( Type==1 ) { 
	  aveBetaCosTheta[0][(int)(Erecon/10.)]+=TMath::Abs(returnBeta(primKE)*TMath::Cos(primTheta));
	  aveBetaCosTheta[1][(int)(Erecon/10.)]+=TMath::Abs(returnBeta(primKE)*TMath::Cos(primTheta));
	  aveBetaCosTheta[2][(int)(Erecon/10.)]+=TMath::Abs(returnBeta(primKE)*TMath::Cos(primTheta));
	  aveBetaCosTheta[5][(int)(Erecon/10.)]+=TMath::Abs(returnBeta(primKE)*TMath::Cos(primTheta));
	  counts[0][(int)(Erecon/10.)]+=1.;
	  counts[1][(int)(Erecon/10.)]+=1.;
	  counts[2][(int)(Erecon/10.)]+=1.;
	  counts[5][(int)(Erecon/10.)]+=1.;
	}
	//Type2
	else if ( Type==2 ) { 
	  aveBetaCosTheta[0][(int)(Erecon/10.)]+=TMath::Abs(returnBeta(primKE)*TMath::Cos(primTheta));
	  aveBetaCosTheta[2][(int)(Erecon/10.)]+=TMath::Abs(returnBeta(primKE)*TMath::Cos(primTheta));
	  aveBetaCosTheta[6][(int)(Erecon/10.)]+=TMath::Abs(returnBeta(primKE)*TMath::Cos(primTheta));
	  aveBetaCosTheta[7][(int)(Erecon/10.)]+=TMath::Abs(returnBeta(primKE)*TMath::Cos(primTheta));
	  aveBetaCosTheta[8][(int)(Erecon/10.)]+=TMath::Abs(returnBeta(primKE)*TMath::Cos(primTheta));
	  counts[0][(int)(Erecon/10.)]+=1.;
	  counts[2][(int)(Erecon/10.)]+=1.;
	  counts[6][(int)(Erecon/10.)]+=1.;
	  counts[7][(int)(Erecon/10.)]+=1.;
	  counts[8][(int)(Erecon/10.)]+=1.;
	}
	//Type3
	else if (Type==3 ) { 
	  aveBetaCosTheta[0][(int)(Erecon/10.)]+=TMath::Abs(returnBeta(primKE)*TMath::Cos(primTheta));
	  aveBetaCosTheta[2][(int)(Erecon/10.)]+=TMath::Abs(returnBeta(primKE)*TMath::Cos(primTheta));
	  aveBetaCosTheta[6][(int)(Erecon/10.)]+=TMath::Abs(returnBeta(primKE)*TMath::Cos(primTheta));
	  aveBetaCosTheta[7][(int)(Erecon/10.)]+=TMath::Abs(returnBeta(primKE)*TMath::Cos(primTheta));
	  aveBetaCosTheta[9][(int)(Erecon/10.)]+=TMath::Abs(returnBeta(primKE)*TMath::Cos(primTheta));
	  counts[0][(int)(Erecon/10.)]+=1.;
	  counts[2][(int)(Erecon/10.)]+=1.;
	  counts[6][(int)(Erecon/10.)]+=1.;
	  counts[7][(int)(Erecon/10.)]+=1.;
	  counts[9][(int)(Erecon/10.)]+=1.;
	} 
	
      }
    }
  
    delete input;

  }

  for (UInt_t i=0;i<anaChoices.size();++i) {
    for (UInt_t j=0;j<counts[0].size();++j) {
      aveBetaCosTheta[i][j]/=counts[i][j];
    }
  }


  //Calculate corrections and write out
  for (UInt_t i=0;i<anaChoices.size();++i) {
    
    std::ofstream ofile(TString::Format("%s_DeltaAngle_anaCh%s.txt",year.Data(),anaChoices[i].Data()));
    ofile << std::setprecision(15);
    for (UInt_t j=0;j<counts[0].size();++j) {
      Double_t en = j*10.+5.;
      Double_t corr = 0.5*returnBeta(en)/aveBetaCosTheta[i][j]-1.;
      ofile << en << "\t" << aveBetaCosTheta[i][j] << "\t" << corr << "\n";
    }
    ofile.close();
  }
}
