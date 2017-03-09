/* 
    Determine the 2/3 separation cut
*/


#include "MBUtils.hh"

#include <vector>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <fstream>

#include <TRandom3.h>
#include <TTree.h>
#include <TLeaf.h>
#include <TFile.h>
#include <TH1F.h>
#include <TF1.h>
#include <TH1D.h>
#include <TChain.h>
#include <TMath.h>
#include <TSQLRow.h>
#include <TString.h>
#include <TStyle.h>
#include <TBranch.h>
#include <TLeaf.h>


using namespace std;

std::vector <Int_t> badOct = {7,9,59,60,61,62,63,64,65,66,70,92}; 


std::vector<Int_t> readOctetFile(int octet) {

  std::vector<Int_t> runs;
  
  char fileName[200];
  sprintf(fileName,"%s/All_Octets/octet_list_%i.dat",getenv("OCTET_LIST"),octet);
  std::ifstream infile(fileName);
  std::string runTypeHold;
  int runNumberHold;
  
  // Populate the map with the runType and runNumber from this octet
  while (infile >> runTypeHold >> runNumberHold) {

    if ( runTypeHold=="A2" || runTypeHold=="A10" || runTypeHold=="B5" || runTypeHold=="B7"
	 || runTypeHold=="B2" || runTypeHold=="B10" || runTypeHold=="A5" || runTypeHold=="A7" )  {
      runs.push_back(runNumberHold);
    }

  }

  infile.close();
 
  std::cout << "Read in octet file for octet " << octet << "\n";
  return runs;
};



void SeparateBackscatters(int octetMin, int octetMax) 
{

  gStyle->SetOptStat(0);
  double fiducialCut = 50.; //mm
  std::vector<Int_t> betaRuns;
 
  cout << "Calculating FG spectra..." << endl;

  //Load all beta decay runs
  
  for ( int i=octetMin ; i<=octetMax ; i++ ) {

    if ( std::find(badOct.begin(), badOct.end(),i) != badOct.end() ) continue;  //Checking if octet should be ignored for data quality reasons
    std::vector<Int_t> octRuns = readOctetFile(i);
    betaRuns.insert(betaRuns.end(),octRuns.begin(),octRuns.end());

  }
   
  char temp[200];

  //Root file to output to
  TFile *outfile = new TFile(TString::Format("Backscatters_%i-%i.root",octetMin,octetMax),"RECREATE");

  //Make all the pertinent histograms

  TH1D *hType2E = new TH1D("hType2E",TString::Format("Type II East Octets %i-%i E_{MWPC}",octetMin,octetMax),100, 0., 20.);
  TH1D *hType3E = new TH1D("hType3E",TString::Format("Type III East Octets %i-%i E_{MWPC}",octetMin,octetMax),100, 0., 20.);
  TH1D *hType2W = new TH1D("hType2W",TString::Format("Type II West Octets %i-%i E_{MWPC}",octetMin,octetMax),100, 0., 20.);
  TH1D *hType3W = new TH1D("hType3W",TString::Format("Type III West Octets %i-%i E_{MWPC}",octetMin,octetMax),100, 0., 20.);
  
 
  // Separated type 2
  TH1D* hScint2E = new TH1D("hScint2E",TString::Format("Type II East Octets %i-%i E_{recon}",octetMin,octetMax),120,0.,1200.);
  TH1D* hScint3E = new TH1D("hScint3E",TString::Format("Type III East Octets %i-%i E_{recon}",octetMin,octetMax),120,0.,1200.);
  TH1D* hScint2W = new TH1D("hScint2W",TString::Format("Type II West Octets %i-%i E_{recon}",octetMin,octetMax),120,0.,1200.);
  TH1D* hScint3W = new TH1D("hScint3W",TString::Format("Type III West Octets %i-%i E_{recon}",octetMin,octetMax),120,0.,1200.);


  //Process all runs
  for ( auto rn : betaRuns ) {
    
    std::string infile;
    TFile *input;
    TTree *Tin;
    
    sprintf(temp,"revCalSim_%i_Beta.root",rn);
    infile = getenv("REVCALSIM")+std::string("/beta/")+std::string(temp);
    input = new TFile(infile.c_str(), "READ");
    Tin = (TTree*)input->Get("revCalSim");

    double mwpcEX, mwpcEY, mwpcWX, mwpcWY, Erecon;
    double EmwpcE=0., EmwpcW=0.;
    int PID, type, side; // basic analysis tags
    double primTheta;
    int hitCountSD[24];

    Tin->SetBranchAddress("primTheta",&primTheta);
    Tin->SetBranchAddress("PID", &PID);
    Tin->SetBranchAddress("type", &type);
    Tin->SetBranchAddress("side", &side); 
    Tin->SetBranchAddress("Erecon",&Erecon);
    Tin->GetBranch("MWPCEnergy")->GetLeaf("MWPCEnergyE")->SetAddress(&EmwpcE);
    Tin->GetBranch("MWPCEnergy")->GetLeaf("MWPCEnergyW")->SetAddress(&EmwpcW);
    Tin->GetBranch("xE")->GetLeaf("center")->SetAddress(&mwpcEX);
    Tin->GetBranch("yE")->GetLeaf("center")->SetAddress(&mwpcEY);
    Tin->GetBranch("xW")->GetLeaf("center")->SetAddress(&mwpcWX);
    Tin->GetBranch("yW")->GetLeaf("center")->SetAddress(&mwpcWY);
    Tin->SetBranchAddress("hitCountSD",hitCountSD);

    

    unsigned int nevents = Tin->GetEntriesFast();

    double r2E = 0.; //position of event squared
    double r2W = 0.;

    for (unsigned int n=0 ; n<nevents ; n++ ) {
      
      
      Tin->GetEvent(n);

      if ( PID==1 && side<2 && type==2 && Erecon>0.) {
	
	//***********************************************************************************************************************
	// Filling rate histograms with "good" events to calculate the corrections
	
	r2E = mwpcEX*mwpcEX + mwpcEY*mwpcEY;
	r2W = mwpcWX*mwpcWX + mwpcWY*mwpcWY;

	if ( r2E<(fiducialCut*fiducialCut) && r2W<(fiducialCut*fiducialCut) )	  {
		
	  if (side==0) {
	    // Type 2 goes west then triggers east for east side
	    if ( primTheta < TMath::Pi()/2. && hitCountSD[5]==1 ) {
	      hType2E->Fill(EmwpcE);
	      hScint2E->Fill(Erecon);
	    }
	    //Type 3 goes east and triggers, then goes west
	    else if ( primTheta > TMath::Pi()/2. && hitCountSD[15]==1 ) {
	      hType3E->Fill(EmwpcE);
	      hScint3E->Fill(Erecon);
	    }
	  }
	  else if (side==1) {
	    // Type 2 goes east then triggers west for west side
	    if ( primTheta > TMath::Pi()/2. && hitCountSD[15]==1 ) {
	      hType2W->Fill(EmwpcW);
	      hScint2W->Fill(Erecon);
	    }
	    //Type 3 goes west and triggers, then goes east
	    else if ( primTheta < TMath::Pi()/2. && hitCountSD[5]==1 ){
	      hType3W->Fill(EmwpcW);
	      hScint3W->Fill(Erecon);
	    }
	  }
	}
      }
    }

    input->Close();
    if (input) delete input;
    cout << "Finished Run " << rn << endl;
  }

  outfile->Write();
  delete outfile;
  cout << endl;
    
};
  

int main(int argc, char *argv[]) {

  int octetMin = atoi(argv[1]);
  int octetMax = atoi(argv[2]);

  SeparateBackscatters(octetMin, octetMax);
  //writeRatesToFile(octetMin, octetMax);

  //tests
  /*UInt_t XePeriod = getXeRunPeriod(atoi(argv[1]));
  vector < vector <Double_t> > triggerFunc = getTriggerFunctionParams(XePeriod,7);
  Double_t triggProbE = triggerProbability(triggerFunc[0],25.);
  Double_t triggProbW = triggerProbability(triggerFunc[1],25.);
  cout << triggProbE << " " << triggProbW << endl;*/


}
  
