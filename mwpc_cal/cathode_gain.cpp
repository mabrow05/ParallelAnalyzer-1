#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <sstream>

// ROOT libraries
#include <TF1.h>
#include <TH1D.h>
#include <TMath.h>
#include <TString.h>

#include "runInfo.h"
#include "DataTree.hh"
#include "MWPCPositionResponse.hh"
#include "peaks.hh"

std::vector <int>  readOctetFile(int octet) {

  std::vector <int> runs;
  
  char fileName[200];
  sprintf(fileName,"%s/All_Octets/octet_list_%i.dat",getenv("OCTET_LIST"),octet);
  std::ifstream infile(fileName);
  std::string runTypeHold;
  int runNumberHold;
  int numRuns = 0;
  // Populate the map with the runType and runNumber from this octet
  while (infile >> runTypeHold >> runNumberHold) {

    if (runTypeHold=="A2" || runTypeHold=="A5" || runTypeHold=="A7" || runTypeHold=="A10" || 
	runTypeHold=="B2" || runTypeHold=="B5" || runTypeHold=="B7" || runTypeHold=="B10" )  {
     
      runs.push_back(runNumberHold);

    }
    numRuns++;
  }

  infile.close();
 
  std::cout << "Read in octet file for octet " << octet << "\n";
  return runs;

};


int main(int argc, char *argv[])
{
  cout.setf(ios::fixed, ios::floatfield);
  cout.precision(12);

  int octetNumber = atoi(argv[1]);
 
  cout << "Octet Number = " << octetNumber << endl;

  std::vector <int> betaRuns = readOctetFile( octetNumber );

  // Output root file
  TString path = TString::Format("%s/gain_cathodes_octet_%i.root",getenv("GAIN_CATHODES"),octetNumber);
  TFile *fileOut = new TFile(path,"RECREATE");


  // Output histograms
  TH1D *hisEX[16];
  TH1D *hisEY[16];
  TH1D *hisWX[16];
  TH1D *hisWY[16];

  for ( int i = 0; i<16; ++i ) {

    hisEX[i] = new TH1D(TString::Format("hisEX%i",i),TString::Format("hisEX%i",i),120,0.,6.);
    hisEY[i] = new TH1D(TString::Format("hisEY%i",i),TString::Format("hisEY%i",i),120,0.,6.);
    hisWX[i] = new TH1D(TString::Format("hisWX%i",i),TString::Format("hisWX%i",i),120,0.,6.);
    hisWY[i] = new TH1D(TString::Format("hisWY%i",i),TString::Format("hisWY%i",i),120,0.,6.);

  }

  MWPCCathodeHandler cathResp;
  double wireSep = TMath::Sqrt(0.6)*10.16;
	
  
  //Loop over all runs in octet
  
  for ( auto rn : betaRuns ) {
    
    // DataTree structure
    DataTree *t = new DataTree();

    // Input ntuple
    char tempIn[500];
    TString infile = TString::Format("%s/replay_pass3_%i.root", getenv("REPLAY_PASS3"),rn);
    t->setupInputTree(std::string(infile.Data()),"pass3");
 
    int nEvents = t->getEntries();
    cout << "... Processing nEvents = " << nEvents << endl;

    // Loop over all electron events 
    for ( int evt=0; evt<nEvents; ++evt ) {

      t->getEvent(evt);

      if ( evt%10000==0 ) std::cout << evt << std::endl;
      
      if ( t->PID != 1 || t->Type != 0 ) continue; // only use type 0 electrons

      // East side first
      if ( t->Side == 0 && t->xE.nClipped==0 && t->yE.nClipped==0 && t->Erecon>100. ) { // Go above 
      
	// EX plane first
	int wire = 0;

	while ( TMath::Abs( t->xE.center - TMath::Sqrt(0.6)*cathResp.getWirePosEX(wire) )  > wireSep/2. ) wire++;
	hisEX[wire]->Fill( t->Cathodes_Ex[wire]/t->AnodeE );

	// EY plane
	wire = 0;
	
	while ( TMath::Abs( t->yE.center - TMath::Sqrt(0.6)*cathResp.getWirePosEY(wire) )  > wireSep/2. ) wire++;
	hisEY[wire]->Fill( t->Cathodes_Ey[wire]/t->AnodeE );
      }

      // Now the west side
      if ( t->Side == 1 && t->xW.nClipped==0 && t->yW.nClipped==0 && t->Erecon>100. ) {
      
	// WX plane first
	int wire = 0;

	while ( TMath::Abs( t->xW.center - TMath::Sqrt(0.6)*cathResp.getWirePosWX(wire) )  > wireSep/2. ) wire++;
	hisWX[wire]->Fill( t->Cathodes_Wx[wire]/t->AnodeW );

	// WY plane
	wire = 0;
	
	while ( TMath::Abs( t->yW.center - TMath::Sqrt(0.6)*cathResp.getWirePosWY(wire) )  > wireSep/2. ) wire++;
	hisWY[wire]->Fill( t->Cathodes_Wy[wire]/t->AnodeW );
      }

    }

    delete t;
  }

  /////////////////////////////////////////////////////////////////

  // Now we want to fit all of the histograms with a gaussian

  std::vector <double> mean_EX(16,0.);
  std::vector <double> mean_EY(16,0.);
  std::vector <double> mean_WX(16,0.);
  std::vector <double> mean_WY(16,0.);
  
  SinglePeakHist *p;
  
  for ( int h = 0; h<16; ++h ) {

    if ( hisEX[h]->GetEntries() > 200 ) {
      p = new SinglePeakHist(hisEX[h],0., 6., true, 5, 0.75,1.25);
      mean_EX[h] = p->ReturnMean();
      delete p;
    }
    else mean_EX[h] = hisEX[h]->GetMean();

    if ( hisEY[h]->GetEntries() > 200 ) {
      p = new SinglePeakHist(hisEY[h],0., 6., true, 5, 0.75,1.25);
      mean_EY[h] = p->ReturnMean();
      delete p;
    }
    else mean_EY[h] = hisEY[h]->GetMean();

    if ( hisWX[h]->GetEntries() > 200 ) {
      p = new SinglePeakHist(hisWX[h],0., 6., true, 5, 0.75,1.25);
      mean_WX[h] = p->ReturnMean();
      delete p;
    }
    else mean_WX[h] = hisWX[h]->GetMean();

    if ( hisWY[h]->GetEntries() > 200 ) {
      p = new SinglePeakHist(hisWY[h],0., 6., true, 5, 0.75,1.25);
      mean_WY[h] = p->ReturnMean();
      delete p;
    }
    else mean_WY[h] = hisWY[h]->GetMean();

  }

  fileOut->Write();
  delete fileOut;

  // Now normalize to wire 7
  std::vector <double> gain_EX(16,0.);
  std::vector <double> gain_EY(16,0.);
  std::vector <double> gain_WX(16,0.);
  std::vector <double> gain_WY(16,0.);
  
  for ( int h = 0; h<16; ++h ) {
    
    gain_EX[h] = mean_EX[h] / mean_EX[7];
    gain_EY[h] = mean_EY[h] / mean_EY[7];
    gain_WX[h] = mean_WX[h] / mean_WX[7];
    gain_WY[h] = mean_WY[h] / mean_WY[7];
    
  }

  std::ofstream gainFile(TString::Format("%s/gain_cathodes_octet_%i.dat",getenv("GAIN_CATHODES"),octetNumber).Data());

  gainFile << "#EX" << std::setw(12)
	   << "EY"  << std::setw(12)
	   << "WX"  << std::setw(12)
	   << "WY\n";
  gainFile << std::setprecision(7);
  
  for ( int h = 0; h<16; ++h ) {
    
    gainFile << gain_EX[h] << std::setw(12)
	     << gain_EY[h] << std::setw(12)
	     << gain_WX[h] << std::setw(12)
	     << gain_WY[h] ;
    
    if ( h<15 ) gainFile << std::endl;
    
  }

  gainFile.close();
  
  return 0;
}
