#include "BetaDecayTools.hh"
#include "DataTree.hh"
#include <TString.h>
#include <TH1D.h>
#include <TGraphErrors.h>
#include <TMath.h>

#include <iostream>
#include <fstream>


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

int main(int argc, char *argv[]) {

  
  //Show that we can do the Kurie fit...

  int octet = atoi(argv[1]);

  int nBins = 100;
  double minRange = 0.;
  double maxRange = 1000.;
  
  std::vector <int> runs = readOctetFile(octet);

  
  // Vectors for creating the individual runs events and errors
  std::vector < std::vector < std::vector <Double_t> > > pmtSpec;
  std::vector < std::vector < std::vector <Double_t> > > pmtSpecErr;
    
  pmtSpec.resize(8,std::vector<std::vector<Double_t>>(8,std::vector<Double_t>(nBins,0.)));
  pmtSpecErr.resize(8,std::vector<std::vector<Double_t>>(8,std::vector<Double_t>(nBins,0.)));

  // Now loop over each run to determine the individual spectra, then fill their appropriate bins in the vector

  TH1D *spec[8]; // All 8 PMTs signals

  int nRun = 0;
  
  for ( auto rn : runs ) {

    for ( int i=0; i<8; ++i ) {
      spec[i] = new TH1D(TString::Format("PMT%i",i),TString::Format("PMT %i",i),
		      nBins, minRange, maxRange);
    }
    
    // DataTree structure
    DataTree t;

    // Input ntuple
    char tempIn[500];
    sprintf(tempIn, "%s/replay_pass3_%i.root", getenv("REPLAY_PASS3"),rn);
    
    t.setupInputTree(std::string(tempIn),"pass3");

    unsigned int nevents = t.getEntries();

    //t.getEvent(nevents-1);
    //totalTimeOFF += t.Time;
    //totalBLINDTimeOFF[0] += t.TimeE;
    //totalBLINDTimeOFF[1] += t.TimeW;

    double r2E = 0.; //position of event squared
    double r2W = 0.;
    
    for (unsigned int n=0 ; n<nevents ; n++ ) {

      t.getEvent(n);

      if ( t.PID==1 && t.Side<2 && t.Type==0 && t.Erecon>0. ) {

	if ( t.Side==0 ) {
	  spec[0]->Fill(t.ScintE.e1);
	  spec[1]->Fill(t.ScintE.e2);
	  spec[2]->Fill(t.ScintE.e3);
	  spec[3]->Fill(t.ScintE.e4);
	}

	if ( t.Side==1 ) {
	  spec[4]->Fill(t.ScintW.e1);
	  spec[5]->Fill(t.ScintW.e2);
	  spec[6]->Fill(t.ScintW.e3);
	  spec[7]->Fill(t.ScintW.e4);
	}

      }
    }

    for ( int p=0; p<8; ++p ) {
      for ( int bin=1; bin<=nBins; ++bin ) {
	pmtSpec[nRun][p][bin-1] = (double)spec[p]->GetBinContent(bin);
	pmtSpecErr[nRun][p][bin-1] = spec[p]->GetBinError(bin);
      }
    }

    for ( int i=0; i<8; ++i ) delete spec[i];
    nRun++;
    
  }


  TFile *fout = new TFile("testingKurie.root","RECREATE");
  
  //TH1D *spec[8]; // All 8 PMTs signals
  TGraphErrors kurie[8]; //These will hold the kurie plots

  for ( int i=0; i<8; ++i ) {

    spec[i] = new TH1D(TString::Format("PMT%i",i),TString::Format("PMT %i",i),
		    nBins, minRange, maxRange);

  }

  //Now loop over all of the runs and create weighted averages in every bin for each
  // PMT

  
  for ( int p=0; p<8; ++p ) {
    for ( int bin=1; bin<=nBins; ++bin ) {

      double numer = 0.;
      double denom = 0.;
      
      for ( int i=0; i<nRun; ++i ) {
	double weight = pmtSpec[i][p][bin]>0. ? 1./(pmtSpecErr[i][p][bin]*pmtSpecErr[i][p][bin]) : 1.;
	numer += pmtSpec[i][p][bin]*weight;
	denom += weight;
      }

      spec[p]->SetBinContent(bin,numer/denom);
      spec[p]->SetBinError(bin,TMath::Sqrt(1./denom));
    }
  }

  //Now do some Kurie Fitting

  KurieFitter kf;

  kf.FitSpectrum(spec[0],150., 635., 0.9);
  
  kurie[0] = kf.returnKuriePlot();

  fout->Write();
  delete fout;

  return 0;
  
}

  
