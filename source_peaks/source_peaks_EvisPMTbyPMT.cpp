#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>

#include <sstream>

// ROOT libraries
#include "TRandom3.h"
#include <TTree.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TF1.h>
#include <TH1D.h>

#include "fullTreeVariables.h"
#include "MWPCGeometry.h"
#include "pedestals.h"
#include "cuts.h"
#include "basic_reconstruction.h"

#include "replay_pass2.h"
#include "replay_pass3.h"
#include "replay_pass4.h"
#include "sourcePeaks.h"
#include "DataTree.hh"

using namespace std;

int main(int argc, char *argv[])
{
  //cout.setf(ios::fixed, ios::floatfield);
  //cout.precision(12);

  // Run number integer
  istringstream ss(argv[1]);
  int runNumber;
  ss >> runNumber;
  cout << "runNumber = " << runNumber << endl;

  // Read source list
  int nSources = 0;
  string sourceName[3];
  char tempList[500];
  sprintf(tempList, "%s/source_list_%s.dat",getenv("SOURCE_LIST"), argv[1]);
  //cout << tempList << endl;
  ifstream fileList(tempList);
  if (fileList.is_open()) fileList >> nSources;
  cout << " ... Number of sources: " << nSources << endl;
  for (int n=0; n<nSources; n++) {
    fileList >> sourceName[n];
    cout << "  " << sourceName[n] << endl;
  }

  //Checking if one of the sources is Bi so we can add in that we will fit the second Bi peak as well
  bool useLowBiPeak=false;
  int BiPeakIndex = 0;
  for (int n=0; n<nSources; n++) {
    if (sourceName[n]=="Bi") {
      useLowBiPeak=true;
      BiPeakIndex=n;
      continue;
    }
  }
      

  //Adding in here that we are only looking for runs with Ce, Sn, or Bi in them
  bool correctSource = false;
  bool useSource[3] = {false,false,false};
  for (int n=0; n<nSources; n++) {
    if (sourceName[n]=="Ce" || sourceName[n]=="Sn" || sourceName[n]=="Bi" || sourceName[n]=="In") {
      correctSource = true;
      useSource[n]=true;
      continue; 
    }
  }
  
  if (!correctSource) {
    cout << "Source run with no Ce, Sn, or Bi\n";
    exit(0);
  }

  // Open output ntuple
  char tempOut[500];
  sprintf(tempOut, "%s/source_peaks_EvisPMTbyPMT_%s.root",getenv("SOURCE_PEAKS"), argv[1]);
  TFile *fileOut = new TFile(tempOut,"RECREATE");


  // Quenched energies from simulation
  double EQ[3];
  /*for (int n=0; n<nSources; n++) {
    if (sourceName[n] == "Ce") EQ[n] = peakCe;//98.2;
    else if (sourceName[n] == "Sn") EQ[n] = peakSn;//331.2;
    else if (sourceName[n] == "Bi") EQ[n] = peakBiHigh;//928.0;
    else {cout << "Unidentified source" << endl; EQ[n] = 1.;} // This is a place holder // exit(0);}
    }*/

  // Read source positions
  double xEast[3], yEast[3], sigmaEast[3];
  double xWest[3], yWest[3], sigmaWest[3];
  char tempPos[500];
  sprintf(tempPos, "%s/source_positions_%s.dat",getenv("SOURCE_POSITIONS"), argv[1]);
  ifstream filePos(tempPos);
  cout << " ... Reading source positions" << endl;
  for (int n=0; n<nSources; n++) {
    filePos >> xEast[n] >> yEast[n] >> sigmaEast[n]
            >> xWest[n] >> yWest[n] >> sigmaWest[n];
    cout << xEast[n] << " " << yEast[n] << " " << sigmaEast[n] << " "
         << xWest[n] << " " << yWest[n] << " " << sigmaWest[n] << endl;
  }

  // Output histograms
  int nBin = 400;

  TH1F *his[3][8];
  his[0][0] = new TH1F("his1_E0", "", nBin,0.0,1200.0);
  his[0][1] = new TH1F("his1_E1", "", nBin,0.0,1200.0);
  his[0][2] = new TH1F("his1_E2", "", nBin,0.0,1200.0);
  his[0][3] = new TH1F("his1_E3", "", nBin,0.0,1200.0);
  his[0][4] = new TH1F("his1_W0", "", nBin,0.0,1200.0);
  his[0][5] = new TH1F("his1_W1", "", nBin,0.0,1200.0);
  his[0][6] = new TH1F("his1_W2", "", nBin,0.0,1200.0);
  his[0][7] = new TH1F("his1_W3", "", nBin,0.0,1200.0);

  his[1][0] = new TH1F("his2_E0", "", nBin,0.0,1200.0);
  his[1][1] = new TH1F("his2_E1", "", nBin,0.0,1200.0);
  his[1][2] = new TH1F("his2_E2", "", nBin,0.0,1200.0);
  his[1][3] = new TH1F("his2_E3", "", nBin,0.0,1200.0);
  his[1][4] = new TH1F("his2_W0", "", nBin,0.0,1200.0);
  his[1][5] = new TH1F("his2_W1", "", nBin,0.0,1200.0);
  his[1][6] = new TH1F("his2_W2", "", nBin,0.0,1200.0);
  his[1][7] = new TH1F("his2_W3", "", nBin,0.0,1200.0);

  his[2][0] = new TH1F("his3_E0", "", nBin,0.0,1200.0);
  his[2][1] = new TH1F("his3_E1", "", nBin,0.0,1200.0);
  his[2][2] = new TH1F("his3_E2", "", nBin,0.0,1200.0);
  his[2][3] = new TH1F("his3_E3", "", nBin,0.0,1200.0);
  his[2][4] = new TH1F("his3_W0", "", nBin,0.0,1200.0);
  his[2][5] = new TH1F("his3_W1", "", nBin,0.0,1200.0);
  his[2][6] = new TH1F("his3_W2", "", nBin,0.0,1200.0);
  his[2][7] = new TH1F("his3_W3", "", nBin,0.0,1200.0);


  // Open input ntuple
  char tempIn[500];
  sprintf(tempIn, "%s/replay_pass4_%s.root",getenv("REPLAY_PASS4"), argv[1]);
  DataTree *t = new DataTree();
  t->setupInputTree(std::string(tempIn),"pass4");

  int nEvents = t->getEntries();
  cout << "Processing " << argv[1] << " ... " << endl;
  cout << "... nEvents = " << nEvents << endl;

  // Loop over events
  for (int i=0; i<nEvents; i++) {
    t->getEvent(i);

    // Use Type 0 events
    if (t->Type != 0 || t->PID!=1) continue;

    if (useSource[0]) {

      if (t->Side == 0) {
	if ( (t->xE.center - xEast[0])*(t->xE.center - xEast[0]) +
	     (t->yE.center - yEast[0])*(t->yE.center - yEast[0]) <
	     (2.*sigmaEast[0])*(2.*sigmaEast[0]) ) {
	  if (t->ScintE.e1>0.) his[0][0]->Fill(t->ScintE.e1);
	  if (t->ScintE.e2>0.) his[0][1]->Fill(t->ScintE.e2);
	  if (t->ScintE.e3>0.) his[0][2]->Fill(t->ScintE.e3);
	  if (t->ScintE.e4>0.) his[0][3]->Fill(t->ScintE.e4);
	}
      }
      if (t->Side == 1) {
	if ( (t->xW.center - xWest[0])*(t->xW.center - xWest[0]) +
	     (t->yW.center - yWest[0])*(t->yW.center - yWest[0]) <
	     (2.*sigmaWest[0])*(2.*sigmaWest[0]) ) {
	  if (t->ScintW.e1>0.) his[0][4]->Fill(t->ScintW.e1);
	  if (t->ScintW.e2>0.) his[0][5]->Fill(t->ScintW.e2);
	  if (t->ScintW.e3>0.) his[0][6]->Fill(t->ScintW.e3);
	  if (t->ScintW.e4>0.) his[0][7]->Fill(t->ScintW.e4);
	}
      }
      
    }

    // Second source (x,y)
    if (nSources > 1 && useSource[1]) {

      if (t->Side == 0) {
	if ( (t->xE.center - xEast[1])*(t->xE.center - xEast[1]) +
	     (t->yE.center - yEast[1])*(t->yE.center - yEast[1]) <
	     (2.*sigmaEast[1])*(2.*sigmaEast[1]) ) {
	  if (t->ScintE.e1>0.) his[1][0]->Fill(t->ScintE.e1);
	  if (t->ScintE.e2>0.) his[1][1]->Fill(t->ScintE.e2);
	  if (t->ScintE.e3>0.) his[1][2]->Fill(t->ScintE.e3);
	  if (t->ScintE.e4>0.) his[1][3]->Fill(t->ScintE.e4);
	}
      }
      if (t->Side == 1) {
	if ( (t->xW.center - xWest[1])*(t->xW.center - xWest[1]) +
	     (t->yW.center - yWest[1])*(t->yW.center - yWest[1]) <
	     (2.*sigmaWest[1])*(2.*sigmaWest[1]) ) {
	  if (t->ScintW.e1>0.) his[1][4]->Fill(t->ScintW.e1);
	  if (t->ScintW.e2>0.) his[1][5]->Fill(t->ScintW.e2);
	  if (t->ScintW.e3>0.) his[1][6]->Fill(t->ScintW.e3);
	  if (t->ScintW.e4>0.) his[1][7]->Fill(t->ScintW.e4);
	}
      }
    }

    // Third source (x,y)
    if (nSources > 2 && useSource[2]) {

      if (t->Side == 0) {
	if ( (t->xE.center - xEast[2])*(t->xE.center - xEast[2]) +
	     (t->yE.center - yEast[2])*(t->yE.center - yEast[2]) <
	     (2.*sigmaEast[2])*(2.*sigmaEast[2]) ) {
	  if (t->ScintE.e1>0.) his[2][0]->Fill(t->ScintE.e1);
	  if (t->ScintE.e2>0.) his[2][1]->Fill(t->ScintE.e2);
	  if (t->ScintE.e3>0.) his[2][2]->Fill(t->ScintE.e3);
	  if (t->ScintE.e4>0.) his[2][3]->Fill(t->ScintE.e4);
	}
      }
      if (t->Side == 1) {
	if ( (t->xW.center - xWest[2])*(t->xW.center - xWest[2]) +
	     (t->yW.center - yWest[2])*(t->yW.center - yWest[2]) <
	     (2.*sigmaWest[2])*(2.*sigmaWest[2]) ) {
	  if (t->ScintW.e1>0.) his[2][4]->Fill(t->ScintW.e1);
	  if (t->ScintW.e2>0.) his[2][5]->Fill(t->ScintW.e2);
	  if (t->ScintW.e3>0.) his[2][6]->Fill(t->ScintW.e3);
	  if (t->ScintW.e4>0.) his[2][7]->Fill(t->ScintW.e4);
	}
      }

    }
    
  }

  delete t; //closes input file
  

  // Find maximum bin
  double maxBin[3][8];
  double maxCounts[3][8];
  for (int n=0; n<nSources; n++) {
    for (int j=0; j<8; j++) {
      maxCounts[n][j] = -1.0;
    }
  }

  double binCenter[3][8];
  double binCenterMax[3][8];
  double binCounts[3][8];
  for (int n=0; n<nSources; n++) {
    for (int i=0.045*nBin; i<nBin; i++) { //Don't want zeroth bin or any bins affected by forcing PMT curves to zero at low energy
      for (int j=0; j<8; j++) {
        binCenter[n][j] = his[n][j]->GetBinCenter(i+1);
        binCounts[n][j] = his[n][j]->GetBinContent(i+1);
        if (binCounts[n][j] > maxCounts[n][j]) {
          maxCounts[n][j] = binCounts[n][j];
          maxBin[n][j] = i;
          binCenterMax[n][j] = binCenter[n][j];
	}
      }
    }
  }
    
  
  // Define histogram fit ranges
  double xLow[4][8], xHigh[4][8];
  for (int n=0; n<nSources; n++) {

    for (int j=0; j<8; j++) {
      for (int i=maxBin[n][j]; i<nBin; i++) {
        if (his[n][j]->GetBinContent(i+1) < 0.33*maxCounts[n][j]) {
          xHigh[n][j] = his[n][j]->GetBinCenter(i+1);
	  if ((i-maxBin[n][j])<5) xHigh[n][j] = binCenterMax[n][j]+50.; //In case the statistics cause a poor determination of FWHM, set one to capture the fit
          break;
        }
	if (i==(nBin-1)) xHigh[n][j] = his[n][j]->GetBinCenter(i);
      }
      for (int i=maxBin[n][j]; i>0; i--) {
        if (his[n][j]->GetBinContent(i+1) < 0.5*maxCounts[n][j]) {
          xLow[n][j] = his[n][j]->GetBinCenter(i+1);
	  if ((maxBin[n][j]-i)<5) xLow[n][j] = binCenterMax[n][j]-30.; //In case the statistics cause a poor determination of FWHM, set one to capture the fit
          break;
        }
      }
    }

  }


  // Fit histograms
  TF1 *gaussian_fit[3][8];
  double fitMean[3][8];

  for (int n=0; n<nSources; n++) {
    if (useSource[n]) {
      for (int j=0; j<8; j++) {
	char fitName[500];
	sprintf(fitName, "gaussian_fit_%i_%i", n, j);

	gaussian_fit[n][j] = new TF1(fitName,
				     "gaus",
				     xLow[n][j], xHigh[n][j]);
	gaussian_fit[n][j]->SetParameter(0,maxCounts[n][j]);
	gaussian_fit[n][j]->SetParameter(1,binCenterMax[n][j]);
	gaussian_fit[n][j]->SetParameter(2,100.0);
	//gaussian_fit[n][j]->SetParLimits(1,xLow[n][j],xHigh[n][j]);

	fitMean[n][j] = 0.;
	double sigma = 0.;
	Int_t ntries = 0;
	while (fitMean[n][j]<xLow[n][j] || fitMean[n][j]>xHigh[n][j] || sigma>225. || sigma<5.) {
	  his[n][j]->Fit(fitName, "RQ");
	  fitMean[n][j] = gaussian_fit[n][j]->GetParameter(1);	 
	  sigma = gaussian_fit[n][j]->GetParameter(2);
	  //cout << fitMean[n][j] << endl;
	  gaussian_fit[n][j]->SetParameter(1,binCenterMax[n][j]-1.);
	  ntries++;
 
	  if (ntries>4) {
	    cout << "Run " << runNumber << " can't converge on " << sourceName[n] << " peak in PMT " << j << endl;
	    break;
	  }
	}
      }
    }
  }
    
  double lowBiFitMean[8] = {0.};
  if (useLowBiPeak) {
    char fitName[500];
    TF1 *lowBiGauss[8];
    for (int j=0;j<8;j++) {
      Int_t nbins = his[BiPeakIndex][j]->GetXaxis()->GetNbins();
      his[BiPeakIndex][j]->GetXaxis()->SetRangeUser(250.,650.);
      Int_t maxBin = his[BiPeakIndex][j]->GetMaximumBin();
      Int_t maxBinContent = his[BiPeakIndex][j]->GetBinContent(maxBin);
      his[BiPeakIndex][j]->GetXaxis()->SetRange(0,nbins);
      Double_t max = his[BiPeakIndex][j]->GetXaxis()->GetBinCenter(maxBin);
      //cout << maxBin << endl;
      Double_t Xmin, Xmax; //Calculate the max and min bins for fitting
      Xmin=250.;
      Xmax=650.;
      for (int i=maxBin; i<nbins; i++) {
	if (his[BiPeakIndex][j]->GetBinContent(i+1)<0.5*maxBinContent) {
	  Xmax = his[BiPeakIndex][j]->GetXaxis()->GetBinCenter(i+1); 
	  if (Xmax>600.) Xmax=560.;
	  break;}
	if (his[BiPeakIndex][j]->GetXaxis()->GetBinCenter(i)>600.) {
	  Xmax = 560.;
	  break;}
      }
      for (int i=maxBin; i>0; i--) {
	if (his[BiPeakIndex][j]->GetBinContent(i-1)<0.25*maxBinContent) {
	  Xmin = his[BiPeakIndex][j]->GetXaxis()->GetBinCenter(i-1); 
	  break;}
	if (his[BiPeakIndex][j]->GetXaxis()->GetBinCenter(i-1)<250.) {
	  Xmin = 250.;
	  break;}
      }
      //cout << Xmin << " " << Xmax << endl;

      sprintf(fitName, "lowBiGauss%i",j);
      lowBiGauss[j] = new TF1(fitName,"gaus",
				 Xmin, Xmax);
    
      lowBiGauss[j]->SetParameter(0,(float)maxBinContent);
      lowBiGauss[j]->SetParameter(1, max);
      lowBiGauss[j]->SetParameter(2, 100.0);
      //lowBiGauss[j]->SetParLimits(1, Xmin, Xmax);

      lowBiFitMean[j] = 0.;
      double sigma = 0.;
      Int_t ntries = 0;
      while (lowBiFitMean[j]<Xmin || lowBiFitMean[j]>Xmax || sigma>225. || sigma<5.) {
	
	his[BiPeakIndex][j]->Fit(fitName, "RQ+");
	lowBiFitMean[j] = lowBiGauss[j]->GetParameter(1);
	sigma = lowBiGauss[j]->GetParameter(2);
	lowBiGauss[j]->SetParameter(1,max-1.);
	ntries++;
	//cout << lowBiFitMean[j] << endl;
      //his[BiPeakIndex][j]->GetXaxis()->SetRange(0,nbins);
	if (ntries>4) {
	  cout << "Run " << runNumber << " can't converge on Bi2 peak in PMT " << j << endl;
	  break;
	}
      }
    }
  }

  // Write results to file
  char tempResults[500];
  sprintf(tempResults, "%s/source_peaks_EvisPMTbyPMT_%s.dat",getenv("SOURCE_PEAKS"), argv[1]);
  ofstream outResults(tempResults);

  for (int n=0; n<nSources; n++) {
    if (useSource[n]) {
      string srcName = sourceName[n];
      if (srcName=="Bi") srcName = srcName+"1";
      outResults << runNumber << " "
		 << srcName << " "
		 << fitMean[n][0] << " "
		 << fitMean[n][1] << " " 
		 << fitMean[n][2] << " "
		 << fitMean[n][3] << " "
		 << fitMean[n][4] << " "
		 << fitMean[n][5] << " "
		 << fitMean[n][6] << " "
		 << fitMean[n][7] << endl;
    }
  }
  if (useLowBiPeak) {
    outResults << runNumber << " "
	       << "Bi2" << " "
	       << lowBiFitMean[0] << " "
	       << lowBiFitMean[1] << " "
	       << lowBiFitMean[2] << " "
	       << lowBiFitMean[3] << " "
	       << lowBiFitMean[4] << " "
	       << lowBiFitMean[5] << " "
	       << lowBiFitMean[6] << " "
	       << lowBiFitMean[7] << endl;
  }
  outResults.close();

  // Write output ntuple
  fileOut->Write();
  fileOut->Close();

  return 0;
}
