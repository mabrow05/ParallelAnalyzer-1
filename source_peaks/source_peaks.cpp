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
#include "sourcePeaks.h"

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

  //Adding in here that we are only looking for runs with Ce, Sn, or Bi in them
  bool correctSource = false;
  bool useSource[3] = {false,false,false};
  for (int n=0; n<nSources; n++) {
    if (sourceName[n]=="Ce" || sourceName[n]=="Sn" || sourceName[n]=="Bi") {
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
  sprintf(tempOut, "%s/source_peaks_%s.root",getenv("SOURCE_PEAKS"), argv[1]);
  TFile *fileOut = new TFile(tempOut,"RECREATE");


  // Quenched energies from simulation
  //double EQ[3];
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
  his[0][0] = new TH1F("his1_E0", "", nBin,0.0,4000.0);
  his[0][1] = new TH1F("his1_E1", "", nBin,0.0,4000.0);
  his[0][2] = new TH1F("his1_E2", "", nBin,0.0,4000.0);
  his[0][3] = new TH1F("his1_E3", "", nBin,0.0,4000.0);
  his[0][4] = new TH1F("his1_W0", "", nBin,0.0,4000.0);
  his[0][5] = new TH1F("his1_W1", "", nBin,0.0,4000.0);
  his[0][6] = new TH1F("his1_W2", "", nBin,0.0,4000.0);
  his[0][7] = new TH1F("his1_W3", "", nBin,0.0,4000.0);

  his[1][0] = new TH1F("his2_E0", "", nBin,0.0,4000.0);
  his[1][1] = new TH1F("his2_E1", "", nBin,0.0,4000.0);
  his[1][2] = new TH1F("his2_E2", "", nBin,0.0,4000.0);
  his[1][3] = new TH1F("his2_E3", "", nBin,0.0,4000.0);
  his[1][4] = new TH1F("his2_W0", "", nBin,0.0,4000.0);
  his[1][5] = new TH1F("his2_W1", "", nBin,0.0,4000.0);
  his[1][6] = new TH1F("his2_W2", "", nBin,0.0,4000.0);
  his[1][7] = new TH1F("his2_W3", "", nBin,0.0,4000.0);

  his[2][0] = new TH1F("his3_E0", "", nBin,0.0,4000.0);
  his[2][1] = new TH1F("his3_E1", "", nBin,0.0,4000.0);
  his[2][2] = new TH1F("his3_E2", "", nBin,0.0,4000.0);
  his[2][3] = new TH1F("his3_E3", "", nBin,0.0,4000.0);
  his[2][4] = new TH1F("his3_W0", "", nBin,0.0,4000.0);
  his[2][5] = new TH1F("his3_W1", "", nBin,0.0,4000.0);
  his[2][6] = new TH1F("his3_W2", "", nBin,0.0,4000.0);
  his[2][7] = new TH1F("his3_W3", "", nBin,0.0,4000.0);

  // Open input ntuple
  char tempIn[500];
  sprintf(tempIn, "%s/replay_pass3_%s.root",getenv("REPLAY_PASS3"), argv[1]);
  TFile *fileIn = new TFile(tempIn, "READ");
  TTree *Tin = (TTree*)(fileIn->Get("pass3"));

  // Variables
  Tin->SetBranchAddress("pmt0_pass3", &pmt_pass3[0]);
  Tin->SetBranchAddress("pmt1_pass3", &pmt_pass3[1]);
  Tin->SetBranchAddress("pmt2_pass3", &pmt_pass3[2]);
  Tin->SetBranchAddress("pmt3_pass3", &pmt_pass3[3]);
  Tin->SetBranchAddress("pmt4_pass3", &pmt_pass3[4]);
  Tin->SetBranchAddress("pmt5_pass3", &pmt_pass3[5]);
  Tin->SetBranchAddress("pmt6_pass3", &pmt_pass3[6]);
  Tin->SetBranchAddress("pmt7_pass3", &pmt_pass3[7]);

  Tin->SetBranchAddress("xE_pass3", &xE_pass3);
  Tin->SetBranchAddress("yE_pass3", &yE_pass3);
  Tin->SetBranchAddress("xW_pass3", &xW_pass3);
  Tin->SetBranchAddress("yW_pass3", &yW_pass3);

  Tin->SetBranchAddress("PID_pass3",  &PID_pass3);
  Tin->SetBranchAddress("type_pass3", &type_pass3);
  Tin->SetBranchAddress("side_pass3", &side_pass3);
  Tin->SetBranchAddress("posError_pass3", &posError_pass3);

  int nEvents = Tin->GetEntries();
  cout << "Processing " << argv[1] << " ... " << endl;
  cout << "... nEvents = " << nEvents << endl;

  // Loop over events
  for (int i=0; i<nEvents; i++) {
    Tin->GetEvent(i);

    // Use Type 0 events
    if (type_pass3 != 0) continue;

    // First source (x,y)
    if (useSource[0]) {

      if (side_pass3 == 0) {
	if ( (xE_pass3 - xEast[0])*(xE_pass3 - xEast[0]) +
	     (yE_pass3 - yEast[0])*(yE_pass3 - yEast[0]) <
	     (2.*sigmaEast[0])*(2.*sigmaEast[0]) ) {
	  for (int p=0; p<4; p++) {
	    his[0][p]->Fill(pmt_pass3[p]);
	  }
	}
      }
      if (side_pass3 == 1) {
	if ( (xW_pass3 - xWest[0])*(xW_pass3 - xWest[0]) +
	     (yW_pass3 - yWest[0])*(yW_pass3 - yWest[0]) <
	     (2.*sigmaWest[0])*(2.*sigmaWest[0]) ) {
	  for (int p=4; p<8; p++) {
	    his[0][p]->Fill(pmt_pass3[p]);
	  }
	}
      }
    }

    // Second source (x,y)
    if (nSources > 1 && useSource[1]) {
      if (side_pass3 == 0) {
	if ( (xE_pass3 - xEast[1])*(xE_pass3 - xEast[1]) +
	     (yE_pass3 - yEast[1])*(yE_pass3 - yEast[1]) <
	     (2.*sigmaEast[1])*(2.*sigmaEast[1]) ) {
	  for (int p=0; p<4; p++) {
	    his[1][p]->Fill(pmt_pass3[p]);
	  }
	}
      }
      if (side_pass3 == 1) {
	if ( (xW_pass3 - xWest[1])*(xW_pass3 - xWest[1]) +
	     (yW_pass3 - yWest[1])*(yW_pass3 - yWest[1]) <
	     (2.*sigmaWest[1])*(2.*sigmaWest[1]) ) {
	  for (int p=4; p<8; p++) {
	    his[1][p]->Fill(pmt_pass3[p]);
	  }
	}
      }
    }

    // Third source (x,y)
    if (nSources > 2 && useSource[2]) {
      if (side_pass3 == 0) {
	if ( (xE_pass3 - xEast[2])*(xE_pass3 - xEast[2]) +
             (yE_pass3 - yEast[2])*(yE_pass3 - yEast[2]) <
             (2.*sigmaEast[2])*(2.*sigmaEast[2]) ) {
          for (int p=0; p<4; p++) {
            his[2][p]->Fill(pmt_pass3[p]);
          }
        }
      }
      if (side_pass3 == 1) {
	if ( (xW_pass3 - xWest[2])*(xW_pass3 - xWest[2]) +
             (yW_pass3 - yWest[2])*(yW_pass3 - yWest[2]) <
             (2.*sigmaWest[2])*(2.*sigmaWest[2]) ) {
          for (int p=4; p<8; p++) {
            his[2][p]->Fill(pmt_pass3[p]);
          }
	}
      }
    }

  }

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
    for (int i=0; i<nBin; i++) {
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
  double xLow[3][8], xHigh[3][8];
  for (int n=0; n<nSources; n++) {

    for (int j=0; j<8; j++) {
      for (int i=maxBin[n][j]; i<nBin; i++) {
        if (his[n][j]->GetBinContent(i+1) < 0.5*maxCounts[n][j]) {
          xHigh[n][j] = his[n][j]->GetBinCenter(i+1);
          break;
        }
      }
      for (int i=maxBin[n][j]; i>0; i--) {
        if (his[n][j]->GetBinContent(i+1) < 0.5*maxCounts[n][j]) {
          xLow[n][j] = his[n][j]->GetBinCenter(i+1);
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
				     "[0]*exp(-(x-[1])*(x-[1])/(2.*[2]*[2]))",
				     xLow[n][j], xHigh[n][j]);
	gaussian_fit[n][j]->SetParameter(0,maxCounts[n][j]);
	gaussian_fit[n][j]->SetParameter(1,binCenterMax[n][j]);
	gaussian_fit[n][j]->SetParameter(2,100.0);
	gaussian_fit[n][j]->SetParLimits(1,xLow[n][j],xHigh[n][j]);

	his[n][j]->Fit(fitName, "LRQ");
	fitMean[n][j] = gaussian_fit[n][j]->GetParameter(1);
      }
    }
  }

  // Write results to file
  char tempResults[500];
  sprintf(tempResults, "%s/source_peaks_%s.dat",getenv("SOURCE_PEAKS"), argv[1]);
  ofstream outResults(tempResults);

  for (int n=0; n<nSources; n++) {
    if (useSource[n]) {
      if (sourceName[n]=="Bi") sourceName[n]=sourceName[n]+"1";
      outResults << runNumber << " "
		 << sourceName[n] << " "
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

  // Write output ntuple
  fileOut->Write();
  fileOut->Close();

  return 0;
}
