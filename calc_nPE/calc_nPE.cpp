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
#include "posMapReader.h"
#include "runInfo.h"
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

  //Adding in here that we are only using runs with Sn
  bool correctSource = false;
  bool useSource[3] = {false,false,false};
  for (int n=0; n<nSources; n++) {
    if (sourceName[n]=="Sn") {
      correctSource = true;
      useSource[n]=true;
      continue; 
    }
  }
  
  if (!correctSource) {
    cout << "Source run with no Sn\n";
    exit(0);
  }

  // Open output ntuple
  char tempOut[500];
  sprintf(tempOut, "%s/nPE_weights_%s.root",getenv("NPE_WEIGHTS"), argv[1]);
  TFile *fileOut = new TFile(tempOut,"RECREATE");

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

  TH1F *his[8];
  his[0] = new TH1F("his1_E0", "", nBin,0.0,4000.0);
  his[1] = new TH1F("his1_E1", "", nBin,0.0,4000.0);
  his[2] = new TH1F("his1_E2", "", nBin,0.0,4000.0);
  his[3] = new TH1F("his1_E3", "", nBin,0.0,4000.0);
  his[4] = new TH1F("his1_W0", "", nBin,0.0,4000.0);
  his[5] = new TH1F("his1_W1", "", nBin,0.0,4000.0);
  his[6] = new TH1F("his1_W2", "", nBin,0.0,4000.0);
  his[7] = new TH1F("his1_W3", "", nBin,0.0,4000.0);

  /*his[1][0] = new TH1F("his2_E0", "", nBin,0.0,4000.0);
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
  his[2][7] = new TH1F("his3_W3", "", nBin,0.0,4000.0);*/

  // We are also going to calculate the average value of the position correction from the 
  // position map for every event. This will be needed later on when calculating an adjusted 
  // number of photoelectrons per event in simulation
  GetPositionMap(getXeRunPeriod(runNumber)); //Reads in the proper position map
  Int_t numDataPoints[2] = {0}; //Holds the number of data points for each side
  Double_t aveEta[8] = {0.}; //Holds the average value of eta for the the data being read in
  Int_t EastXgrid, EastYgrid, WestXgrid, WestYgrid; // holds the grid point of the data being read in
  EastXgrid=EastYgrid=WestXgrid=WestYgrid=0;


  // Open input ntuple
  char tempIn[500];
  sprintf(tempIn, "%s/replay_pass3_%s.root",getenv("REPLAY_PASS3"), argv[1]);
  DataTree *t = new DataTree();
  t->setupInputTree(std::string(tempIn),"pass3");

  int nEvents = t->getEntries();
  cout << "Processing " << argv[1] << " ... " << endl;
  cout << "... nEvents = " << nEvents << endl;

  // Loop over events
  for (int i=0; i<nEvents; i++) {
    t->getEvent(i);

    // Use Type 0 events
    if (t->Type != 0) continue;
    

    // First source (x,y)
    if (useSource[0]) {
      vector < vector <int> > gridPoint = getGridPoint(t->xE.center, t->yE.center, t->xW.center, t->yW.center);
      EastXgrid = gridPoint[0][0];
      EastYgrid = gridPoint[0][1];
      WestXgrid = gridPoint[1][0];
      WestYgrid = gridPoint[1][1];
      //cout << EastXgrid << endl;
      if (t->Side == 0) {
	if ( (t->xE.center - xEast[0])*(t->xE.center - xEast[0]) +
	     (t->yE.center - yEast[0])*(t->yE.center - yEast[0]) <
	     (2.*sigmaEast[0])*(2.*sigmaEast[0]) ) {

	  his[0]->Fill(t->ScintE.q1);
	  his[1]->Fill(t->ScintE.q2);
	  his[2]->Fill(t->ScintE.q3);
	  his[3]->Fill(t->ScintE.q4);

	  numDataPoints[0]++;
	  aveEta[0]+=positionMap[0][EastXgrid][EastYgrid];
	  aveEta[1]+=positionMap[1][EastXgrid][EastYgrid];
	  aveEta[2]+=positionMap[2][EastXgrid][EastYgrid];
	  aveEta[3]+=positionMap[3][EastXgrid][EastYgrid];
	  
	}
      }
      if (t->Side == 1) {
	if ( (t->xW.center - xWest[0])*(t->xW.center - xWest[0]) +
	     (t->yW.center - yWest[0])*(t->yW.center - yWest[0]) <
	     (2.*sigmaWest[0])*(2.*sigmaWest[0]) ) {
	 
	  his[4]->Fill(t->ScintW.q1);
	  his[5]->Fill(t->ScintW.q2);
	  his[6]->Fill(t->ScintW.q3);
	  his[7]->Fill(t->ScintW.q4);

	  numDataPoints[1]++;
	  aveEta[4]+=positionMap[4][WestXgrid][WestYgrid];
	  aveEta[5]+=positionMap[5][WestXgrid][WestYgrid];
	  aveEta[6]+=positionMap[6][WestXgrid][WestYgrid];
	  aveEta[7]+=positionMap[7][WestXgrid][WestYgrid];
	}
      }
    }

    // Second source (x,y)
    if (nSources > 1 && useSource[1]) {
      vector < vector <int> > gridPoint = getGridPoint(t->xE.center, t->yE.center, t->xW.center, t->yW.center);
      EastXgrid = gridPoint[0][0];
      EastYgrid = gridPoint[0][1];
      WestXgrid = gridPoint[1][0];
      WestYgrid = gridPoint[1][1];
      //cout << EastXgrid << endl;
      if (t->Side == 0) {
	if ( (t->xE.center - xEast[1])*(t->xE.center - xEast[1]) +
	     (t->yE.center - yEast[1])*(t->yE.center - yEast[1]) <
	     (2.*sigmaEast[1])*(2.*sigmaEast[1]) ) {
	  
	  his[0]->Fill(t->ScintE.q1);
	  his[1]->Fill(t->ScintE.q2);
	  his[2]->Fill(t->ScintE.q3);
	  his[3]->Fill(t->ScintE.q4);
	
	  numDataPoints[0]++;
	  aveEta[0]+=positionMap[0][EastXgrid][EastYgrid];
	  aveEta[1]+=positionMap[1][EastXgrid][EastYgrid];
	  aveEta[2]+=positionMap[2][EastXgrid][EastYgrid];
	  aveEta[3]+=positionMap[3][EastXgrid][EastYgrid];
	}
      }
      if (t->Side == 1) {
	if ( (t->xW.center - xWest[1])*(t->xW.center - xWest[1]) +
	     (t->yW.center - yWest[1])*(t->yW.center - yWest[1]) <
	     (2.*sigmaWest[1])*(2.*sigmaWest[1]) ) {
	  
	  his[4]->Fill(t->ScintW.q1);
	  his[5]->Fill(t->ScintW.q2);
	  his[6]->Fill(t->ScintW.q3);
	  his[7]->Fill(t->ScintW.q4);
	  
	  numDataPoints[1]++;
	  aveEta[4]+=positionMap[4][WestXgrid][WestYgrid];
	  aveEta[5]+=positionMap[5][WestXgrid][WestYgrid];
	  aveEta[6]+=positionMap[6][WestXgrid][WestYgrid];
	  aveEta[7]+=positionMap[7][WestXgrid][WestYgrid];
	}
      }
    }

    // Third source (x,y)
    if (nSources > 2 && useSource[2]) {
      vector < vector <int> > gridPoint = getGridPoint(t->xE.center, t->yE.center, t->xW.center, t->yW.center);
      EastXgrid = gridPoint[0][0];
      EastYgrid = gridPoint[0][1];
      WestXgrid = gridPoint[1][0];
      WestYgrid = gridPoint[1][1];
      //cout << EastXgrid << endl;
      if (t->Side == 0) {
	if ( (t->xE.center - xEast[2])*(t->xE.center - xEast[2]) +
             (t->yE.center - yEast[2])*(t->yE.center - yEast[2]) <
             (2.*sigmaEast[2])*(2.*sigmaEast[2]) ) {
         
	  his[0]->Fill(t->ScintE.q1);
	  his[1]->Fill(t->ScintE.q2);
	  his[2]->Fill(t->ScintE.q3);
	  his[3]->Fill(t->ScintE.q4);

	  numDataPoints[0]++;
	  aveEta[0]+=positionMap[0][EastXgrid][EastYgrid];
	  aveEta[1]+=positionMap[1][EastXgrid][EastYgrid];
	  aveEta[2]+=positionMap[2][EastXgrid][EastYgrid];
	  aveEta[3]+=positionMap[3][EastXgrid][EastYgrid];
        }
      }
      if (t->Side == 1) {
	if ( (t->xW.center - xWest[2])*(t->xW.center - xWest[2]) +
             (t->yW.center - yWest[2])*(t->yW.center - yWest[2]) <
             (2.*sigmaWest[2])*(2.*sigmaWest[2]) ) {
          
	  his[4]->Fill(t->ScintW.q1);
	  his[5]->Fill(t->ScintW.q2);
	  his[6]->Fill(t->ScintW.q3);
	  his[7]->Fill(t->ScintW.q4);
	  
	  numDataPoints[1]++;
	  aveEta[4]+=positionMap[4][WestXgrid][WestYgrid];
	  aveEta[5]+=positionMap[5][WestXgrid][WestYgrid];
	  aveEta[6]+=positionMap[6][WestXgrid][WestYgrid];
	  aveEta[7]+=positionMap[7][WestXgrid][WestYgrid];
	}
      }
    }

  }

  // Find maximum bin
  double maxBin[8];
  double maxCounts[8];
  for (int j=0; j<8; j++) {
    maxCounts[j] = -1.0; 
  }

  double binCenter[8];
  double binCenterMax[8];
  double binCounts[8];
  for (int i=0; i<nBin; i++) {
    for (int j=0; j<8; j++) {
      binCenter[j] = his[j]->GetBinCenter(i+1);
      binCounts[j] = his[j]->GetBinContent(i+1);
      if (binCounts[j] > maxCounts[j]) {
	maxCounts[j] = binCounts[j];
	maxBin[j] = i;
	binCenterMax[j] = binCenter[j];
      }
    }
  }
  

  // Define histogram fit ranges
  double xLow[8], xHigh[8];

  for (int j=0; j<8; j++) {
    for (int i=maxBin[j]; i<nBin; i++) {
      if (his[j]->GetBinContent(i+1) < 0.5*maxCounts[j]) {
	xHigh[j] = his[j]->GetBinCenter(i+1);
	break;
      }
    }
    for (int i=maxBin[j]; i>0; i--) {
      if (his[j]->GetBinContent(i+1) < 0.5*maxCounts[j]) {
	xLow[j] = his[j]->GetBinCenter(i+1);
	break;
      }
    }
  }



  // Fit histograms
  TF1 *gaussian_fit[8];
  double fitMean[8], fitSigma[8], nPE[8], nPE_per_channel[8];

  for (int j=0; j<8; j++) {
    char fitName[500];
    sprintf(fitName, "gaussian_fit_%i", j);
    
    gaussian_fit[j] = new TF1(fitName,
				 "gaus",
				 xLow[j], xHigh[j]);
    gaussian_fit[j]->SetParameter(0,maxCounts[j]);
    gaussian_fit[j]->SetParameter(1,binCenterMax[j]);
    gaussian_fit[j]->SetParameter(2,100.0);
    gaussian_fit[j]->SetParLimits(1,xLow[j],xHigh[j]);
    
    his[j]->Fit(fitName, "LRQ");
    fitMean[j] = gaussian_fit[j]->GetParameter(1);
    fitSigma[j] = gaussian_fit[j]->GetParameter(2);
    nPE[j] = (fitMean[j]*fitMean[j])/(fitSigma[j]*fitSigma[j]);
    nPE_per_channel[j] = nPE[j]/fitMean[j];
  }

  // Write results to file
  char tempResults[500];
  sprintf(tempResults, "%s/nPE_weights_%s.dat",getenv("NPE_WEIGHTS"), argv[1]);
  ofstream outResults(tempResults);

  for (int j=0; j<8; j++) {
    outResults << fitMean[j] << " " 
	       << fitSigma[j] << " " 
	       << nPE[j] << " "
	       << nPE_per_channel[j] << endl;
  }
  outResults.close();

  //Writing average eta values to file
  sprintf(tempResults, "%s/nPE_meanEtaVal_%s.dat",getenv("NPE_WEIGHTS"), argv[1]);
  outResults.open(tempResults);
  for (int j=0; j<4; j++) {
    outResults << aveEta[j]/(double)numDataPoints[0] << endl;
  }
  for (int j=4; j<8; j++) {
    outResults << aveEta[j]/(double)numDataPoints[1] << endl;
  }
  outResults.close();

  // Write output ntuple
  fileOut->Write();
  fileOut->Close();

  return 0;
}
