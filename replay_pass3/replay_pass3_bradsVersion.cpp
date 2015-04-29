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

using namespace std;

int main(int argc, char *argv[])
{
  cout.setf(ios::fixed, ios::floatfield);
  cout.precision(12);

  // Position bins
  int nPMT = 8;
  int nPosBinsX = 11;
  int nPosBinsY = 11;  
  double xBinWidth = 10.0;
  double yBinWidth = 10.0;
  double xBinLower[nPosBinsX];
  double xBinUpper[nPosBinsX];
  double xBinCenter[nPosBinsX];
  double yBinLower[nPosBinsY];
  double yBinUpper[nPosBinsY];
  double yBinCenter[nPosBinsY];
  int intXBinCenter[nPosBinsX];
  int intYBinCenter[nPosBinsY];

  for (int k=0; k<nPosBinsX; k++) {
    xBinLower[k]     = -55.0 + ((double) k)*xBinWidth;
    xBinUpper[k]     = -55.0 + ((double) k)*xBinWidth + xBinWidth;
    xBinCenter[k]    = (xBinLower[k] + xBinUpper[k])/2.;
    intXBinCenter[k] = (int) xBinCenter[k];
    cout << xBinLower[k] << " " << intXBinCenter[k] << " " << xBinUpper[k] << endl;
  }

  for (int k=0; k<nPosBinsY; k++) {
    yBinLower[k]     = -55.0 + ((double) k)*yBinWidth;
    yBinUpper[k]     = -55.0 + ((double) k)*yBinWidth + yBinWidth;
    yBinCenter[k]    = (yBinLower[k] + yBinUpper[k])/2.;
    intYBinCenter[k] = (int) yBinCenter[k];
    cout << yBinLower[k] << " " << intYBinCenter[k] << " " << yBinUpper[k] << endl;
  }

  // Run number integer
  cout << "Run " << argv[1] << " ..." << endl;
  cout << "... Applying Xe position map ..." << endl;
  istringstream ss(argv[1]);
  int runNumber;
  ss >> runNumber;

  // Determine position map to use
  char tempFileXePositionMap[500];
  if (runNumber < 18081) {
    sprintf(tempFileXePositionMap, "../position_map/position_map_2.dat");
  }
  else if (runNumber < 18390) {
    sprintf(tempFileXePositionMap, "../position_map/position_map_3.dat");
  }
  else if (runNumber < 18712) {
    sprintf(tempFileXePositionMap, "../position_map/position_map_4.dat");
  }
  else if (runNumber < 19589) {
    sprintf(tempFileXePositionMap, "../position_map/position_map_5.dat");
  }
  else if (runNumber < 19873) {
    sprintf(tempFileXePositionMap, "../position_map/position_map_6.dat");
  }
  else if (runNumber < 19967) {
    sprintf(tempFileXePositionMap, "../position_map/position_map_7.dat");
  }
  cout << "... Reading: " << tempFileXePositionMap << endl;

  // Read position map
  double x, y;
  double positionMap[nPMT][nPosBinsX][nPosBinsY];
  ifstream fileXePositionMap(tempFileXePositionMap);
  for (int i=0; i<nPosBinsX; i++) {
    for (int j=0; j<nPosBinsY; j++) {
      fileXePositionMap >> x >> y
                        >> positionMap[0][i][j]
			>> positionMap[1][i][j]
			>> positionMap[2][i][j]
			>> positionMap[3][i][j]
			>> positionMap[4][i][j]
			>> positionMap[5][i][j]
			>> positionMap[6][i][j]
			>> positionMap[7][i][j];
    }
  }

  // Calculate (x,y) correction factor: eta = 1/positionMap
  double eta[nPMT][nPosBinsX][nPosBinsY];
  for (int i=0; i<nPosBinsX; i++) {
    for (int j=0; j<nPosBinsY; j++) {
      eta[0][i][j] = 1. / positionMap[0][i][j];
      eta[1][i][j] = 1. / positionMap[1][i][j];
      eta[2][i][j] = 1. / positionMap[2][i][j];
      eta[3][i][j] = 1. / positionMap[3][i][j];
      eta[4][i][j] = 1. / positionMap[4][i][j];
      eta[5][i][j] = 1. / positionMap[5][i][j];
      eta[6][i][j] = 1. / positionMap[6][i][j];
      eta[7][i][j] = 1. / positionMap[7][i][j];
    }
  }

  // Open output ntuple
  char tempOut[500];
  sprintf(tempOut, "/extern/UCNA/replay_pass3/replay_pass3_%s.root", argv[1]);
  TFile *fileOut = new TFile(tempOut,"RECREATE");
  TTree *Tout = new TTree("pass3", "pass3");

  // Variables
  Tout->Branch("pmt0_pass3", &pmt_pass3[0], "pmt0_pass3/D");
  Tout->Branch("pmt1_pass3", &pmt_pass3[1], "pmt1_pass3/D");
  Tout->Branch("pmt2_pass3", &pmt_pass3[2], "pmt2_pass3/D");
  Tout->Branch("pmt3_pass3", &pmt_pass3[3], "pmt3_pass3/D");
  Tout->Branch("pmt4_pass3", &pmt_pass3[4], "pmt4_pass3/D");
  Tout->Branch("pmt5_pass3", &pmt_pass3[5], "pmt5_pass3/D");
  Tout->Branch("pmt6_pass3", &pmt_pass3[6], "pmt6_pass3/D");
  Tout->Branch("pmt7_pass3", &pmt_pass3[7], "pmt7_pass3/D");

  Tout->Branch("xE_pass3", &xE_pass3, "xE_pass3/D");
  Tout->Branch("yE_pass3", &yE_pass3, "yE_pass3/D");
  Tout->Branch("xW_pass3", &xW_pass3, "xW_pass3/D");
  Tout->Branch("yW_pass3", &yW_pass3, "yW_pass3/D");

  Tout->Branch("PID_pass3",  &PID_pass3,  "PID_pass3/I");
  Tout->Branch("type_pass3", &type_pass3, "type_pass3/I");
  Tout->Branch("side_pass3", &side_pass3, "side_pass3/I");
  Tout->Branch("posError_pass3", &posError_pass3, "posError_pass3/I");

  // Open input ntuple
  char tempIn[500];
  sprintf(tempIn, "/extern/UCNA/replay_pass2/replay_pass2_%s.root", argv[1]);

  TFile *fileIn = new TFile(tempIn, "READ");
  TTree *Tin = (TTree*)(fileIn->Get("pass2"));

  // Variables
  Tin->SetBranchAddress("pmt0_pass2", &pmt_pass2[0]);
  Tin->SetBranchAddress("pmt1_pass2", &pmt_pass2[1]);
  Tin->SetBranchAddress("pmt2_pass2", &pmt_pass2[2]);
  Tin->SetBranchAddress("pmt3_pass2", &pmt_pass2[3]);
  Tin->SetBranchAddress("pmt4_pass2", &pmt_pass2[4]);
  Tin->SetBranchAddress("pmt5_pass2", &pmt_pass2[5]);
  Tin->SetBranchAddress("pmt6_pass2", &pmt_pass2[6]);
  Tin->SetBranchAddress("pmt7_pass2", &pmt_pass2[7]);

  Tin->SetBranchAddress("xE_pass2", &xE_pass2);
  Tin->SetBranchAddress("yE_pass2", &yE_pass2);
  Tin->SetBranchAddress("xW_pass2", &xW_pass2);
  Tin->SetBranchAddress("yW_pass2", &yW_pass2);

  Tin->SetBranchAddress("PID_pass2",  &PID_pass2);
  Tin->SetBranchAddress("type_pass2", &type_pass2);
  Tin->SetBranchAddress("side_pass2", &side_pass2);
  Tin->SetBranchAddress("posError_pass2", &posError_pass2);

  int nEvents = Tin->GetEntries();
  cout << "... Processing nEvents = " << nEvents << endl;

  // Loop over events
  for (int i=0; i<nEvents; i++) {
    Tin->GetEvent(i);

    // Determine (x,y) bin center
    int intEastBinX = -1;
    int intEastBinY = -1;
    int intWestBinX = -1;
    int intWestBinY = -1;

    for (int m=0; m<nPosBinsX; m++) {
      if ( (xE_pass2 >= xBinLower[m]) && (xE_pass2 < xBinUpper[m]) ) intEastBinX = m;
      if ( (xW_pass2 >= xBinLower[m]) && (xW_pass2 < xBinUpper[m]) ) intWestBinX = m;
    }

    for (int m=0; m<nPosBinsY; m++) {
      if ( (yE_pass2 >= yBinLower[m]) && (yE_pass2 < yBinUpper[m]) ) intEastBinY = m;
      if ( (yW_pass2 >= yBinLower[m]) && (yW_pass2 < yBinUpper[m]) ) intWestBinY = m;
    }

    // Determine (x,y) bin limits
    int intEastBinXCenterLow  = -1;
    int intEastBinXCenterHigh = -1;
    if (xE_pass2 < xBinCenter[0]) {
      intEastBinXCenterLow  = -1;
      intEastBinXCenterHigh = 0;
    }
    for (int m=0; m<nPosBinsX-1; m++) {
      if ( (xE_pass2 >= xBinCenter[m]) && (xE_pass2 < xBinCenter[m+1]) ) {
        intEastBinXCenterLow  = m;
        intEastBinXCenterHigh = m+1;
      }
    }
    if (xE_pass2 > xBinCenter[nPosBinsX-1]) {
      intEastBinXCenterLow  = nPosBinsX - 1;
      intEastBinXCenterHigh = nPosBinsX;
    }

    int intEastBinYCenterLow  = -1;
    int intEastBinYCenterHigh = -1;
    if (yE_pass2 < yBinCenter[0]) {
      intEastBinYCenterLow  = -1;
      intEastBinYCenterHigh = 0;
    }
    for (int m=0; m<nPosBinsY-1; m++) {
      if ( (yE_pass2 >= yBinCenter[m]) && (yE_pass2 < yBinCenter[m+1]) ) {
        intEastBinYCenterLow  = m;
        intEastBinYCenterHigh = m+1;
      }
    }
    if (yE_pass2 > yBinCenter[nPosBinsY-1]) {
      intEastBinYCenterLow  = nPosBinsY - 1;
      intEastBinYCenterHigh = nPosBinsY;
    }

    int intWestBinXCenterLow  = -1;
    int intWestBinXCenterHigh = -1;
    if (xW_pass2 < xBinCenter[0]) {
      intWestBinXCenterLow  = -1;
      intWestBinXCenterHigh = 0;
    }
    for (int m=0; m<nPosBinsX-1; m++) {
      if ( (xW_pass2 >= xBinCenter[m]) && (xW_pass2 < xBinCenter[m+1]) ) {
        intWestBinXCenterLow  = m;
        intWestBinXCenterHigh = m+1;
      }
    }
    if (xW_pass2 > xBinCenter[nPosBinsX-1]) {
      intWestBinXCenterLow  = nPosBinsX - 1;
      intWestBinXCenterHigh = nPosBinsX;
    }

    int intWestBinYCenterLow  = -1;
    int intWestBinYCenterHigh = -1;
    if (yW_pass2 < yBinCenter[0]) {
      intWestBinYCenterLow  = -1;
      intWestBinYCenterHigh = 0;
    }
    for (int m=0; m<nPosBinsY-1; m++) {
      if ( (yW_pass2 >= yBinCenter[m]) && (yW_pass2 < yBinCenter[m+1]) ) {
        intWestBinYCenterLow  = m;
        intWestBinYCenterHigh = m+1;
      }
    }
    if (yW_pass2 > yBinCenter[nPosBinsY-1]) {
      intWestBinYCenterLow  = nPosBinsY - 1;
      intWestBinYCenterHigh = nPosBinsY;
    }

    // Determine (x,y) correction factor
    double x1, x2, y1, y2;
    double eta11[8], eta21[8], eta12[8], eta22[8];
    double etaE[8], etaW[8];

    // East
    if ( (intEastBinXCenterLow > -1 && intEastBinXCenterHigh < nPosBinsX+1) &&
         (intEastBinYCenterLow > -1 && intEastBinYCenterHigh < nPosBinsX+1) ) {
      x1 = xBinCenter[intEastBinXCenterLow];
      x2 = xBinCenter[intEastBinXCenterHigh];
      y1 = yBinCenter[intEastBinYCenterLow];
      y2 = yBinCenter[intEastBinYCenterHigh];
      for (int p=0; p<4; p++) {
        eta11[p] = eta[p][intEastBinXCenterLow][intEastBinYCenterLow];
        eta21[p] = eta[p][intEastBinXCenterHigh][intEastBinYCenterLow];
        eta12[p] = eta[p][intEastBinXCenterLow][intEastBinYCenterHigh];
        eta22[p] = eta[p][intEastBinXCenterHigh][intEastBinYCenterHigh];

        etaE[p] = (1./((x2-x1)*(y2-y1))) * (eta11[p]*(x2-xE_pass2)*(y2-yE_pass2) +
                                            eta21[p]*(xE_pass2-x1)*(y2-yE_pass2) +
                                            eta12[p]*(x2-xE_pass2)*(yE_pass2-y1) +
                                            eta22[p]*(xE_pass2-x1)*(yE_pass2-y1));

        pmt_pass3[p] = pmt_pass2[p] * etaE[p];
      }
    }
    else {
      for (int p=0; p<4; p++) {
        pmt_pass3[p] = pmt_pass2[p] * eta[p][intEastBinX][intEastBinY];
      }
    }

    // West
    if ( (intWestBinXCenterLow > -1 && intWestBinXCenterHigh < nPosBinsX+1) &&
         (intWestBinYCenterLow > -1 && intWestBinYCenterHigh < nPosBinsX+1) ) {
      x1 = xBinCenter[intWestBinXCenterLow];
      x2 = xBinCenter[intWestBinXCenterHigh];
      y1 = yBinCenter[intWestBinYCenterLow];
      y2 = yBinCenter[intWestBinYCenterHigh];
      for (int p=4; p<8; p++) {
        eta11[p] = eta[p][intWestBinXCenterLow][intWestBinYCenterLow];
        eta21[p] = eta[p][intWestBinXCenterHigh][intWestBinYCenterLow];
        eta12[p] = eta[p][intWestBinXCenterLow][intWestBinYCenterHigh];
        eta22[p] = eta[p][intWestBinXCenterHigh][intWestBinYCenterHigh];

        etaW[p] = (1./((x2-x1)*(y2-y1))) * (eta11[p]*(x2-xW_pass2)*(y2-yW_pass2) +
                                            eta21[p]*(xW_pass2-x1)*(y2-yW_pass2) +
                                            eta12[p]*(x2-xW_pass2)*(yW_pass2-y1) +
                                            eta22[p]*(xW_pass2-x1)*(yW_pass2-y1));

        pmt_pass3[p] = pmt_pass2[p] * etaW[p];
      }
    }
    else {
      for (int p=4; p<8; p++) {
        pmt_pass3[p] = pmt_pass2[p] * eta[p][intWestBinX][intWestBinY];
      }
    }

    // Pass other variables from pass2 to pass3
    xE_pass3 = xE_pass2;
    yE_pass3 = yE_pass2;
    xW_pass3 = xW_pass2;
    yW_pass3 = yW_pass2;

    PID_pass3 = PID_pass2;
    type_pass3 = type_pass2;
    side_pass3 = side_pass2;
    posError_pass3 = posError_pass2;

    Tout->Fill();
  }

  // Write output ntuple
  fileOut->Write();
  fileOut->Close();

  return 0;
}
