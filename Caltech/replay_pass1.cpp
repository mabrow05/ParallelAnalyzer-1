//////////////////////////////////////////////////////////////////
// Code to apply cuts and do pedestal subtraction. 
// Ultimately produces first trees with branches of physical 
// meaning to be further processed in later replays.
//////////////////////////////////////////////////////////////////

#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>

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

//swank's included file, 1 of 1. 
#include "WireChamberResponse.h"

using namespace std;

int main(int argc, char *argv[])
{
  cout.setf(ios::fixed, ios::floatfield);
  cout.precision(12);

  cout << "Run " << argv[1] << " ..." << endl;

  // Read cuts file
  char tempFileCuts[500];
  sprintf(tempFileCuts, "%s/cuts_MB/cuts_%s.dat", getenv("PARALLEL_DATA_PATH"),argv[1]);
  cout << "... Reading: " << tempFileCuts << endl;

  ifstream fileCuts(tempFileCuts);
  fileCuts >> cutBeamBurstTime >> comment;
  fileCuts >> nCutsTimeWindows >> comment;
  if (nCutsTimeWindows > 0) {
    for (int i=0; i<nCutsTimeWindows; i++) {
      fileCuts >> cutTimeWindowLower[i] >> cutTimeWindowUpper[i];
    }
  }
  fileCuts >> cutEastAnode >> comment;
  fileCuts >> cutWestAnode >> comment;
  fileCuts >> cutEastTwoFold >> comment;
  fileCuts >> cutWestTwoFold >> comment;
  fileCuts >> cutEastTopVetoQADC >> comment;
  fileCuts >> cutEastTopVetoTDC >> comment;
  fileCuts >> cutEastDriftTubeTAC >> comment;
  fileCuts >> cutWestDriftTubeTAC >> comment;
  fileCuts >> cutEastBackingVetoQADC >> comment;
  fileCuts >> cutEastBackingVetoTDC >> comment;
  fileCuts >> cutWestBackingVetoQADC >> comment;
  fileCuts >> cutWestBackingVetoTDC >> comment;

  cout << "... Beam Burst T0 Cut: " << cutBeamBurstTime << endl;
  cout << "... Number of Time Windows Cuts: " << nCutsTimeWindows << endl;
  if (nCutsTimeWindows > 0) {
    for (int i=0; i<nCutsTimeWindows; i++) {
      cout << "        [" << cutTimeWindowLower[i] << ", " << cutTimeWindowUpper[i] << "]" << endl;
    }
  }
  cout << "... East MWPC Anode Cut: " << cutEastAnode << endl;
  cout << "... West MWPC Anode Cut: " << cutWestAnode << endl;
  cout << "... East Scintillator Two-Fold Trigger Cut: " << cutEastTwoFold << endl;
  cout << "... West Scintillator Two-Fold Trigger Cut: " << cutWestTwoFold << endl;
  cout << "... East Top Veto QADC Cut: " << cutEastTopVetoQADC << endl;
  cout << "... East Top Veto TDC Cut: " << cutEastTopVetoTDC << endl;
  cout << "... East Drift Tube TAC Cut: " << cutEastDriftTubeTAC << endl;
  cout << "... West Drift Tube TAC Cut: " << cutWestDriftTubeTAC << endl;
  cout << "... East Backing Veto QADC Cut: " << cutEastBackingVetoQADC << endl;
  cout << "... East Backing Veto TDC Cut: " << cutEastBackingVetoTDC << endl;
  cout << "... West Backing Veto QADC Cut: " << cutWestBackingVetoQADC << endl;
  cout << "... West Backing Veto TDC Cut: " << cutWestBackingVetoTDC << endl;

  // Read pedestals file
  char tempFilePed[500];
  int iRun;
  sprintf(tempFilePed, "/extern/UCNA/pedestals_MB/pedestals_%s.dat", argv[1]);
  cout << "... Reading: " << tempFilePed << endl;

  ifstream filePed(tempFilePed);
  for (int i=0; i<8; i++) {
    filePed >> iRun >> pedQadc[i];
  }
  for (int i=0; i<32; i++) {
    filePed >> iRun >> pedPdc2[i];
  }
  for (int i=0; i<32; i++) {
    filePed >> iRun >> pedPadc[i];
  }
  //filePed >> iRun >> pedPdc30;
  //filePed >> iRun >> pedPdc34;

  //cout << iRun << " " << pedPdc30 << endl;
  //cout << iRun << " " << pedPdc34 << endl;
  
  // Open output ntuple
  char tempOut[500];
  sprintf(tempOut, "%s/replay_pass1_%s.root",getenv("REPLAY_PASS1"), argv[1]);
  TFile *fileOut = new TFile(tempOut,"RECREATE");
  TTree *Tout = new TTree("pass1", "pass1");

  // Variables
  Tout->Branch("pmt0", &pmt[0], "pmt0/D");
  Tout->Branch("pmt1", &pmt[1], "pmt1/D");
  Tout->Branch("pmt2", &pmt[2], "pmt2/D");
  Tout->Branch("pmt3", &pmt[3], "pmt3/D");
  Tout->Branch("pmt4", &pmt[4], "pmt4/D");
  Tout->Branch("pmt5", &pmt[5], "pmt5/D");
  Tout->Branch("pmt6", &pmt[6], "pmt6/D");
  Tout->Branch("pmt7", &pmt[7], "pmt7/D");

  Tout->Branch("cathE0",  &cathodeEast[0],  "cathE0/D");
  Tout->Branch("cathE1",  &cathodeEast[1],  "cathE1/D");
  Tout->Branch("cathE2",  &cathodeEast[2],  "cathE2/D");
  Tout->Branch("cathE3",  &cathodeEast[3],  "cathE3/D");
  Tout->Branch("cathE4",  &cathodeEast[4],  "cathE4/D");
  Tout->Branch("cathE5",  &cathodeEast[5],  "cathE5/D");
  Tout->Branch("cathE6",  &cathodeEast[6],  "cathE6/D");
  Tout->Branch("cathE7",  &cathodeEast[7],  "cathE7/D");
  Tout->Branch("cathE8",  &cathodeEast[8],  "cathE8/D");
  Tout->Branch("cathE9",  &cathodeEast[9],  "cathE9/D");
  Tout->Branch("cathE10", &cathodeEast[10], "cathE10/D");
  Tout->Branch("cathE11", &cathodeEast[11], "cathE11/D");
  Tout->Branch("cathE12", &cathodeEast[12], "cathE12/D");
  Tout->Branch("cathE13", &cathodeEast[13], "cathE13/D");
  Tout->Branch("cathE14", &cathodeEast[14], "cathE14/D");
  Tout->Branch("cathE15", &cathodeEast[15], "cathE15/D");
  Tout->Branch("cathE16", &cathodeEast[16], "cathE16/D");
  Tout->Branch("cathE17", &cathodeEast[17], "cathE17/D");
  Tout->Branch("cathE18", &cathodeEast[18], "cathE18/D");
  Tout->Branch("cathE19", &cathodeEast[19], "cathE19/D");
  Tout->Branch("cathE20", &cathodeEast[20], "cathE20/D");
  Tout->Branch("cathE21", &cathodeEast[21], "cathE21/D");
  Tout->Branch("cathE22", &cathodeEast[22], "cathE22/D");
  Tout->Branch("cathE23", &cathodeEast[23], "cathE23/D");
  Tout->Branch("cathE24", &cathodeEast[24], "cathE24/D");
  Tout->Branch("cathE25", &cathodeEast[25], "cathE25/D");
  Tout->Branch("cathE26", &cathodeEast[26], "cathE26/D");
  Tout->Branch("cathE27", &cathodeEast[27], "cathE27/D");
  Tout->Branch("cathE28", &cathodeEast[28], "cathE28/D");
  Tout->Branch("cathE29", &cathodeEast[29], "cathE29/D");
  Tout->Branch("cathE30", &cathodeEast[30], "cathE30/D");
  Tout->Branch("cathE31", &cathodeEast[31], "cathE31/D");

  Tout->Branch("cathW0",  &cathodeWest[0],  "cathW0/D");
  Tout->Branch("cathW1",  &cathodeWest[1],  "cathW1/D");
  Tout->Branch("cathW2",  &cathodeWest[2],  "cathW2/D");
  Tout->Branch("cathW3",  &cathodeWest[3],  "cathW3/D");
  Tout->Branch("cathW4",  &cathodeWest[4],  "cathW4/D");
  Tout->Branch("cathW5",  &cathodeWest[5],  "cathW5/D");
  Tout->Branch("cathW6",  &cathodeWest[6],  "cathW6/D");
  Tout->Branch("cathW7",  &cathodeWest[7],  "cathW7/D");
  Tout->Branch("cathW8",  &cathodeWest[8],  "cathW8/D");
  Tout->Branch("cathW9",  &cathodeWest[9],  "cathW9/D");
  Tout->Branch("cathW10", &cathodeWest[10], "cathW10/D");
  Tout->Branch("cathW11", &cathodeWest[11], "cathW11/D");
  Tout->Branch("cathW12", &cathodeWest[12], "cathW12/D");
  Tout->Branch("cathW13", &cathodeWest[13], "cathW13/D");
  Tout->Branch("cathW14", &cathodeWest[14], "cathW14/D");
  Tout->Branch("cathW15", &cathodeWest[15], "cathW15/D");
  Tout->Branch("cathW16", &cathodeWest[16], "cathW16/D");
  Tout->Branch("cathW17", &cathodeWest[17], "cathW17/D");
  Tout->Branch("cathW18", &cathodeWest[18], "cathW18/D");
  Tout->Branch("cathW19", &cathodeWest[19], "cathW19/D");
  Tout->Branch("cathW20", &cathodeWest[20], "cathW20/D");
  Tout->Branch("cathW21", &cathodeWest[21], "cathW21/D");
  Tout->Branch("cathW22", &cathodeWest[22], "cathW22/D");
  Tout->Branch("cathW23", &cathodeWest[23], "cathW23/D");
  Tout->Branch("cathW24", &cathodeWest[24], "cathW24/D");
  Tout->Branch("cathW25", &cathodeWest[25], "cathW25/D");
  Tout->Branch("cathW26", &cathodeWest[26], "cathW26/D");
  Tout->Branch("cathW27", &cathodeWest[27], "cathW27/D");
  Tout->Branch("cathW28", &cathodeWest[28], "cathW28/D");
  Tout->Branch("cathW29", &cathodeWest[29], "cathW29/D");
  Tout->Branch("cathW30", &cathodeWest[30], "cathW30/D");
  Tout->Branch("cathW31", &cathodeWest[31], "cathW31/D");

  Tout->Branch("AnodeE", &AnodeE, "AnodeE/D");
  Tout->Branch("AnodeW", &AnodeW, "AnodeW/D");
  
  Tout->Branch("timeE", &timeE, "timeE/D");
  Tout->Branch("timeW", &timeW, "timeW/D");
  Tout->Branch("timeE_BB", &timeE_BB, "timeE_BB/D");
  Tout->Branch("timeW_BB", &timeW_BB, "timeW_BB/D");
  Tout->Branch("UBtime", &UBtime, "UBtime/D");
  Tout->Branch("UBtime_BB", &UBtime_BB, "UBtime_BB/D");
  Tout->Branch("twoFoldE", &twoFoldE, "twoFoldE/D");
  Tout->Branch("twoFoldW", &twoFoldW, "twoFoldW/D");
  
  Tout->Branch("xE", &xE, "xE/D");
  Tout->Branch("yE", &yE, "yE/D");
  Tout->Branch("xW", &xW, "xW/D");
  Tout->Branch("yW", &yW, "yW/D");

// Swank's addition   (1 of 2)  
  int xeRC,yeRC,xwRC,ywRC;
  WireChamberResponse * WCR = new WireChamberResponse();	
  Tout->Branch("xeRC", &xeRC, "xeRC/I"); //x east response class. 
  Tout->Branch("yeRC", &yeRC, "yeRC/I"); //y east response class... 
  Tout->Branch("xwRC", &xwRC, "xwRC/I");
  Tout->Branch("ywRC", &ywRC, "ywRC/I");
//end of Swank's addition (1 of 2)

  Tout->Branch("PID",  &PID,  "PID/I");
  Tout->Branch("type", &type, "type/I");
  Tout->Branch("side", &side, "side/I");
  Tout->Branch("posError", &posError, "posError/I");

  // East MWPC y                                                                  
  posPdc2[0]  =  76.20;
  posPdc2[1]  =  66.04;
  posPdc2[2]  =  55.88;
  posPdc2[3]  =  45.72;
  posPdc2[4]  =  35.56;
  posPdc2[5]  =  25.40;
  posPdc2[6]  =  15.24;
  posPdc2[7]  =   5.08;
  posPdc2[8]  =  -5.08;
  posPdc2[9]  = -15.24;
  posPdc2[10] = -25.40;
  posPdc2[11] = -35.56;
  posPdc2[12] = -45.72;
  posPdc2[13] = -55.88;
  posPdc2[14] = -66.04;
  posPdc2[15] = -76.20;

  // East MWPC x                                                                  
  posPdc2[16] =  76.20;
  posPdc2[17] =  66.04;
  posPdc2[18] =  55.88;
  posPdc2[19] =  45.72;
  posPdc2[20] =  35.56;
  posPdc2[21] =  25.40;
  posPdc2[22] =  15.24;
  posPdc2[23] =   5.08;
  posPdc2[24] =  -5.08;
  posPdc2[25] = -15.24;
  posPdc2[26] = -25.40;
  posPdc2[27] = -35.56;
  posPdc2[28] = -45.72;
  posPdc2[29] = -55.88;
  posPdc2[30] = -66.04;
  posPdc2[31] = -76.20;

  // West MWPC y                                                                  
  posPadc[0]  =  76.20;
  posPadc[1]  =  66.04;
  posPadc[2]  =  55.88;
  posPadc[3]  =  45.72;
  posPadc[4]  =  35.56;
  posPadc[5]  =  25.40;
  posPadc[6]  =  15.24;
  posPadc[7]  =   5.08;
  posPadc[8]  =  -5.08;
  posPadc[9]  = -15.24;
  posPadc[10] = -25.40;
  posPadc[11] = -35.56;
  posPadc[12] = -45.72;
  posPadc[13] = -55.88;
  posPadc[14] = -66.04;
  posPadc[15] = -76.20;

  // West MWPC x                                                                  
  posPadc[16] = -76.20;
  posPadc[17] = -66.04;
  posPadc[18] = -55.88;
  posPadc[19] = -45.72;
  posPadc[20] = -35.56;
  posPadc[21] = -25.40;
  posPadc[22] = -15.24;
  posPadc[23] =  -5.08;
  posPadc[24] =   5.08;
  posPadc[25] =  15.24;
  posPadc[26] =  25.40;
  posPadc[27] =  35.56;
  posPadc[28] =  45.72;
  posPadc[29] =  55.88;
  posPadc[30] =  66.04;
  posPadc[31] =  76.20;

  // Open input ntuple
  char tempIn[500];
  sprintf(tempIn, "/extern/mabrow05/ucna/rawdata/full%s.root", argv[1]);

  TFile *fileIn = new TFile(tempIn, "READ");
  TTree *Tin = (TTree*)(fileIn->Get("h1"));

  // Variables
  Tin->SetBranchAddress("Pdc30",  &Pdc30);
  Tin->SetBranchAddress("Pdc34",  &Pdc34);

  Tin->SetBranchAddress("Tdc016", &Tdc016);
  Tin->SetBranchAddress("Tdc017", &Tdc017);

  Tin->SetBranchAddress("Sis00",  &Sis00);

  Tin->SetBranchAddress("Qadc0",  &Qadc[0]);
  Tin->SetBranchAddress("Qadc1",  &Qadc[1]);
  Tin->SetBranchAddress("Qadc2",  &Qadc[2]);
  Tin->SetBranchAddress("Qadc3",  &Qadc[3]);
  Tin->SetBranchAddress("Qadc4",  &Qadc[4]);
  Tin->SetBranchAddress("Qadc5",  &Qadc[5]);
  Tin->SetBranchAddress("Qadc6",  &Qadc[6]);
  Tin->SetBranchAddress("Qadc7",  &Qadc[7]);

  Tin->SetBranchAddress("Pdc20",  &Pdc2[0]);
  Tin->SetBranchAddress("Pdc21",  &Pdc2[1]);
  Tin->SetBranchAddress("Pdc22",  &Pdc2[2]);
  Tin->SetBranchAddress("Pdc23",  &Pdc2[3]);
  Tin->SetBranchAddress("Pdc24",  &Pdc2[4]);
  Tin->SetBranchAddress("Pdc25",  &Pdc2[5]);
  Tin->SetBranchAddress("Pdc26",  &Pdc2[6]);
  Tin->SetBranchAddress("Pdc27",  &Pdc2[7]);
  Tin->SetBranchAddress("Pdc28",  &Pdc2[8]);
  Tin->SetBranchAddress("Pdc29",  &Pdc2[9]);
  Tin->SetBranchAddress("Pdc210", &Pdc2[10]);
  Tin->SetBranchAddress("Pdc211", &Pdc2[11]);
  Tin->SetBranchAddress("Pdc212", &Pdc2[12]);
  Tin->SetBranchAddress("Pdc213", &Pdc2[13]);
  Tin->SetBranchAddress("Pdc214", &Pdc2[14]);
  Tin->SetBranchAddress("Pdc215", &Pdc2[15]);
  Tin->SetBranchAddress("Pdc216", &Pdc2[16]);
  Tin->SetBranchAddress("Pdc217", &Pdc2[17]);
  Tin->SetBranchAddress("Pdc218", &Pdc2[18]);
  Tin->SetBranchAddress("Pdc219", &Pdc2[19]);
  Tin->SetBranchAddress("Pdc220", &Pdc2[20]);
  Tin->SetBranchAddress("Pdc221", &Pdc2[21]);
  Tin->SetBranchAddress("Pdc222", &Pdc2[22]);
  Tin->SetBranchAddress("Pdc223", &Pdc2[23]);
  Tin->SetBranchAddress("Pdc224", &Pdc2[24]);
  Tin->SetBranchAddress("Pdc225", &Pdc2[25]);
  Tin->SetBranchAddress("Pdc226", &Pdc2[26]);
  Tin->SetBranchAddress("Pdc227", &Pdc2[27]);
  Tin->SetBranchAddress("Pdc228", &Pdc2[28]);
  Tin->SetBranchAddress("Pdc229", &Pdc2[29]);
  Tin->SetBranchAddress("Pdc230", &Pdc2[30]);
  Tin->SetBranchAddress("Pdc231", &Pdc2[31]);

  Tin->SetBranchAddress("Padc0",  &Padc[0]);
  Tin->SetBranchAddress("Padc1",  &Padc[1]);
  Tin->SetBranchAddress("Padc2",  &Padc[2]);
  Tin->SetBranchAddress("Padc3",  &Padc[3]);
  Tin->SetBranchAddress("Padc4",  &Padc[4]);
  Tin->SetBranchAddress("Padc5",  &Padc[5]);
  Tin->SetBranchAddress("Padc6",  &Padc[6]);
  Tin->SetBranchAddress("Padc7",  &Padc[7]);
  Tin->SetBranchAddress("Padc8",  &Padc[8]);
  Tin->SetBranchAddress("Padc9",  &Padc[9]);
  Tin->SetBranchAddress("Padc10", &Padc[10]);
  Tin->SetBranchAddress("Padc11", &Padc[11]);
  Tin->SetBranchAddress("Padc12", &Padc[12]);
  Tin->SetBranchAddress("Padc13", &Padc[13]);
  Tin->SetBranchAddress("Padc14", &Padc[14]);
  Tin->SetBranchAddress("Padc15", &Padc[15]);
  Tin->SetBranchAddress("Padc16", &Padc[16]);
  Tin->SetBranchAddress("Padc17", &Padc[17]);
  Tin->SetBranchAddress("Padc18", &Padc[18]);
  Tin->SetBranchAddress("Padc19", &Padc[19]);
  Tin->SetBranchAddress("Padc20", &Padc[20]);
  Tin->SetBranchAddress("Padc21", &Padc[21]);
  Tin->SetBranchAddress("Padc22", &Padc[22]);
  Tin->SetBranchAddress("Padc23", &Padc[23]);
  Tin->SetBranchAddress("Padc24", &Padc[24]);
  Tin->SetBranchAddress("Padc25", &Padc[25]);
  Tin->SetBranchAddress("Padc26", &Padc[26]);
  Tin->SetBranchAddress("Padc27", &Padc[27]);
  Tin->SetBranchAddress("Padc28", &Padc[28]);
  Tin->SetBranchAddress("Padc29", &Padc[29]);
  Tin->SetBranchAddress("Padc30", &Padc[30]);
  Tin->SetBranchAddress("Padc31", &Padc[31]);

  Tin->SetBranchAddress("S83028", &S83028);
  Tin->SetBranchAddress("S8200", &S8200);
  Tin->SetBranchAddress("Clk0", &Clk0);
  Tin->SetBranchAddress("Clk1", &Clk1);
  Tin->SetBranchAddress("Clk2", &Clk2);
  Tin->SetBranchAddress("Clk3", &Clk3);

  Tin->SetBranchAddress("Pdc38",  &Pdc38);
  Tin->SetBranchAddress("Pdc39",  &Pdc39);
  Tin->SetBranchAddress("Pdc310", &Pdc310);
  Tin->SetBranchAddress("Pdc311", &Pdc311);

  Tin->SetBranchAddress("Qadc9",  &Qadc9);
  Tin->SetBranchAddress("Tdc019", &Tdc019);
  Tin->SetBranchAddress("Pdc313", &Pdc313);
  Tin->SetBranchAddress("Pdc315", &Pdc315);
  Tin->SetBranchAddress("Qadc8",  &Qadc8);
  Tin->SetBranchAddress("Tdc018", &Tdc018);
  Tin->SetBranchAddress("Qadc10", &Qadc10);
  Tin->SetBranchAddress("Tdc020", &Tdc020);

  int nEvents = Tin->GetEntries();
  cout << "... Processing nEvents = " << nEvents << endl;

  // Loop over events
  for (int i=0; i<nEvents; i++) {
    Tin->GetEvent(i);
    int iSis00 = (int) Sis00;

    // Calculate pedestal-subtracted PMT QADC values
    for (int j=0; j<8; j++) {
      pmt[j] = ((double) Qadc[j]) - pedQadc[j];
    }

    // Calculate pedestal-subtracted MWPC cathode PADC values
    for (int j=0; j<32; j++) {
      cathodeEast[j] = ((double) Pdc2[j]) - pedPdc2[j];
      cathodeWest[j] = ((double) Padc[j]) - pedPadc[j];
    }

    // Calculate pedestal-subtracted MWPC Anode PADC values
    //Took this out for now. We are comparing against a cut instead...
    //AnodeE = ((double) Pdc30) - pedPdc30;
    //AnodeW = ((double) Pdc34) - pedPdc34;

    // UCN monitor events
    bool UCNMonitorTrigger = false;
    if ( (iSis00 == 260) || (iSis00 == 516) || (iSis00 == 1028) || (iSis00 == 2052) ) {
      UCNMonitorTrigger = true;
    }

    // Events with a muon hit
    bool muonHitEast = false;
    bool muonHitWest = false;
    if ( (((double) Qadc8)  > cutEastBackingVetoQADC) ||
         (((double) Tdc018) > cutEastBackingVetoTDC)  ||
         (((double) Qadc9)  > cutEastTopVetoQADC)     ||
         (((double) Tdc019) > cutEastTopVetoTDC)      ||
         (((double) Pdc313) > cutEastDriftTubeTAC) ) {
      muonHitEast = true;
    }
    if ( (((double) Qadc10) > cutWestBackingVetoQADC) ||
         (((double) Tdc020) > cutWestBackingVetoTDC)  ||
         (((double) Pdc315) > cutWestDriftTubeTAC) ) {
      muonHitWest = true;
    }

    // LED trigger events
    bool LEDTrigger = false;
    if ( (iSis00 == 128) || (iSis00 == 129) || (iSis00 == 130) || (iSis00 == 131) || (iSis00 == 163) ) {
      LEDTrigger = true;
    }

    // Bi pulser trigger events
    bool bismuthPulser = false;
    if (iSis00 == 32) bismuthPulser = true;

    // Scintillator events
    bool mwpcHitEast = false;
    bool mwpcHitWest = false;
    bool scintillatorHitEast = false;
    bool scintillatorHitWest = false;
    bool scintillatorHitBoth = false;
    bool coincidenceEast = false;
    bool coincidenceWest = false;
    bool scintillatorHitFirstEast = false;
    bool scintillatorHitFirstWest = false;
    bool scintillatorHitBothBad = false;
    bool triggerEast = false;
    bool triggerWest = false;

    if ( ((double) Pdc30)  > cutEastAnode) mwpcHitEast = true;
    if ( ((double) Pdc34)  > cutWestAnode) mwpcHitWest = true;

    double timeEastTwoFold = ((double) Tdc016)*tdcChannelToTime;
    double timeWestTwoFold = ((double) Tdc017)*tdcChannelToTime;
    if (timeEastTwoFold > 0.*tdcChannelToTime) scintillatorHitEast = true;
    if (timeWestTwoFold > 0.*tdcChannelToTime) scintillatorHitWest = true;
    if (scintillatorHitEast && scintillatorHitWest) {
      scintillatorHitBoth = true;
      if ( (timeEastTwoFold > cutEastTwoFold) && (timeWestTwoFold < cutWestTwoFold) ) {
        scintillatorHitFirstEast = true;
      }
      if ( (timeEastTwoFold < cutEastTwoFold) && (timeWestTwoFold > cutWestTwoFold) ) {
        scintillatorHitFirstWest = true;
      }
      if ( (timeEastTwoFold > cutEastTwoFold) && (timeWestTwoFold > cutWestTwoFold) ) {
        scintillatorHitBothBad = true;
      }
      if ( (timeEastTwoFold < cutEastTwoFold) && (timeWestTwoFold < cutWestTwoFold) ) {
        scintillatorHitBothBad = true;
      }
    }

    if (mwpcHitEast && scintillatorHitEast) coincidenceEast = true;
    if (mwpcHitWest && scintillatorHitWest) coincidenceWest = true;

    // Event PID logic
    if (UCNMonitorTrigger) PID = 5;
    else if (LEDTrigger) PID = 3;
    else if (muonHitEast || muonHitWest) PID = 2;
    else if (bismuthPulser) PID = 4;
    else if ( (scintillatorHitEast && !mwpcHitEast && !scintillatorHitWest && !mwpcHitWest) ||
              (!scintillatorHitEast && !mwpcHitEast && scintillatorHitWest && !mwpcHitWest) ) PID = 0;
    else if ( (coincidenceEast && !scintillatorHitWest && !mwpcHitWest) ||
              (!scintillatorHitEast && !mwpcHitEast && coincidenceWest) ||
              (coincidenceEast && coincidenceWest && !scintillatorHitBothBad) ||
              (coincidenceEast && !scintillatorHitWest && mwpcHitWest) ||
              (!scintillatorHitEast && mwpcHitEast && coincidenceWest) ) PID = 1;
    else PID = 6;

    type = -1;
    side = -1;
    if (PID == 1) {
      if (coincidenceEast && !scintillatorHitWest && !mwpcHitWest) {
        type = 0;
        side = 0;
      }
      if (!scintillatorHitEast && !mwpcHitEast && coincidenceWest) {
        type = 0;
        side = 1;
      }
      if (coincidenceEast && coincidenceWest && !scintillatorHitBothBad) {
        type = 1;
        if (scintillatorHitFirstEast) side = 0;
        if (scintillatorHitFirstWest) side = 1;
      }
      if (coincidenceEast && !scintillatorHitWest && mwpcHitWest) {
        type = 2;
        side = 0;
      }
      if (!scintillatorHitEast && mwpcHitEast && coincidenceWest) {
        type = 2;
        side = 1;
      }
    }

    // Calculate MWPC positions for electron events
    xE = -100.;
    yE = -100.;
    xW = -100.;
    yW = -100.;
    posError = 0;

    if (PID == 1) {

      /*
      // Find cathode wire with maximum signal on each plane
      double cathodeEastMax_x = -1.0;
      double cathodeEastMax_y = -1.0;
      double cathodeWestMax_x = -1.0;
      double cathodeWestMax_y = -1.0;
      int intCathodeEastMax_x = -1;
      int intCathodeEastMax_y = -1;
      int intCathodeWestMax_x = -1;
      int intCathodeWestMax_y = -1;
      for (int j=0; j<16; j++) {
        if (cathodeEast[j+16] > cathodeEastMax_x) {
          cathodeEastMax_x = cathodeEast[j+16];
          intCathodeEastMax_x = j+16;
        }
        if (cathodeEast[j] > cathodeEastMax_y) {
          cathodeEastMax_y = cathodeEast[j];
          intCathodeEastMax_y = j;
        }
        if (cathodeWest[j+16] > cathodeWestMax_x) {
          cathodeWestMax_x = cathodeWest[j+16];
          intCathodeWestMax_x = j+16;
        }
        if (cathodeWest[j] > cathodeWestMax_y) {
          cathodeWestMax_y = cathodeWest[j];
          intCathodeWestMax_y = j;
        }
      }

      // East x
      double xMWPCEast = 0.;
      double xMWPCEastSum = 0.;
      if (mwpcHitEast) {
        if (intCathodeEastMax_x == 16) {
          xMWPCEast    = cathodeEast[16]*posPdc2[16] + cathodeEast[17]*posPdc2[17];
          xMWPCEastSum = cathodeEast[16] + cathodeEast[17];
        }
        else if (intCathodeEastMax_x == 31) {
          xMWPCEast    = cathodeEast[30]*posPdc2[30] + cathodeEast[31]*posPdc2[31];
          xMWPCEastSum = cathodeEast[30] + cathodeEast[31];
        }
        else {
          int n = intCathodeEastMax_x;
          xMWPCEast    = cathodeEast[n-1]*posPdc2[n-1] + cathodeEast[n]*posPdc2[n] +
                         cathodeEast[n+1]*posPdc2[n+1];
          xMWPCEastSum = cathodeEast[n-1] + cathodeEast[n] + cathodeEast[n+1];
        }
      }
      if (xMWPCEastSum > 0.) {
        xMWPCEast = xMWPCEast / xMWPCEastSum;
      }
      else {
        xMWPCEast = -100.;
        posError = 1;
      }

      // East y      
      double yMWPCEast = 0.;
      double yMWPCEastSum = 0.;
      if (mwpcHitEast) {
        if (intCathodeEastMax_y == 0) {
          yMWPCEast    = cathodeEast[0]*posPdc2[0] + cathodeEast[1]*posPdc2[1];
          yMWPCEastSum = cathodeEast[0] + cathodeEast[0];
        }
        else if (intCathodeEastMax_y == 15) {
          yMWPCEast    = cathodeEast[14]*posPdc2[14] + cathodeEast[15]*posPdc2[15];
          yMWPCEastSum = cathodeEast[14] + cathodeEast[15];
        }
        else {
          int n = intCathodeEastMax_y;
          yMWPCEast    = cathodeEast[n-1]*posPdc2[n-1] + cathodeEast[n]*posPdc2[n] +
                         cathodeEast[n+1]*posPdc2[n+1];
          yMWPCEastSum = cathodeEast[n-1] + cathodeEast[n] + cathodeEast[n+1];
        }
      }
      if (yMWPCEastSum > 0.) {
        yMWPCEast = yMWPCEast / yMWPCEastSum;
      }
      else {
	yMWPCEast = -100.;
        posError = 1;
      }

      // West x
      double xMWPCWest = 0.;
      double xMWPCWestSum = 0.;
      if (mwpcHitWest) {
        if (intCathodeWestMax_x == 16) {
          xMWPCWest    = cathodeWest[16]*posPadc[16] + cathodeWest[17]*posPadc[17];
          xMWPCWestSum = cathodeWest[16] + cathodeWest[17];
        }
        else if (intCathodeWestMax_x == 31) {
          xMWPCWest    = cathodeWest[30]*posPadc[30] + cathodeWest[31]*posPadc[31];
          xMWPCWestSum = cathodeWest[30] + cathodeWest[31];
        }
        else {
          int n = intCathodeWestMax_x;
          xMWPCWest    = cathodeWest[n-1]*posPadc[n-1] + cathodeWest[n]*posPadc[n] +
	                 cathodeWest[n+1]*posPadc[n+1];
          xMWPCWestSum = cathodeWest[n-1] + cathodeWest[n] + cathodeWest[n+1];
        }
      }
      if (xMWPCWestSum > 0.) {
        xMWPCWest = xMWPCWest / xMWPCWestSum;
      }
      else {
	xMWPCWest = -100.;
	posError = 1;
      }

      // West y
      double yMWPCWest = 0.;
      double yMWPCWestSum = 0.;
      if (mwpcHitWest) {
        if (intCathodeWestMax_y == 0) {
          yMWPCWest    = cathodeWest[0]*posPadc[0] + cathodeWest[1]*posPadc[1];
          yMWPCWestSum = cathodeWest[0] + cathodeWest[1];
        }
        else if (intCathodeWestMax_x == 15) {
          yMWPCWest    = cathodeWest[14]*posPadc[14] + cathodeWest[15]*posPadc[15];
          yMWPCWestSum = cathodeWest[14] + cathodeWest[15];
        }
        else {
          int n = intCathodeWestMax_y;
          yMWPCWest    = cathodeWest[n-1]*posPadc[n-1] + cathodeWest[n]*posPadc[n] +
                         cathodeWest[n+1]*posPadc[n+1];
          yMWPCWestSum = cathodeWest[n-1] + cathodeWest[n] + cathodeWest[n+1];
        }
      }
      if (yMWPCWestSum > 0.) {
        yMWPCWest = yMWPCWest / yMWPCWestSum;
      }
      else {
        yMWPCWest = -100.;
        posError = 1;
      }
      */

/*  Swank's addition to the UK PA.  (2 of 2)
	 this seems like a good spot since its around the position calculation. 
*/
		//changing the length 32 double array to length 16  float array.  using variables defined in WCR. 
		for(int j = 0; j<16; j++)
		{
		  WCR->cathex[j]=(float)cathodeEast[16+j];
		  WCR->cathey[j]=(float)cathodeEast[j];
		  WCR->cathwx[j]=(float)cathodeWest[16+j];
		  WCR->cathwy[j]=(float)cathodWest[j];
		} 
		   		  
		xeRC=WCR->ResponseType(WCR->cathex);   //response class x co-ordinate, East. 
    	yeRC=WCR->ResponseType(WCR->cathey);   //response class y co-ordinate, East. 
		xwRC=WCR->ResponseType(WCR->cathwx);   //response class x co-ordinate, West. 
		ywRC=WCR->ResponseType(WCR->cathwy);   //response class y co-ordinate, West.
/*
     End of swank's addtion (2 of 2)
*/

      double xMWPCEast = 0.;
      double yMWPCEast = 0.;
      double xMWPCWest = 0.;
      double yMWPCWest = 0.;
      double xMWPCEastSum = 0.;
      double yMWPCEastSum = 0.;
      double xMWPCWestSum = 0.;
      double yMWPCWestSum = 0.;

	  

		
      for (int j=0; j<16; j++) {

        if (cathodeEast[j+16] > 100.) {
          xMWPCEast += cathodeEast[j+16]*posPdc2[j+16];
          xMWPCEastSum += cathodeEast[j+16];
        }
        if (cathodeEast[j] > 100.) {
          yMWPCEast += cathodeEast[j]*posPdc2[j];
          yMWPCEastSum += cathodeEast[j];
        }
        if (cathodeWest[j+16] > 100.) {
          xMWPCWest += cathodeWest[j+16]*posPadc[j+16];
          xMWPCWestSum += cathodeWest[j+16];
        }
        if (cathodeWest[j] > 100.) {
          yMWPCWest += cathodeWest[j]*posPadc[j];
          yMWPCWestSum += cathodeWest[j];
        }

      }

      if (xMWPCEastSum > 0.)
        xMWPCEast = xMWPCEast / xMWPCEastSum;
      if (yMWPCEastSum > 0.)
        yMWPCEast = yMWPCEast / yMWPCEastSum;
      if (xMWPCWestSum > 0.)
        xMWPCWest = xMWPCWest / xMWPCWestSum;
      if (yMWPCWestSum > 0.)
        yMWPCWest = yMWPCWest / yMWPCWestSum;

      if (mwpcHitEast && xMWPCEastSum > 0.)
        xE = xMWPCEast * positionProjection;
      if (mwpcHitEast && yMWPCEastSum > 0.)
        yE = yMWPCEast * positionProjection;
      if (mwpcHitWest && xMWPCWestSum > 0.)
        xW = xMWPCWest * positionProjection;
      if (mwpcHitWest && yMWPCWestSum > 0.)
        yW = yMWPCWest * positionProjection;

      // Pass timing information to next root file
      timeE = Clk0*scalerCountsToTime;
      timeW = Clk1*scalerCountsToTime;
      timeE_BB = Clk2*scalerCountsToTime;
      timeW_BB = Clk3*scalerCountsToTime;
      UBtime = S83028*scalerCountsToTime;
      UBtime_BB = S8200*scalerCountsToTime;
      twoFoldE = Tdc016;
      twoFoldW = Tdc017;
    }

    Tout->Fill();
	

  }

  // Write output ntuple
  fileOut->Write();
  fileOut->Close();

  return 0;
}
