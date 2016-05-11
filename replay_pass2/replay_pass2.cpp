//////////////////////////////////////////////////////////
// Applies Bi pulser gain monitoring corrections to 
// Scintillator PMT signals
//////////////////////////////////////////////////////////

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
#include "DataTree.hh"

#include "replay_pass2.h"


using namespace std;

int main(int argc, char *argv[])
{
  cout.setf(ios::fixed, ios::floatfield);
  cout.precision(12);

  cout << "Run " << argv[1] << " ..." << endl;
  cout << "... Applying Bi pulser gain corrections ..." << endl;

  // Read gain corrections file
  char tempFileGain[500];
  sprintf(tempFileGain, "%s/gain_bismuth_%s.dat",getenv("GAIN_BISMUTH"), argv[1]);
  cout << "... Reading: " << tempFileGain << endl;

  double fitMean[8], gainCorrection[8];
  ifstream fileGain(tempFileGain);
  for (int i=0; i<8; i++) {
    fileGain >> fitMean[i] >> gainCorrection[i];
  }
  cout << "...   PMT E1: " << gainCorrection[0] << endl;
  cout << "...   PMT E2: " << gainCorrection[1] << endl;
  cout << "...   PMT E3: " << gainCorrection[2] << endl;
  cout << "...   PMT E4: " << gainCorrection[3] << endl;
  cout << "...   PMT W1: " << gainCorrection[4] << endl;
  cout << "...   PMT W2: " << gainCorrection[5] << endl;
  cout << "...   PMT W3: " << gainCorrection[6] << endl;
  cout << "...   PMT W4: " << gainCorrection[7] << endl;

  // Open output ntuple
  char tempOut[500];
  sprintf(tempOut, "%s/replay_pass2_%s.root",getenv("REPLAY_PASS2"), argv[1]);
  //sprintf(tempOut, "replay_pass2_%s.root", argv[1]);
  DataTree *t = new DataTree();
  t->makeOutputTree(std::string(tempOut),"pass2");

  char tempIn[500];
  sprintf(tempIn, "%s/replay_pass1_%s.root", getenv("REPLAY_PASS1"),argv[1]);
  //sprintf(tempIn, "../replay_pass1/replay_pass1_%s.root", argv[1]);
  t->setupInputTree(std::string(tempIn),"pass1");

  int nEvents = t->getEntries();
  cout << "... Processing nEvents = " << nEvents << endl;

  // Loop over events
  for (Int_t i=0; i<nEvents; i++) {
    t->getEvent(i);

    // Apply gain correction factors
    t->ScintE.q1 = t->ScintE.q1*gainCorrection[0];
    t->ScintE.q2 = t->ScintE.q2*gainCorrection[1];
    t->ScintE.q3 = t->ScintE.q3*gainCorrection[2];
    t->ScintE.q4 = t->ScintE.q4*gainCorrection[3];

    t->ScintW.q1 = t->ScintW.q1*gainCorrection[4];
    t->ScintW.q2 = t->ScintW.q2*gainCorrection[5];
    t->ScintW.q3 = t->ScintW.q3*gainCorrection[6];
    t->ScintW.q4 = t->ScintW.q4*gainCorrection[7]; 

    t->fillOutputTree();
  }

  // Write output ntuple
  t->writeOutputFile();
  delete t;

  return 0;
}

