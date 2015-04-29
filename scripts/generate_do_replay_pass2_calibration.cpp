#include <iostream>
#include <fstream>
using namespace std;

int main()
{

  // Output command file
  char tempOut[500];
  sprintf(tempOut, "do_replay_pass2_calibration");
  ofstream outFile(tempOut);
  outFile << "#!/bin/csh" << endl;
  outFile << endl;

  // Source calibration runs
  int runNumber;
  char tempIn[500];
  for (int i=0; i<11; i++) {
    sprintf(tempIn, "../run_lists/Source_Calibration_Run_Period_%i.dat", i+1);
    cout << tempIn << endl;

    ifstream inFile(tempIn);
    while (!inFile.eof()) {
      inFile >> runNumber;
      if (inFile.fail()) break;
      cout << runNumber << endl;
      outFile << "~plaster/ucna/Analyzer/replay_pass2/replay_pass2.exe " << runNumber << endl;
    }
  }

  // Xenon calibration runs
  for (int i=0; i<7; i++) {
    sprintf(tempIn, "../run_lists/Xenon_Calibration_Run_Period_%i.dat", i+1);
    cout << tempIn << endl;

    ifstream inFile(tempIn);
    while (!inFile.eof()) {
      inFile >> runNumber;
      if (inFile.fail()) break;
      cout << runNumber << endl;
      outFile << "~plaster/ucna/Analyzer/replay_pass2/replay_pass2.exe " << runNumber << endl;
    }
  }

  outFile.close();

  return 0;
}
