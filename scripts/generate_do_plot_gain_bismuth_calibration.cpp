#include <iostream>
#include <fstream>
using namespace std;

int main()
{

  // Output command file
  char tempOut[500];
  sprintf(tempOut, "do_plot_gain_bismuth_calibration");
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
      outFile << "root -b -q '../gain_bismuth/plot_gain_bismuth.C(\"" << runNumber << "\")'" << endl;
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
      outFile << "root -b -q '../gain_bismuth/plot_gain_bismuth.C(\"" << runNumber << "\")'" << endl;
    }
  }

  outFile.close();

  return 0;
}
