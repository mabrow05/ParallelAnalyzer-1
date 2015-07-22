//Holds the smeared source peaks which can be changed in single location

#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdio>
#include <vector>
#include <cstdlib>
#include <string>


using namespace std;

const static double peakCe = 80.5;
const static double peakSn = 317.8;
const static double peakBiLow = 448.8;
const static double peakBiHigh = 926.;

vector < vector <double> > returnPeaks(int srcPeriod) {
  vector < vector <double> > peaks;
  peaks.resize( 8, vector <double> (4,0.));
  ifstream peakfile;
  char filename[500];
  sprintf(filename, "../smeared_peaks/smear_%i.dat",srcPeriod);
  peakfile.open(filename);
  
  for (int pmt = 0; pmt<8; pmt++) {
    peakfile >> peaks[pmt][0] >> peaks[pmt][1] >> peaks[pmt][2] >> peaks[pmt][3];
  }
  return peaks;
}


