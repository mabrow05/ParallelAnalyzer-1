//Holds the smeared source peaks which can be changed in single location

#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdio>
#include <vector>
#include <cstdlib>
#include <string>


using namespace std;

// These are the weighted true peaks from the KLM CE sources
const static double peakCe = 131.5956;//80.5;
const static double peakSn = 368.4938;//317.8;
const static double peakBiLow = 501.420;//448.8;
const static double peakBiHigh = 993.789;//926.;


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

vector < vector <double> > EQ2EtrueFit(int srcPeriod) {
  vector < vector <double> > fits;
  fits.resize( 8, vector <double> (3,0.));
  ifstream sourcefile;
  char filename[500];
  sprintf(filename, "../smeared_peaks/EQ2Etrue_srcPeriod_%i.dat",srcPeriod);
  sourcefile.open(filename);
  
  for (int pmt = 0; pmt<8; pmt++) {
    sourcefile >> fits[pmt][0] >> fits[pmt][1] >> fits[pmt][2];
  }
  return fits;
}


