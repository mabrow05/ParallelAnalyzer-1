#include "triggerMap.hh"

#include <string>
#include <fstream>
#include <cstdlib>

triggerMap::triggerMap(int bin_width) : binWidth(bin_width) {
  
  int nBinsHold = 120./binWidth;
  nBinsXY = ((int)nBinsHold%2)==0 ? (int)nBinsHold : (int)nBinsHold+1;
  
  xBinLower.resize(nBinsXY,0.);
  xBinUpper.resize(nBinsXY,0.);
  xBinCenter.resize(nBinsXY,0.);
  yBinLower.resize(nBinsXY,0.);
  yBinUpper.resize(nBinsXY,0.);
  yBinCenter.resize(nBinsXY,0.);

  triggMap.resize(nBinsXY, std::vector <double> (nBinsXY,0.));
  
  setBinValues();

};

triggerMap::~triggerMap() {
  
};

void triggerMap::setBinValues() {
  for (int k=0; k<nBinsXY; k++) {
    xyBinLower[k]     = -(double)nBinsXY*binWidth/2. + ((double) k)*binWidth;
    xyBinUpper[k]     = -(double)nBinsXY*binWidth/2. + ((double) k)*binWidth + binWidth;
    xyBinCenter[k]    = (xBinLower[k] + xBinUpper[k])/2.;
    
    //cout << xBinLower[k] << " " << intXBinCenter[k] << " " << xBinUpper[k] << endl;
  }
};

void triggerMap::readTriggerMap(int XeRunPeriod) {
  XePeriod = XeRunPeriod;
  std::string file = getenv



