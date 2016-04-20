#include "triggerMap.hh"

#include <string>
#include <fstream>
#include <cstdlib>
#include "MBUtils.hh"

#include <TString.h>

TriggerMap::TriggerMap(double bin_width) : binWidth(bin_width), nParams(8); {
  
  int nBinsHold = (int)(120./binWidth);
  nBinsXY = ((int)nBinsHold%2)==0 ? (int)nBinsHold : (int)nBinsHold+1;
  nBinsTotal = nBinsXY*nBinsXY;

  xBinLower.resize(nBinsXY,0.);
  xBinUpper.resize(nBinsXY,0.);
  xBinCenter.resize(nBinsXY,0.);
  yBinLower.resize(nBinsXY,0.);
  yBinUpper.resize(nBinsXY,0.);
  yBinCenter.resize(nBinsXY,0.);

  triggMap.resize(nBinsXY, std::vector <double> (nBinsXY, std::vector <double> (nParams,0.)));
  
  setBinValues();

};

TriggerMap::~TriggerMap() {
  
};

void TriggerMap::setBinValues() {
  for (int k=0; k<nBinsXY; k++) {
    xyBinLower[k]     = -(double)nBinsXY*binWidth/2. + ((double) k)*binWidth;
    xyBinUpper[k]     = -(double)nBinsXY*binWidth/2. + ((double) k)*binWidth + binWidth;
    xyBinCenter[k]    = (xBinLower[k] + xBinUpper[k])/2.;
    
    //cout << xBinLower[k] << " " << intXBinCenter[k] << " " << xBinUpper[k] << endl;
  }
};

void TriggerMap::readTriggerMap(int XeRunPeriod) {
  XePeriod = XeRunPeriod;
  std::string file = getenv("TRIGGER_FUNCTIONS") + "/trigger_functions_XePeriod_"+itos(XePeriod) + ".dat";
  ifstream infile(file);

  int binValx, binValy;
  std::vector <double> p(8,0.);
  
  while (infile >> binValx >> binValy >> p[0] >> p[1] >> p[2] >> p[3] >> p[4] >> p[5] >> p[6] >> p[7])
    {
      setTriggerMapPoint(binValx, binValy, p);
    }

  std::cout << "Read in Trigger Map from Xe Period " << XePeriod << std::endl;
};

int TriggerMap::getBinNumber(double pos) {
  for (int m=0; m<nBinsXY; m++) {
    if ( (pos >= xyBinLower[m]) && (pos < xyBinUpper[m]) )  return = m;
  }
};

void TriggerMap::setTriggerMapPoint(int xBin, int yBin, std::vector <double> params) {
  triggMap[xBin][yBin] = params;
};

