#include "triggerMap.hh"

#include <string>
#include <fstream>
#include <cstdlib>
#include "MBUtils.hh"

#include <TString.h>

TriggerMap::TriggerMap(double bin_width, Int_t nparams) : binWidth(bin_width), nParams(nparams), func(NULL) {
  
  int nBinsHold = (int)(110./binWidth);
  nBinsXY = ((int)nBinsHold%2)==0 ? (int)nBinsHold+1 : (int)nBinsHold;
  nBinsTotal = nBinsXY*nBinsXY;

  xyBinLower.resize(nBinsXY,0.);
  xyBinUpper.resize(nBinsXY,0.);
  xyBinCenter.resize(nBinsXY,0.);
  

  triggMap.resize(nBinsXY, std::vector <double> (nBinsXY, std::vector <double> (nParams,0.)));
  
  setBinValues();

};

TriggerMap::~TriggerMap() {
  delete func;
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
  std::string file = getenv("TRIGGER_FUNCTIONS") + "/trigger_functions_XePeriod_"+itos(XePeriod)+"_"ftos(binWidth)+"mm" + ".dat";
  ifstream infile(file);

  int binValx, binValy;
  std::vector <double> p(8,0.);
  
  while (infile >> binValx >> binValy >> p[0] >> p[1] >> p[2] >> p[3] >> p[4] >> p[5] >> p[6] >> p[7])
    {
      setTriggerMapPoint(getBinNumber(binValx), getBinNumber(binValy), p);
    }
  infile.close();

  std::cout << "Read in Trigger Map from Xe Period " << XePeriod << std::endl;
};

void TriggerMap::setTriggerFunc(TF1* f) {
  if (func) delete func;
  func = new TF1(*f); 
  
};
void TriggerMap::setTriggerFunc(TF1 f) {
  if (func) delete func;
  func = new TF1(f);
};

double TriggerMap::returnTriggerProbability(double x, double y, double En) {
  
  if (!func) { std::cout << "No Function Set..\n"; return 0.; }

  else {
    std::vector <Double_t> pars = triggMap[getBinNumber(x)][getBinNumber(y)];
    return func->EvalPar(En, &pars[0]);
  }
}

int TriggerMap::getBinNumber(double pos) {
  for (int m=0; m<nBinsXY; m++) {
    if ( (pos >= xyBinLower[m]) && (pos < xyBinUpper[m]) )  return = m;
  }
};

void TriggerMap::setTriggerMapPoint(int xBin, int yBin, std::vector <double> params) {
  triggMap[xBin][yBin] = params;
};


void TriggerMap::writeTriggerMap(int XeRunPeriod) {
  
  std::string file = getenv("TRIGGER_FUNCTIONS") + "/trigger_functions_XePeriod_"+itos(XePeriod)+"_"ftos(binWidth)+"mm" + ".dat";
  ofstream ofile(file);
  
  for (Int_t xbin=0; xbin<nBinsXY; xbin++) {
    for (Int_t ybin=0; ybin<nBinsXY; ybin++) {

      ofile << xyBinCenter[xbin] << " " 
	    << xyBinCenter[ybin] << " ";
      
      for (Int_t n=0; n<nParams; n++) {
	
	ofile << triggMap[xbin][ybin][n];
	if (n!=nParams-1) ofile << " ";
      }
      
      if (xbin!=nBinsXY-1) ofile << std::endl;
    }	   
  }
}
