#include "positionMapHandler.hh"

#include <string>
#include <fstream>
#include <cstdlib>
#include "MBUtils.hh"




PositionMap::PositionMap(Double_t bin_width) : binWidth(bin_width), nParams(8); {
  
  Int_t nBinsHold = (int)(120./binWidth);
  nBinsXY = ((int)nBinsHold%2)==0 ? (int)nBinsHold : (int)nBinsHold+1;
  nBinsTotal = nBinsXY*nBinsXY;

  xBinLower.resize(nBinsXY,0.);
  xBinUpper.resize(nBinsXY,0.);
  xBinCenter.resize(nBinsXY,0.);
  yBinLower.resize(nBinsXY,0.);
  yBinUpper.resize(nBinsXY,0.);
  yBinCenter.resize(nBinsXY,0.);

  triggMap.resize(nBinsXY, std::vector <Double_t> (nBinsXY, std::vector <Double_t> (nParams,0.)));
  
  setBinValues();
  
  funcs.push_back(new TF1("func0", "-0.5(1.-x)*(1.-x)*x", 0., 1.));
  funcs.push_back(new TF1("func1", "(1.+x-3./2.*x*x)*(1.-x)", 0., 1.));
  funcs.push_back(new TF1("func2", "(1+(1-x)-3./2.*(1.-x)*(1-.x))*x", 0., 1.));
  funcs.push_back(new TF1("func3", "-0.5(1.-x)*x*x", 0., 1.));
};

PositionMap::~PositionMap() {
  for (UInt_t i=0; i<funcs.size(); i++) {
    if (funcs[i]) delete funcs[i];
  }
};

void PositionMap::setBinValues() {
  for (Int_t k=0; k<nBinsXY; k++) {
    xyBinLower[k]     = -(Double_t)nBinsXY*binWidth/2. + ((Double_t) k)*binWidth;
    xyBinUpper[k]     = -(Double_t)nBinsXY*binWidth/2. + ((Double_t) k)*binWidth + binWidth;
    xyBinCenter[k]    = (xBinLower[k] + xBinUpper[k])/2.;
    
    //cout << xBinLower[k] << " " << intXBinCenter[k] << " " << xBinUpper[k] << endl;
  }
};

void PositionMap::readPositionMap(Int_t XeRunPeriod) {
  XePeriod = XeRunPeriod;
  std::string file = getenv("POSITION_MAPS") + "/position_map_" + itos(XePeriod) + "_" + "RC_123_" + ftos(binWidth) +  "mm.dat";
  ifstream infile(file);
  
  Int_t binValx, binValy;
  std::vector <Double_t> p(8,0.);
  
  while (infile >> binValx >> binValy >> p[0] >> p[1] >> p[2] >> p[3] >> p[4] >> p[5] >> p[6] >> p[7])
    {
      setPositionMapPoint(binValx, binValy, p);
    }

  std::cout << "Read in Position Map from Xe Period " << XePeriod << std::endl;
};

Int_t PositionMap::getBinNumber(Double_t pos) {
  for (Int_t m=0; m<nBinsXY; m++) {
    if ( (pos >= xyBinLower[m]) && (pos < xyBinUpper[m]) )  return = m;
  }
};

void PositionMap::setPositionMapPoint(Int_t xBin, Int_t yBin, std::vector <Double_t> params) {
  posMap[xBin][yBin] = params;
};




Double_t PositionMap::getInterpolatedEta(Double_t x, Double_t y) {

  std::vector <Double_t> xInterp;
  std::vector <Double_t> yInterp;

  Int_t xBin = x > getBinCenter(getBinNumber(x)) ? getBinNumber(x) : getBinNumber(x)-1;
  Int_t yBin = y > getBinCenter(getBinNumber(y)) ? getBinNumber(y) : getBinNumber(y)-1;

  Double_t xp = (x-xBin)/binWidth;
  Double_t yp = (y-yBin)/binWidth;

  Int_t xx=0, yy=0;

  Double_t val = 0.;

  for (int xb=xBin-1; xb<=xBin+2; xb++) {
    for (int yb=yBin-1; yb<=yBin+2; yb++) {
      val+=funcs[xx]->Eval(xp)*funcs[yy]->Eval(yp)*posMap[xb][yb];
      yy++;
    }
    xx++;
  }

  return val;

};

