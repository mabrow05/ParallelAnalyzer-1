#include <vector>
#include <TF1.h>

class PositionMap {

public:

  PositionMap(Int_t bin_width); //Constructor
  ~PositionMap();

  void readPositionMap(Int_t XeRunPeriod); //Read in a trigger map
  void writePositionMap(Int_t XeRunPeriod); //Write trigger map
  void setPositionMapPoint(Int_t xBin, Int_t yBin, Double_t val); //Set a trigger map val

  Double_t getInterpolatedEta(Double_t x, Double_t y); //return interpolated value from position map

  Int_t getBinNumber(Double_t pos); //returns bin number for position

  Double_t getBinWidth() {return binWidth;}
  Int_t getNbins() {return nBinsTotal;}
  Int_t getNbinsXY() {return nBinsXY;}
  Double_t getBinUpper(Int_t bin) {return xyBinUpper[bin];}
  Double_t getBinLower(Int_t bin) {return xyBinLower[bin];}
  Double_t getBinCenter(Int_t bin) {return xyBinCenter[bin];}
  Int_t getCurrentXeRunPeriod() {return XePeriod;}

private:

  void setBinValues();

  Double_t binWidth; //x and y width of position bins
  Int_t nBinsXY, nBinsTotal;
  Int_t XePeriod;
  Int_t nParams;
  std::vector <Double_t> xyBinLower;
  std::vector <Double_t> xyBinUpper;
  std::vector <Double_t> xyBinCenter;
  std::vector <TF1*> funcs;
  //std::vector <Double_t> yBinLower;
  //std::vector <Double_t> yBinUpper;
  //std::vector <Double_t> yBinCenter;

  std::vector < std::vector < std::vector <Double_t> > > posMap;
 
};
