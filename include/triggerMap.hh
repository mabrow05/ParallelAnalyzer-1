//
#include <vector>
#include <TF1.h>

class TriggerMap {

public:

  TriggerMap(int bin_width); //Constructor
  ~TriggerMap();

  void readTriggerMap(int XeRunPeriod, Double_t nparams); //Read in a trigger map
  double returnTriggerProbability(double x, double y, double En);
  std::vector <double> getTriggerFunction(double x, double y); //returns the trigger function for the given event position
  void writeTriggerMap(int XeRunPeriod); //Write trigger map
  void setTriggerMapPoint(int xBin, int yBin, double val); //Set a trigger map val

  void setTriggerFunc(TF1* f);
  void setTriggerFunc(TF1 f);

  int getBinNumber(double pos); //returns bin number for position

  double getBinWidth() {return binWidth;}
  int getNbins() {return nBinsTotal;}
  int getNbinsXY() {return nBinsXY;}
  Double_t getBinUpper(Int_t bin) {return (bin >= 0 && bin < nBinsTotal) ? xyBinUpper[bin] : -1000. ;}
  Double_t getBinLower(Int_t bin) {return (bin >= 0 && bin < nBinsTotal) ? xyBinLower[bin] : -1000.;}
  Double_t getBinCenter(Int_t bin) {return (bin >= 0 && bin < nBinsTotal) ? xyBinCenter[bin] : -1000.;}
  Int_t getCurrentXeRunPeriod() {return XePeriod;}

private:

  void setBinValues();

  double binWidth; //x and y width of position bins
  int nBinsXY, nBinsTotal;
  int XePeriod;
  int nParams;
  std::vector <double> xyBinLower;
  std::vector <double> xyBinUpper;
  std::vector <double> xyBinCenter;
  //std::vector <double> yBinLower;
  //std::vector <double> yBinUpper;
  //std::vector <double> yBinCenter;

  TF1 *func;

  std::vector < std::vector < std::vector <double> > > triggMap;
 
};

  
