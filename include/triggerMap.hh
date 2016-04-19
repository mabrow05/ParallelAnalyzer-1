//
#include <vector>

class triggerMap {

public:

  triggerMap(int bin_width); //Constructor
  ~triggerMap();

  void readTriggerMap(int XeRunPeriod); //Read in a trigger map
  double returnTriggerProbability(double x, double y, double En);
  std::vector <double> getTriggerFunction(double x, double y); //returns the trigger function for the given event position
  void writeTriggerMap(int XeRunPeriod); //Write trigger map
  void setTriggerMapPoint(int xBin, int yBin, double val); //Set a trigger map val

  int getBinNumber(double pos); //returns bin number for position

  int getBinWidth() {return binWidth;}
  int getNbins() {return nBinsTotal;}
  double getBinUpper(int bin) {return xyBinUpper[bin];}
  double getBinLower(int bin) {return xyBinLower[bin];}
  double getBinCenter(int bin) {return xyBinCenter[bin];}
  int getCurrentXeRunPeriod() {return XePeriod;}

private:

  void setBinValues();

  int binWidth; //x and y width of position bins
  int nBinsXY, nBinsTotal;
  int XePeriod;
  int nParams;
  std::vector <double> xyBinLower;
  std::vector <double> xyBinUpper;
  std::vector <double> xyBinCenter;
  //std::vector <double> yBinLower;
  //std::vector <double> yBinUpper;
  //std::vector <double> yBinCenter;

  std::vector < std::vector < std::vector <double> > > triggMap;
 
};

  
