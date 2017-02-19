#ifndef MWPCPOSRESPONSE_HH
#define MWPCPOSRESPONSE_HH

#include <vector>



class MWPCCathodeHandler {

  
  
public:

  MWPCCathodeHandler(double *ex,double *ey,double *wx,double *wy,
		     double *pedex,double *pedey,double *pedwx,double *pedwy);
  ~MWPCCathodeHandler() {};
  void PrintSignals();
  void findAllPositions(); // You only need to call this to return all position.
  //void setCathResponse();
  //void setPedestalValues(double *ex,double *ey,double *wx,double *wy);
 
  void doPedestalSubtraction(); // Subtracts pedestals and fills pedSubtr** 
  void doThresholdCheck(); // Sets all of the wires that were above threshold in signal**
  
  
  std::vector<double> getPosEX() { return posEX; }
  std::vector<double> getPosEY() { return posEY; }
  std::vector<double> getPosWX() { return posWX; }
  std::vector<double> getPosWY() { return posWY; }

  int getnClippedEX() { return clippedEX.size(); }
  int getnClippedEY() { return clippedEY.size(); }
  int getnClippedWX() { return clippedWX.size(); }
  int getnClippedWY() { return clippedWY.size(); }

  int getMultEX() { return signalEX.size(); }
  int getMultEY() { return signalEY.size(); }
  int getMultWX() { return signalWX.size(); }
  int getMultWY() { return signalWY.size(); }
  
  int getMaxWireEX() { return getMaxWire(signalEX,pedSubtrEX); }
  int getMaxWireEY() { return getMaxWire(signalEX,pedSubtrEY); }
  int getMaxWireWX() { return getMaxWire(signalEX,pedSubtrWX); }
  int getMaxWireWY() { return getMaxWire(signalEX,pedSubtrWY); }
  




private:
  double *cathEX, *cathEY, *cathWX, *cathWY;
  double *pedEX, *pedEY, *pedWX, *pedWY;
  double pedSubtrEX[16],pedSubtrEY[16],pedSubtrWX[16],pedSubtrWY[16];
  std::vector<double> posEX, posEY, posWX, posWY; // Final position of the event (mean,width,height)
  //int nClippedEX, nClippedEY, nClippedWX, nClippedWY;
  std::vector <int> clippedEX, clippedEY, clippedWX, clippedWY; // Vector of clipped wires
  std::vector <int> signalEX, signalEY, signalWX, signalWY; // Holds wires with signal above threshold

  std::vector <int>  doClipping(std::vector<int> wires, double *sig); //Checks for clipped wires
  std::vector<double> fitCathResponse(std::vector <int> wires, std::vector<int> clip, double *sig, const double *pos); // // returns gaussian mean and width and height of a signal in a three element vector
  int getMaxWire(std::vector <int> wires, double *sig); // returns max wire
  int getMaxWireNotClipped(std::vector <int> wires,std::vector<int> clipped, double *sig); // returns max wire that was not clipped
  std::vector <int> getNonClippedSorted( const std::vector<int>& wires, const std::vector<int>& clipWires, double *sig);

  std::vector<double> fitGaus(std::vector<double> x, std::vector<double> y); // returns gaussian mean and width and height of data arrays in a 3 element vector

  
  bool boolPedSubtr{false}; // Whether or not the pedestal subtraction has been calculated
  bool boolThresholdCheck{false}; // Whether or not the threshold has been checked has been calculated
 
  const double cathodeThreshold {100.};
  const double wireposEX[16] {76.20, 66.04, 55.88, 45.72, 35.56, 25.40, 15.24, 5.08, 
      -5.08, -15.24, -25.40, -35.56, -45.72, -55.88, -66.04, -76.20};
  const double wireposEY[16] {76.20, 66.04, 55.88, 45.72, 35.56, 25.40, 15.24, 5.08, 
      -5.08, -15.24, -25.40, -35.56, -45.72, -55.88, -66.04, -76.20};
  const double wireposWY[16] {76.20, 66.04, 55.88, 45.72, 35.56, 25.40, 15.24, 5.08, 
      -5.08, -15.24, -25.40, -35.56, -45.72, -55.88, -66.04, -76.20};
  const double wireposWX[16] {-76.20, -66.04, -55.88, -45.72, -35.56, -25.40, -15.24, -5.08,
      5.08, 15.24, 25.40, 35.56, 45.72, 55.88, 66.04, 76.20};
  
};


#endif


  
