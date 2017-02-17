#ifndef MWPCPOSRESPONSE_HH
#define MWPCPOSRESPONSE_HH

#include <vector>



class MWPCCathodeHandler {

public:

  MWPCCathodeHandler();
  ~MWPCCathodeHandler() {};
  
  void setCathResponse(double *ex,double *ey,double *wx,double *wy);
  void setPedestalValues(double *ex,double *ey,double *wx,double *wy);
  void findAllPositions();
  
  
  
  double getPosEX() { return posEX; }
  double getPosEY() { return posEY; }
  double getPosWX() { return posWX; }
  double getPosWY() { return posWY; }

  int getnClippedEX() { return clippedEX.size(); }
  int getnClippedEY() { return clippedEY.size(); }
  int getnClippedWX() { return clippedWX.size(); }
  int getnClippedWY() { return clippedWY.size(); }
  




private:
  double *cathEX, *cathEY, *cathWX, *cathWY;
  double *pedEX, *pedEY, *pedWX, *pedWY;
  double pedSubtrEX[16],pedSubtrEY[16],pedSubtrWX[16],pedSubtrWY[16];
  double posEX, posEY, posWX, posWY; // Final position of the event
  //int nClippedEX, nClippedEY, nClippedWX, nClippedWY;
  std::vector <int> clippedEX, clippedEY, clippedWX, clippedWY; // Vector of clipped wires
  std::vector <int> signalEX, signalEY, signalWX, signalWY; // Holds wires with signal above threshold
  std::vector <int> fitWiresEX, fitWiresEY, fitWiresWX, fitWiresWY; // The wires to be fitted

  void doThresholdCheck(); // Sets all of the wires that were above threshold in signal**
  
  std::vector <int> chooseWires(double *cath, std::vector <int> clipped);
  std::vector <int>  doClipping(std::vector<int> wires, double *sig); //Checks for clipped wires
  double fitCathResponse(std::vector <int> wires, std::vector<int> clip, double *sig, double *pos); // Returns the position of the cathode signal in a plane
  int getMaxWire(std::vector <int> wires, double *sig); // returns max wire

  
  bool boolPedSubtr; // Whether or not the pedestal subtraction has been calculated
  bool boolThresholdCheck; // Whether or not the threshold has been checked has been calculated
  
  const double cathodeThreshold = 100.;
  const double wireposEX[16] = {76.20, 66.04, 55.88, 45.72, 35.56, 25.40, 15.24, 5.08, 
				-5.08, -15.24, -25.40, -35.56, -45.72, -55.88, -66.04, -76.20};
  const double wireposEY[16] = {76.20, 66.04, 55.88, 45.72, 35.56, 25.40, 15.24, 5.08, 
				-5.08, -15.24, -25.40, -35.56, -45.72, -55.88, -66.04, -76.20};
  const double wireposWY[16] = {76.20, 66.04, 55.88, 45.72, 35.56, 25.40, 15.24, 5.08, 
				-5.08, -15.24, -25.40, -35.56, -45.72, -55.88, -66.04, -76.20};
  const double wireposWX[16] = {-76.20, -66.04, -55.88, -45.72, -35.56, -25.40, -15.24, -5.08,
				5.08, 15.24, 25.40, 35.56, 45.72, 55.88, 66.04, 76.20};
}


#endif


  
