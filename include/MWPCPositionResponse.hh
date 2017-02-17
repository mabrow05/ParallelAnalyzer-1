#ifndef MWPCPOSRESPONSE_HH
#define MWPCPOSRESPONSE_HH

#include <vector>



class MWPCCathodeHandler {

public:

  MWPCCathodeHandler();
  ~MWPCCathodeHandler() {};
  
  void setCathResponse(double *ex,double *ey,double *wx,double *wy);
  void setPedestalValues(double *ex,double *ey,double *wx,double *wy);
  void doClipping();
  void fitCathResponse();
  
  double getPosEX() { return posEX; }
  double getPosEY() { return posEY; }
  double getPosWX() { return posWX; }
  double getPosWY() { return posWY; }

  int getnClippedEX() { return nClippedEX; }
  int getnClippedEY() { return nClippedEY; }
  int getnClippedWX() { return nClippedWX; }
  int getnClippedWY() { return nClippedWY; }
  




private:
  double *cathEX, *cathEY, *cathWX, *cathWY;
  double *pedEX, *pedEY, *pedWX, *pedWY;
  double pedSubtrEX[16],pedSubtrEY[16],pedSubtrWX[16],pedSubtrWY[16];
  double posEX, posEY, posWX, posWY;
  int nClippedEX, nClippedEY, nClippedWX, nClippedWY;
  std::vector <int> clippedEX, clippedEY, clippedWX, clippedWY;

  std::vector <int> goodWiresEX,goodWiresEY,goodWiresWX,goodWiresWY;

  std::vector <int> chooseWires(double *cath, std::vector <int> clipped);

  bool boolClipped;  // Whether or not the clipping has been calculated
  bool boolPedSubtr; // Whether or not the pedestal subtraction has been calculated

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

