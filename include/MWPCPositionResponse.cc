#include "MWPCPositionResponse.hh" 



MWPCCathodeHandler::MWPCCathodeHandler() {

  posEX = posEY = posWX = posWY = 0.;
  nClippedEX = nClippedEY = nClippedWX = nClippedWY = 0;
  goodWiresEX = goodWiresEY = goodWiresWX = goodWiresWY = {0.};
  clippedEX.resize(16,0);
  clippedEY.resize(16,0);
  clippedWX.resize(16,0);
  clippedWY.resize(16,0);
  boolClipped = false;
  boolPedSubtr = false;

};

void MWPCCathodeHandler::setCathResponse(double *ex,double *ey,double *wx,double *wy) {
  cathEX = ex;
  cathEY = ey;
  cathWX = wx;
  cathWY = wy;
};
  
void MWPCCathodeHandler::setPedestalValues(double *ex,double *ey,double *wx,double *wy) {
  pedEX = ex;
  pedEY = ey;
  pedWX = wx;
  pedWY = wy;
};


void MWPCCathodeHandler::doClipping() {

  for ( int i=0; i<16; ++i ) {
    if ( cathEX[i]>4095. ) { clippedEX[i]=1; nClippedEX++; }
    if ( cathEY[i]>4095. ) { clippedEY[i]=1; nClippedEY++; }
    if ( cathWX[i]>4095. ) { clippedWX[i]=1; nClippedWX++; }
    if ( cathWY[i]>4095. ) { clippedWY[i]=1; nClippedWY++; }
 }
  boolClipped = true;
};

void MWPCCathodeHandler::doPedestalSubtraction() {

  for ( int i=0; i<16; i++ ) {   
    pedSubtrEX[i] = cathEX[i] - pedEX[i];
    pedSubtrEY[i] = cathEY[i] - pedEY[i];
    pedSubtrWX[i] = cathWX[i] - pedWX[i];
    pedSubtrWY[i] = cathWY[i] - pedWY[i];   
  }
  boolPedSubtr = true;
};



void MWPCCathodeHandler::fitCathResponse() {

  if ( !boolClipped )  doClipping();
  if ( !boolPedSubtr ) doPedestalSubtraction();
  
  if ( nClippedEX > 2 ) /////////////////////////////////////

  chooseWires(pedSubtrEX, clippedEX);
  chooseWires(pedSubtrEY, clippedEY);
  chooseWires(pedSubtrWX, clippedWX);
  chooseWires(pedSubtrWY, clippedWY);


  posEX = fitGaus(goodWiresEX);
  posEY = fitGaus(goodWiresEY);
  
  posWX = fitGaus(goodWiresWX);
  posWY = fitGaus(goodWiresWY);

};

std::vector <int> MWPCCathodeHandler::chooseWires(double *cath, std::vector <int> clipped) {
  
  int maxWire = 100;
  double maxVal = cathodeThreshold;
  int wiresAboveThresh = 0;

  std::vector <int> wires;

  // Choose max wire
  for ( int i=0; i<16; ++i ) {
    if ( cath[i] > cathodeThreshold ) wiresAboveThresh++;

    if ( cath[i] > maxVal && !clipped[i] ) {
      maxWire = i;
      maxVal = cath[i];
    }
  }

  
  
  if ( maxWire!=100 && maxWire!=0 && maxWire!=15 ) { // So there is a max wire not on the edge
                                                     // and it is not a clipped wire   
    int lowWire = -100;
    int highWire = 100;
    
    for ( int i = maxWire; i>-1; --i )
      if ( !clipped[i] ) { lowWire == i; break; }
    
    for ( int i = maxWire; i<16; ++i )
      if ( !clipped[i] ) { highWire == i; break; }
    
    if ( lowWire!=-100 ) wires.push_back(lowWire);
    wires.push_back(maxWire);
    if ( highWire!=100 ) wires.push_back(highWire);
    
    return wires; 
  }
  
  else if ( maxWire==0 ) {

    int lowWire = -1;
    int highWire = 100;

  }

    
    
      
};

void MWPCCathodeHandler::fitGaus(double *wires) {
  //make sure all three points are positive
};
