#include "MWPCPositionResponse.hh" 



MWPCCathodeHandler::MWPCCathodeHandler() {

  posEX = posEY = posWX = posWY = 0.;
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


void MWPCCathodeHandler::doPedestalSubtraction() {

  for ( int i=0; i<16; i++ ) {   
    pedSubtrEX[i] = cathEX[i] - pedEX[i];
    pedSubtrEY[i] = cathEY[i] - pedEY[i];
    pedSubtrWX[i] = cathWX[i] - pedWX[i];
    pedSubtrWY[i] = cathWY[i] - pedWY[i];   
  }
  boolPedSubtr = true;
};


std::vector<int> MWPCCathodeHandler::doClipping(std::vector<int> wires, double* signal) {

  std::vector <int> clipped;
  for ( unsigned int i=0; i<wires.size(); ++i ) {
    if ( signal[wires[i]] > 4090. ) clipped.push_back(wires[i]); 
  }
  return clipped;
};

void MWPCCathodeHandler::doThresholdCheck() {
  
  for ( int i=0; i<16; i++ ) {   
    if ( pedSubtrEX[i] > cathodeThreshold ) signalEX.push_back(i);
    if ( pedSubtrEY[i] > cathodeThreshold ) signalEY.push_back(i);
    if ( pedSubtrWX[i] > cathodeThreshold ) signalWX.push_back(i);
    if ( pedSubtrWY[i] > cathodeThreshold ) signalWY.push_back(i);
  }
  
};



void MWPCCathodeHandler::findAllPositions() {

  if ( !boolPedSubtr ) doPedestalSubtraction();
  if ( !boolThresholdCheck ) doThresholdCheck();

  //East X
  clippedEX = doClipping(signalEX, cathEX);
  posEX = fitCathResponse(signalEX, clippedEX, pedSubtrEX,wireposEX);
  
  //East Y
  clippedEY = doClipping(signalEY, cathEY);
  posEY = fitCathResponse(signalEY, clippedEY, pedSubtrEY, wireposEY);
  
  //West X
  clippedWX = doClipping(signalWX, cathWX);
  posWX = fitCathResponse(signalWX, clippedWX, pedSubtrWX,wireposWX);
  
  //West Y
  clippedWY = doClipping(signalWY, cathWY);
  posWY = fitCathResponse(signalWY, clippedWY, pedSubtrWY,wireposWY);

};

int getMaxWire(std::vector <int> wires, double *sig) {
  int maxWire = 0;
  double max = 0.;
  for ( unsigned int i=0; i<wires.size(); ++i ) {
    if ( sig[wires[i]]>max ) {
      max = sig[wires[i]];
      maxWire = wires[i];
    }
  }
  return maxWire;
};

int MWPCCathodeHandler::getMaxWireNotClipped(std::vector <int> wires, std::vector <int> clipped, double *sig) {
  
  int maxWire = 0;
  double max = 0.;
  // Choose max wire
  for ( unsigned int i=0; i<wires.size(); ++i ) {
    for ( unsigned int j=0; j<clipped.size(); ++j ) {
      if ( sig[wires[i]]>max && wires[i]!=clipped[j] ) {
	max = sig[wires[i]];
	maxWire = wires[i];
      }
    }
  }
  return maxWire;
};
  

double MWPCCathodeHandler::fitCathResponse(std::vector <int> wires, std::vector<int> clip, double *sig, double *pos) {

  std::vector <double> posFit(3,0.);
  std::vector <double> valFit(3,0.);
  double posInc = 10.16;

  unsigned int numWires = wires.size();
  unsigned int nClipped = clip.size();
  
  if ( numWires==0 ) return 0.;
  else if ( numWires==1 ) return pos[wires[0]];
  else if ( numWires==2 ) {
    
    if ( nClipped==2 ) return ( pos[wires[0]] + pos[wires[1]] ) / 2.;
    
    int max = getMaxWire(wires, sig);
    
    if ( max==0 ) { // Handle edge wires
      valFit[0] = cathodeThreshold; posFit[0] = pos[0]<0. ? (pos[0] - posInc) : (pos[0] + posInc) ;
      valFit[1] = sig[wires[0]]; posFit[1] = pos[wires[0]];
      valFit[2] = sig[wires[1]]; posFit[2] = pos[wires[1]];
    }
    else if ( max==15 ) {
      valFit[0] = cathodeThreshold; posFit[0] = pos[15]<0. ? (pos[15] - posInc) : (pos[15] + posInc) ;
      valFit[1] = sig[wires[0]]; posFit[1] = pos[wires[0]];
      valFit[2] = sig[wires[1]]; posFit[2] = pos[wires[1]];
    }
    else { // pick one more wire 

      if ( wires[0]==max ) { posFit[0] = pos[wires[0]-1]; valFit[0] = cathodeThreshold; }
      else { posFit[0] = pos[wires[1]+1]; valFit[0] = cathodeThreshold; }

      valFit[1] = sig[wires[0]]; posFit[1] = pos[wires[0]];
      valFit[2] = sig[wires[1]]; posFit[2] = pos[wires[1]];
    }
    return fitGaus(posFit, valFit);
  }

  else { // At least 3 wires above threshold

    int max = getMaxWire(wires, sig);
    
    if ( nClipped < 2 ) { // Less than 2 clipped, use max wire, and maxWire-1 and +1
      if ( max==0 ) {
	valFit[0] = sig[wires[0]]; posFit[0] = pos[wires[0]];
	valFit[1] = sig[wires[1]]; posFit[1] = pos[wires[1]];
	valFit[2] = sig[wires[2]]; posFit[2] = pos[wires[2]];
      }
      else if ( max==15 ) {
	valFit[0] = sig[wires[numWires-3]]; posFit[0] = pos[wires[numWires-3]];
	valFit[1] = sig[wires[numWires-2]]; posFit[1] = pos[wires[numWires-2]];
	valFit[2] = sig[wires[numWires-1]]; posFit[2] = pos[wires[numWires-1]];
      }
      else {
	valFit[0] = sig[max-1]; posFit[0] = pos[max-1];
	valFit[1] = sig[max]; posFit[1] = pos[max];
	valFit[2] = sig[max+1]; posFit[2] = pos[max+1];
      }
      return fitGaus(posFit, valFit);
    }

    else if ( nClipped==2 ) { // 2 clipped and at least 1 other signal 

      if ( numWires>3 ) { // 2 or more other wires over threshold

	if (clip[0]==0) {
	  valFit[0] = sig[wires[2]]; posFit[0] = pos[wires[2]];
	  valFit[1] = 4096.*1.25; posFit[1] = (pos[clip[0]] + pos[clip[1]])/2.;
	  valFit[2] = sig[wires[3]]; posFit[2] = pos[wires[3]];
	}
	else if (clip[1]==15) {
	  valFit[0] = sig[wires[numWires-4]]; posFit[0] = pos[wires[numWires-4]];
	  valFit[1] = 4096.*1.25; posFit[1] = (pos[clip[0]] + pos[clip[1]])/2.;
	  valFit[2] = sig[wires[numWires-3]]; posFit[2] = pos[wires[numWires-3]];
	}
	else {
	  valFit[0] = sig[clip[0]-1]; posFit[0] = pos[clip[0]-1];
	  valFit[1] = 4096.*1.25; posFit[1] = (pos[clip[0]] + pos[clip[1]])/2.;
	  valFit[2] = sig[clip[0]+2]; posFit[2] = pos[clip[0]+2];
	}
	return fitGaus(posFit, valFit);
      }

      else { // 1 other wire over threshold
	
	if ( clip[0]==wires[0] ) {
	  valFit[0] = sig[wires[2]]; posFit[0] = pos[wires[2]];
	  valFit[1] = 4096. + sig[wires[2]]; posFit[1] = pos[clip[1]];
	  valFit[2] = 4096.;  posFit[2] = pos[clip[0]];
	}
	else {
	  valFit[0] = sig[wires[0]]; posFit[0] = pos[wires[0]];
	  valFit[1] = 4096. + sig[wires[0]]; posFit[1] = pos[clip[0]];
	  valFit[2] = 4096.;  posFit[2] = pos[clip[1]];
	}
	return fitGaus(posFit, valFit);
      }
    }
      
    else {
      if ( (numWires - nClipped) > 2 ) { // Check that there are at least 3 points that aren't overflow. Use average of
                              // clipped*nClipped and the other wires above thresh
	double ave = 0.;
	for ( unsigned int i=0; i<nClipped; ++i ) {
	  ave += pos[clip[i]];
	}
	posFit[0] = ave/nClipped; valFit[0] = 4096.*nClipped;
	posFit[1] = pos[wires[0]]; valFit[1] = sig[wires[0]];
	posFit[2] = pos[wires[1]]; valFit[2] = sig[wires[1]];
	for ( unsigned int i=2; i<(numWires-nClipped); ++i ) {
	  posFit.push_back(pos[wires[i]]); valFit.push_back(sig[wires[i]]);
	}
	return fitGaus(posFit, valFit);
      }
    
    
      else { // Else you just take the average and call it a day
	double ave = 0.;
	for ( unsigned int i=0; i<nClipped; ++i ) {
	  ave += pos[clip[i]];
	}
	return ave/nClipped;
      }
    }
  }
};
  
  
void MWPCCathodeHandler::fitGaus(double *wires) {
  //make sure all three points are positive
};
