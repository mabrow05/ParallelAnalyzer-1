#include "MWPCPositionResponse.hh" 
#include <TF1.h>
#include <TGraphErrors.h>
#include <cmath>
#include <TMath.h>


MWPCCathodeHandler::MWPCCathodeHandler(double *ex,double *ey,double *wx,double *wy,
				       double *pedex,double *pedey,double *pedwx,double *pedwy) {

  cathodeThreshold = 100.; //Default for data
  clipThresholdEX = clipThresholdEY = clipThresholdWX = clipThresholdWY = 4090.; //Default for data

  cathEX = ex;
  cathEY = ey;
  cathWX = wx;
  cathWY = wy;

  pedEX = pedex;
  pedEY = pedey;
  pedWX = pedwx;
  pedWY = pedwy;

  _bGaus = false;
  _bGausAllEvents = false;
  _sigma = 6.8;

  doPedestalSubtraction();

};

MWPCCathodeHandler::MWPCCathodeHandler(double *ex,double *ey,double *wx,double *wy) {
  
  cathodeThreshold = 100.; //Default for data
  clipThresholdEX = clipThresholdEY = clipThresholdWX = clipThresholdWY = 4090.; //Default for data

  cathEX = ex;
  cathEY = ey;
  cathWX = wx;
  cathWY = wy;

  double pedHold[]{0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};

  pedEX = &pedHold[0];
  pedEY = &pedHold[0];
  pedWX = &pedHold[0];
  pedWY = &pedHold[0];
  /*pedEX = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
  pedEX = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
  pedEX = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
  pedEX = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};*/
  
  doPedestalSubtraction(); 

};

void MWPCCathodeHandler::PrintSignals() {
  
  for (int i=0; i<16; i++) {
    std::cout << pedSubtrEX[i] << "\t"
	      << pedSubtrEY[i] << "\t"
	      << pedSubtrWX[i] << "\t"
	      << pedSubtrWY[i] <<  std::endl;
  }
  std::cout << "//////////////////////////" << std::endl;
};


void MWPCCathodeHandler::doPedestalSubtraction() {

  for ( int i=0; i<16; i++ ) {   
    pedSubtrEX[i] = cathEX[i]<clipThresholdEX ? cathEX[i] - pedEX[i] : clipThresholdEX;
    pedSubtrEY[i] = cathEY[i]<clipThresholdEY ? cathEY[i] - pedEY[i] : clipThresholdEY;
    pedSubtrWX[i] = cathWX[i]<clipThresholdWX ? cathWX[i] - pedWX[i] : clipThresholdWX;
    pedSubtrWY[i] = cathWY[i]<clipThresholdWY ? cathWY[i] - pedWY[i] : clipThresholdWY;   
  }
  boolPedSubtr = true;
};


void MWPCCathodeHandler::doThresholdCheck() {
  
  for ( int i=0; i<16; i++ ) {   
    if ( pedSubtrEX[i] > cathodeThreshold ) signalEX.push_back(i);
    if ( pedSubtrEY[i] > cathodeThreshold ) signalEY.push_back(i);
    if ( pedSubtrWX[i] > cathodeThreshold ) signalWX.push_back(i);
    if ( pedSubtrWY[i] > cathodeThreshold ) signalWY.push_back(i);
  }  
  boolThresholdCheck = true;
};


std::vector<int> MWPCCathodeHandler::doClipping(std::vector<int> wires, 
						double* signal, double clipThresh) {

  std::vector <int> clipped;
  for ( unsigned int i=0; i<wires.size(); ++i ) {
    if ( signal[wires[i]] > clipThresh ) clipped.push_back(wires[i]); 
 }
  
  return clipped;
};





void MWPCCathodeHandler::findAllPositions(bool gaus, bool gausAllEvents) {

  _bGaus = gaus;
  _bGausAllEvents = gausAllEvents;
  
  if ( !boolPedSubtr ) doPedestalSubtraction();
  if ( !boolThresholdCheck ) doThresholdCheck();

  //East X
  clippedEX = doClipping(signalEX, cathEX, clipThresholdEX);
  //std::cout << clippedEX.size() << "\n";
  posEX = fitCathResponse(signalEX, clippedEX, pedSubtrEX, wireposEX);
  
  //East Y
  clippedEY = doClipping(signalEY, cathEY, clipThresholdEY);
  posEY = fitCathResponse(signalEY, clippedEY, pedSubtrEY, wireposEY);
  
  //West X
  clippedWX = doClipping(signalWX, cathWX, clipThresholdWX);
  posWX = fitCathResponse(signalWX, clippedWX, pedSubtrWX, wireposWX);
  
  //West Y
  clippedWY = doClipping(signalWY, cathWY, clipThresholdWY);
  posWY = fitCathResponse(signalWY, clippedWY, pedSubtrWY, wireposWY);

};

int MWPCCathodeHandler::getMaxWire(std::vector <int> wires, double *sig) {
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
  

std::vector<double> MWPCCathodeHandler::fitCathResponse(std::vector <int> wires, std::vector<int> clip, double *sig, const double *pos) {

  std::vector <double> finalPos(3,0.);

  std::vector <double> posFit(3,0.);
  std::vector <double> valFit(3,0.);
  double posInc = 10.16;

  unsigned int numWires = wires.size();
  unsigned int nClipped = clip.size();

  if ( numWires==0 ) return finalPos; // If there are no wires over threshold, then we just have the origin

  if ( ( nClipped == 0 && !(_bGaus && _bGausAllEvents) ) || !_bGaus ) { // This should be executed if _bGaus=false and if ( _bGaus=true
                                                                        // && _bGausAllEvents=false && nClipped==0 ) and it should be bypassed  
    int maxWire = getMaxWire(wires, sig);                               // if ( _bGaus=true && _bGausAllEvents=true )
    //std::cout << " Max Wire = " << maxWire << std::endl;
    double ave = 0.; double denom = 0.;
    for ( unsigned int i=0; i<numWires; ++i ) {
      ave += pos[wires[i]]*sig[wires[i]];
      denom += sig[wires[i]];
    }
    finalPos[0] = ave/denom;
    finalPos[1] = fabs( pos[wires[0]] - pos[wires[numWires-1]] )/2.;
    finalPos[2] = sig[maxWire];
    return finalPos;
  }
  
  // Here we handle any clipped wires... As per what everyone suggested
  
  std::vector <int> nonClipped = getNonClippedSequential(wires,clip,sig);
  int maxWireNC = getMaxWire(nonClipped, sig);
  //std::cout << "Max wire not clipped = " << maxWireNC << std::endl;
  //std::cout << "Num non clipped = " << nonClipped.size() << std::endl;

  //First determine which one is the max wire
    unsigned int maxIndex = 0;
    for ( auto w : nonClipped ) {
      if ( w < maxWireNC ) maxIndex++;
      //std::cout << w << " " << pos[w] << " " << sig[w] << std::endl;
    }
    //std::cout << "max index: " << maxIndex << std::endl;
    //exit(0);
  
  if ( nonClipped.size() >= 3 ) { // If we have 3 wires that are not clipped use the max wire and two others to determine gaussian
    
    // Check if max index is on the front edge of the vector
    if ( maxIndex == 0 ) {
      for ( unsigned int i=0; i<3; ++i ) { posFit[i] = pos[nonClipped[i]]; valFit[i] = sig[nonClipped[i]]; }
    }

    // Check if max index is on the back edge of the vector
    else if ( maxIndex == ( nonClipped.size()-1 ) ) {
      for ( unsigned int i= (nonClipped.size()-3); i<nonClipped.size(); ++i ) { posFit[i-nonClipped.size()+3] = pos[nonClipped[i]]; valFit[i-nonClipped.size()+3] = sig[nonClipped[i]]; }
    }

    else {
      for ( unsigned int i=(maxIndex-1); i<(maxIndex+2); ++i ) { posFit[i-maxIndex+1] = pos[nonClipped[i]]; valFit[i-maxIndex+1] = sig[nonClipped[i]]; }
    }
    
    //std::cout << posFit[0] << " " << posFit[1] << " " << posFit[2] << std::endl;
    //std::cout << valFit[0] << " " << valFit[1] << " " << valFit[2] << std::endl;

    //Check that the 3 wires being used are reasonably close to each other
    if ( fabs(posFit[0]-posFit[1]) < 5.1*posInc && fabs(posFit[2]-posFit[1]) < 5.1*posInc ) return fitGaus(posFit, valFit);
    else { // Just take weighted average because there must be wires triggering somewhere else..
      int maxWire = getMaxWire(wires, sig);
      //std::cout << " Max Wire = " << maxWire << std::endl;
      double ave = 0.; double denom = 0.;
      for ( unsigned int i=0; i<numWires; ++i ) {
	ave += pos[wires[i]]*sig[wires[i]];
	denom += sig[wires[i]];
      }
      finalPos[0] = ave/denom;
      finalPos[1] = fabs( pos[wires[0]] - pos[wires[numWires-1]] )/2.;
      finalPos[2] = sig[maxWire];
      return finalPos;
    }
  }

  else if ( nonClipped.size() == 2 ) { // If we have 2 wires that are not clipped use the special method fitGaus2wires

    std::vector <double> specPos { pos[ nonClipped[0] ], pos[ nonClipped[1] ] };
    std::vector <double> specVal { sig[ nonClipped[0] ], sig[ nonClipped[1] ] };
    return fitGaus2Points(specPos, specVal);
  }

  
  // Just take the weighted average for anything less than two non-clipped above threshold wires
  else {
    int maxWire = getMaxWire(wires, sig);
    //std::cout << " Max Wire = " << maxWire << std::endl;
    double ave = 0.; double denom = 0.;
    for ( unsigned int i=0; i<numWires; ++i ) {
      ave += pos[wires[i]]*sig[wires[i]];
      denom += sig[wires[i]];
    }
    finalPos[0] = ave/denom;
    finalPos[1] = fabs( pos[wires[0]] - pos[wires[numWires-1]] )/2.;
    finalPos[2] = sig[maxWire];
    return finalPos;
  }
};
  

 


std::vector <int> MWPCCathodeHandler::getNonClippedSorted( const std::vector<int>& wires, const std::vector<int>& clipWires, double *sig) {
  
  std::vector <int> sorted;
  std::vector <double> sortedSignals;

  bool clip = false;

  for ( unsigned int i=0; i<wires.size(); ++i ) {
    
    clip=false;
    
    for ( unsigned int j=0; j<clipWires.size(); ++j ) {
      if (wires[i]==clipWires[j]) { clip = true; continue; }
    }
    
    if ( !clip ) sorted.push_back(wires[i]), sortedSignals.push_back(sig[wires[i]]);
    
  }
  
  unsigned int veclen = sortedSignals.size();
  
  for ( unsigned int i=1; i<veclen; ++i ) {
    
    for ( unsigned int j=0; j<i; ++j ) {

      if ( sortedSignals[i]>sortedSignals[j] ) {
	sorted.insert(sorted.begin()+j,1,sorted[i]);
	sorted.erase(sorted.begin()+i+1);
	sortedSignals.insert(sortedSignals.begin()+j,1,sortedSignals[i]);
	sortedSignals.erase(sortedSignals.begin()+i+1);
	continue;
      }
    }
  }

  //for ( unsigned int j=0; j<clipWires.size(); ++j ) std::cout << "Clip: " << clipWires[j] << " " << std::endl;
  //for ( unsigned int j=0; j<wires.size(); ++j ) std::cout << "Wires: " << wires[j] << " " << std::endl;
  //for ( unsigned int j=0; j<sorted.size(); ++j ) std::cout << "No Clips Sorted: " << sorted[j] << " " << sortedSignals[j] << std::endl;

  //exit(0);
  return sorted;
  
};

std::vector <int> MWPCCathodeHandler::getNonClippedSequential( const std::vector<int>& wires, const std::vector<int>& clipWires, double *sig) {
  
  std::vector <int> ncwires;
  std::vector <double> ncsignals;

  bool clip = false;

  for ( unsigned int i=0; i<wires.size(); ++i ) {
    
    clip=false;
    
    for ( unsigned int j=0; j<clipWires.size(); ++j ) {
      if (wires[i]==clipWires[j]) { clip = true; continue; }
    }
    
    if ( !clip ) ncwires.push_back(wires[i]), ncsignals.push_back(sig[wires[i]]);
    
  }
  

  //for ( unsigned int j=0; j<clipWires.size(); ++j ) std::cout << "Clip: " << clipWires[j] << " " << std::endl;
  //for ( unsigned int j=0; j<wires.size(); ++j ) std::cout << "Wires: " << wires[j] << " " << std::endl;
  //for ( unsigned int j=0; j<sorted.size(); ++j ) std::cout << "No Clips Sorted: " << sorted[j] << " " << sortedSignals[j] << std::endl;

  //exit(0);
  return ncwires;
  
};
  
std::vector<double> MWPCCathodeHandler::fitGaus(std::vector<double> x, std::vector<double> y) {
 
  std::vector <double> params;

  
  // Figure out 
  double p_minus = x[0]-x[1];
  double p_plus  = x[2]-x[1];
  double s_minus = TMath::Log( y[0] );
  double s_0     = TMath::Log( y[1] );
  double s_plus  = TMath::Log( y[2] );

  //First calculate sigma
  double sigma2 = p_plus*p_minus*(p_minus - p_plus) / ( 2.*( p_plus*(s_0 - s_minus) - p_minus*(s_0 - s_plus) ) );

  // Calculate mean
  double mu = ( p_plus*p_plus*(s_0 - s_minus) - p_minus*p_minus*(s_0 - s_plus) ) / ( 2.* ( p_plus*(s_0 - s_minus) - p_minus*(s_0 - s_plus) ) ) ;

  // Calculate amplitude

  double amp = ( p_minus - p_plus )*p_plus*p_minus*s_0 / ( p_minus*p_plus*( p_minus - p_plus ) ) + mu*mu / ( 2.*sigma2 );

  params.push_back(mu + x[1]);
  params.push_back( sqrt( fabs(sigma2) ) );
  params.push_back(amp);

  return params;
  
};

std::vector<double> MWPCCathodeHandler::fitParabola(std::vector<double> x, std::vector<double> y) {
 
  std::vector <double> params;

  
  // Figure out 
  double p_minus = x[0]-x[1];
  double p_plus = x[2]-x[1];
  double s_minus = y[0];
  double s_0 = y[1];
  double s_plus = y[2];
  
  //First calculate sigma
  double sigma2 = p_plus*p_minus*(p_minus - p_plus) / ( 2.*( p_plus*(s_0 - s_minus) - p_minus*(s_0 - s_plus) ) );

  // Calculate mean
  double mu = ( p_plus*p_plus*(s_0 - s_minus) - p_minus*p_minus*(s_0 - s_plus) ) / ( 2.* ( p_plus*(s_0 - s_minus) - p_minus*(s_0 - s_plus) ) ) ;

  // Calculate amplitude

  double amp = ( p_minus - p_plus )*p_plus*p_minus*s_0 / ( p_minus*p_plus*( p_minus - p_plus ) ) + mu*mu / ( 2.*sigma2 );

  params.push_back(mu + x[1]);
  params.push_back( TMath::Sqrt( sigma2 ) );
  params.push_back(amp);

  return params;
  
};


std::vector<double> MWPCCathodeHandler::fitGaus2Points(std::vector<double> x, std::vector<double> y) {
 
  std::vector <double> params;

  
  // Figure out 
  double p_minus = x[0];
  double p_plus = x[1];
  double s_minus = y[0];
  double s_plus = y[1];
  
  //We will use the average width as seen from studying the normal fit to gaussian
  double sigma = _sigma / TMath::Sqrt( 0.6 );
  // Calculate mean
  double mu = (p_minus + p_plus) / 2. + sigma*sigma / (p_plus - p_minus) * TMath::Log(s_plus/s_minus) ;

  // Calculate amplitude

  double amp = ( p_plus*TMath::Log(s_minus) - p_minus*TMath::Log(s_plus) ) / ( p_plus - p_minus ) + ( mu*mu - p_plus*p_minus) / ( 2.*sigma*sigma );

  params.push_back( mu );
  params.push_back( sigma );
  params.push_back( amp );

  return params;
  
};
