#include "MWPCPositionResponse.hh" 
#include <TF1.h>
#include <TGraphErrors.h>
#include <cmath>


MWPCCathodeHandler::MWPCCathodeHandler(double *ex,double *ey,double *wx,double *wy,
				       double *pedex,double *pedey,double *pedwx,double *pedwy) {

  cathEX = ex;
  cathEY = ey;
  cathWX = wx;
  cathWY = wy;

  pedEX = pedex;
  pedEY = pedey;
  pedWX = pedwx;
  pedWY = pedwy;

};

void MWPCCathodeHandler::PrintSignals() {
  
  for (int i=0; i<16; i++) {
    std::cout << pedSubtrEX[i] << std::endl;
  }
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


void MWPCCathodeHandler::doThresholdCheck() {
  
  for ( int i=0; i<16; i++ ) {   
    if ( pedSubtrEX[i] > cathodeThreshold ) signalEX.push_back(i);
    if ( pedSubtrEY[i] > cathodeThreshold ) signalEY.push_back(i);
    if ( pedSubtrWX[i] > cathodeThreshold ) signalWX.push_back(i);
    if ( pedSubtrWY[i] > cathodeThreshold ) signalWY.push_back(i);
  }  
  boolThresholdCheck = true;
};


std::vector<int> MWPCCathodeHandler::doClipping(std::vector<int> wires, double* signal) {

  std::vector <int> clipped;
  for ( unsigned int i=0; i<wires.size(); ++i ) {
    if ( signal[wires[i]] > 4090. ) clipped.push_back(wires[i]); 
  }
  return clipped;
};





void MWPCCathodeHandler::findAllPositions() {

  if ( !boolPedSubtr ) doPedestalSubtraction();
  if ( !boolThresholdCheck ) doThresholdCheck();

  //East X
  clippedEX = doClipping(signalEX, cathEX);
  //std::cout << clippedEX.size() << "\n";
  posEX = fitCathResponse(signalEX, clippedEX, pedSubtrEX, wireposEX);
  
  //East Y
  clippedEY = doClipping(signalEY, cathEY);
  posEY = fitCathResponse(signalEY, clippedEY, pedSubtrEY, wireposEY);
  
  //West X
  clippedWX = doClipping(signalWX, cathWX);
  posWX = fitCathResponse(signalWX, clippedWX, pedSubtrWX, wireposWX);
  
  //West Y
  clippedWY = doClipping(signalWY, cathWY);
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
  

std::vector<double> MWPCCathodeHandler::fitCathResponse(std::vector <int> wires, std::vector<int> clip, double *sig, const double *pos) {

  std::vector <double> finalPos(3,0.);

  std::vector <double> posFit(3,0.);
  std::vector <double> valFit(3,0.);
  double posInc = 10.16;

  unsigned int numWires = wires.size();
  unsigned int nClipped = clip.size();


  //////////// For time being, I'm just using brad's method until we straighten this out

   if ( numWires==0 ) return finalPos; // If there wasn't a signal...

   else {
     int max = getMaxWire(wires, sig);
     double ave = 0.; double denom = 0.;
     for ( unsigned int i=0; i<numWires; ++i ) {
       ave += pos[wires[i]]*sig[wires[i]];
       denom += sig[wires[i]];
     }
     finalPos[0] = ave/denom;
     finalPos[1] = fabs( pos[wires[0]] - pos[wires[numWires-1]] )/2.;
     finalPos[2] = sig[max];
     return finalPos;
   }
     

  // This is a sorted list of the non clipped, threshold passed wires
  std::vector<int> goodWires = getNonClippedSorted(wires,clip,sig);
  
  if ( goodWires.size()>2 ) { // If there are more than 2 wires which were good, fit them and go on.
    posFit[0] = pos[goodWires[1]]; valFit[0] = sig[goodWires[1]];
    posFit[1] = pos[goodWires[0]]; valFit[1] = sig[goodWires[0]];
    posFit[2] = pos[goodWires[2]]; valFit[2] = sig[goodWires[2]];
    return fitGaus(posFit, valFit);
  }

  if ( nClipped > 2 ) { // If there are more than 2 clipped, just take the weighted average
    double ave = 0.; double denom = 0.;
    for ( unsigned int i=0; i<numWires; ++i ) {
      ave += pos[wires[i]]*sig[wires[i]];
      denom += sig[wires[i]];
    }
    finalPos[0] = ave/denom;
    finalPos[1] = fabs( pos[clip[0]] - pos[clip[nClipped-1]] )/2.;
    finalPos[2] = nClipped*4096./2.;
    return finalPos;
  }

  if ( goodWires.size()==0 ) { // no good wires that didn't clip

    if ( numWires==0 ) return finalPos; // If there wasn't a signal...

    else { // Else you just take the average of the clipped wires and call it a day.. This won't happen often
      double ave = 0.; double denom = 0.;
      for ( unsigned int i=0; i<numWires; ++i ) {
	ave += pos[wires[i]]*sig[wires[i]];
	denom += sig[wires[i]];
      }
      finalPos[0] = ave/denom;
      finalPos[1] = fabs( pos[clip[0]] - pos[clip[nClipped-1]] )/2.;
      finalPos[2] = nClipped*4096.;
      return finalPos;
    }
  }

  else if ( goodWires.size()==1 ) { //One good wire
    
    if ( numWires==1 ) { // That one good wire is the only wire
      finalPos[0] = pos[goodWires[0]];
      finalPos[1] = 0.;
      finalPos[2] = sig[goodWires[0]];
      return finalPos;
    }
    
    else if ( numWires==2 ) { // Then we have one clipped and one not clipped. We pick another wire on the other side
      
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
    
    else if ( numWires==3 ) { // So one good wire, 2 clipped wires, anything above this was taken by the nClipped > 2 up above
      
      if ( pos[clip[0]] > pos[goodWires[0]] ) {
	valFit[0] = sig[goodWires[0]]; posFit[0] = pos[goodWires[0]];
	valFit[1] = 4096. + sig[goodWires[0]]; posFit[1] = pos[clip[0]];
	valFit[2] = 4096.;  posFit[2] = pos[clip[1]];
      }
      else {
	valFit[0] = sig[goodWires[0]]; posFit[0] = pos[goodWires[0]];
	valFit[1] = 4096. + sig[goodWires[0]]; posFit[1] = pos[clip[1]];
	valFit[2] = 4096.;  posFit[2] = pos[clip[0]];
      }
      return fitGaus(posFit, valFit);
    }
  }

  else if ( goodWires.size()==2 ) { // 2 good wires! Almost there

    if ( numWires==2 ) { // Only 2, no clipped, Then we pick another wire on the other side
      
      int max = getMaxWire(wires, sig);
      
      if ( max==0 ) { // Handle edge wires
	valFit[0] = cathodeThreshold; posFit[0] = pos[0]<0. ? (pos[0] - posInc) : (pos[0] + posInc) ;
	valFit[1] = sig[goodWires[0]]; posFit[1] = pos[goodWires[0]];
	valFit[2] = sig[goodWires[1]]; posFit[2] = pos[goodWires[1]];
      }
      else if ( max==15 ) {
	valFit[0] = cathodeThreshold; posFit[0] = pos[15]<0. ? (pos[15] - posInc) : (pos[15] + posInc) ;
	valFit[1] = sig[goodWires[0]]; posFit[1] = pos[goodWires[0]];
	valFit[2] = sig[goodWires[1]]; posFit[2] = pos[goodWires[1]];
      }
      else { // pick one more wire 
	
	if ( goodWires[0]==max ) { posFit[0] = pos[goodWires[0]-1]; valFit[0] = cathodeThreshold; }
	else { posFit[0] = pos[goodWires[1]+1]; valFit[0] = cathodeThreshold; }
	
	valFit[1] = sig[goodWires[0]]; posFit[1] = pos[goodWires[0]];
	valFit[2] = sig[goodWires[1]]; posFit[2] = pos[goodWires[1]];
      }
      return fitGaus(posFit, valFit);
    }
    
    else if ( numWires==3 ) { // 2 good, 1 clipped, Then we add the smaller wire to the clipped wire and fit
      
      valFit[0] = sig[goodWires[1]]; posFit[0] = pos[goodWires[1]]; // this is the smaller signal
      valFit[1] = 4096. + sig[goodWires[1]]; posFit[1] = pos[clip[0]]; // this is the clipped one plus the smaller wire
      valFit[2] = sig[goodWires[0]]; posFit[2] = pos[goodWires[0]]; // this is the larger signal
      return fitGaus(posFit, valFit);
    }
    
    else if ( numWires==4 ) { // 2 good, 2 clipped, then we take the center of the clipped wires as their position, and take 2
                              // 1.66 * 4096 as the signal
      
      valFit[0] = sig[goodWires[1]]; posFit[0] = pos[goodWires[1]]; // this is the smaller signal
      valFit[1] = 4096.*1.66 ; posFit[1] = ( pos[clip[0]] + pos[clip[1]] ) / 2.; // these are the clipped wires
      valFit[2] = sig[goodWires[0]]; posFit[2] = pos[goodWires[0]]; // this is the larger signal
      return fitGaus(posFit, valFit);
    }   

  }
  else {
    std::cout << "There is a hole in your logic\n";
    exit(0);
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
  
std::vector<double> MWPCCathodeHandler::fitGaus(std::vector<double> x, std::vector<double> y) {
 
  std::vector <double> params;
  //First making sure all points in y are above threshold. If they aren't for some odd reason,
  // I'll set them equal to threshold..
  std::vector <double> yerr;
  for ( auto &i : y ) {
    if ( i<cathodeThreshold ) i=cathodeThreshold;
    yerr.push_back(sqrt(i));
  }
  //for ( auto &i : yerr) std::cout << yerr << std::endl;
  //std::cout << x.size() << std::endl;
  //std::cout << y[0] << "\t" << y[1] << "\t" << y[2] <<"\n";
  TGraphErrors *g = new TGraphErrors(x.size(), &x[0], &y[0], 0, &yerr[0]);
  TF1 *func = new TF1("func","[0]*TMath::Exp(-(x-[1])*(x-[1])/(2.*[2]*[2]))", -90.,90.);
  //TF1 *func = new TF1("func","gaus(0)", -90.,90.);
  func->SetParameters(y[1], x[1], 10.);
  g->Fit(func,"RQ");

  if ( fabs(func->GetParameter(1)) < 87. && func->GetParameter(0)>0. && func->GetParameter(0)<1.e10) {
    params.push_back(func->GetParameter(1)); //mean
    params.push_back(fabs(func->GetParameter(2))); //width
    params.push_back(func->GetParameter(0)); //height
    return params;
  }

  else  {
    std::cout << "////////////////////////\n";
    for (UInt_t i=0; i<x.size(); ++i ) {
      std::cout << x[i] << "\t" << y[i] << "\t" << yerr[i] << "\t" << func->GetParameter(0) << "\t" << func->GetParameter(1) << "\n";
    }
    std::cout << "////////////////////////\n";
    
    func->SetParameter(0,1000.);
    func->SetParameter(1, x[1]);
    func->SetParameter(2, 20.);
    func->SetParLimits(1,-87., 87.);
    func->SetParLimits(0,0.,100000.);
    g->Fit(func,"RQ");

    if ( fabs(func->GetParameter(2)) < 100. ) {

      params.push_back(func->GetParameter(1)); //mean
      params.push_back(fabs(func->GetParameter(2))); //width
      params.push_back(func->GetParameter(0)); //height
    }
    else { 
      params.push_back(0.);
      params.push_back(0.);
      params.push_back(0.);
    }

    delete func; delete g;
    return params;
  }
  
};
