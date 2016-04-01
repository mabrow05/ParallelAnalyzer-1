/*
Defines different types of Asymmetries to be calculated within an Octet.
ppair, quartet, octet
*/

#include "Asymmetries.hh"
#include "MBUtils.hh"

#include <fstream>
#include <cstdlib>
#include <cmath>

AsymmetryBase::AsymmetryBase(int oct, double enBinWidth, double fidCut, bool ukdata, bool simulation) : UKdata(ukdata), Simulation(simulation), octet(oct), energyBinWidth(enBinWidth), fiducialCut(fidCut), boolAnaChRtVecs(false), runsInOctet(0), analysisChoice(1) {
  unsigned int numBins = (unsigned int)(1200./energyBinWidth);
  A2.resize(4,std::vector < std::vector <double> > (2,std::vector <double> (numBins,0.)));
  B2.resize(4,std::vector < std::vector <double> > (2,std::vector <double> (numBins,0.)));
  A5.resize(4,std::vector < std::vector <double> > (2,std::vector <double> (numBins,0.)));
  B5.resize(4,std::vector < std::vector <double> > (2,std::vector <double> (numBins,0.)));
  A7.resize(4,std::vector < std::vector <double> > (2,std::vector <double> (numBins,0.)));
  B7.resize(4,std::vector < std::vector <double> > (2,std::vector <double> (numBins,0.)));
  A10.resize(4,std::vector < std::vector <double> > (2,std::vector <double> (numBins,0.)));
  B10.resize(4,std::vector < std::vector <double> > (2,std::vector <double> (numBins,0.)));
  A2err.resize(4,std::vector < std::vector <double> > (2,std::vector <double> (numBins,0.)));
  B2err.resize(4,std::vector < std::vector <double> > (2,std::vector <double> (numBins,0.)));
  A5err.resize(4,std::vector < std::vector <double> > (2,std::vector <double> (numBins,0.)));
  B5err.resize(4,std::vector < std::vector <double> > (2,std::vector <double> (numBins,0.)));
  A7err.resize(4,std::vector < std::vector <double> > (2,std::vector <double> (numBins,0.)));
  B7err.resize(4,std::vector < std::vector <double> > (2,std::vector <double> (numBins,0.)));
  A10err.resize(4,std::vector < std::vector <double> > (2,std::vector <double> (numBins,0.)));
  B10err.resize(4,std::vector < std::vector <double> > (2,std::vector <double> (numBins,0.)));
  anaChoice_A2.resize(2,std::vector <double> (numBins,0.));
  anaChoice_B2.resize(2,std::vector <double> (numBins,0.));
  anaChoice_A5.resize(2,std::vector <double> (numBins,0.));
  anaChoice_B5.resize(2,std::vector <double> (numBins,0.));
  anaChoice_A7.resize(2,std::vector <double> (numBins,0.));
  anaChoice_B7.resize(2,std::vector <double> (numBins,0.));
  anaChoice_A10.resize(2,std::vector <double> (numBins,0.));
  anaChoice_B10.resize(2,std::vector <double> (numBins,0.));
  anaChoice_A2err.resize(2,std::vector <double> (numBins,0.));
  anaChoice_B2err.resize(2,std::vector <double> (numBins,0.));
  anaChoice_A5err.resize(2,std::vector <double> (numBins,0.));
  anaChoice_B5err.resize(2,std::vector <double> (numBins,0.));
  anaChoice_A7err.resize(2,std::vector <double> (numBins,0.));
  anaChoice_B7err.resize(2,std::vector <double> (numBins,0.));
  anaChoice_A10err.resize(2,std::vector <double> (numBins,0.));
  anaChoice_B10err.resize(2,std::vector <double> (numBins,0.));
  numEvtsEastByTypeByBin.resize(4,std::vector<double>(numBins,0.));
  numEvtsWestByTypeByBin.resize(4,std::vector<double>(numBins,0.));
  A2len.resize(2,std::vector <double> (2,0.));
  A5len.resize(2,std::vector <double> (2,0.));
  A7len.resize(2,std::vector <double> (2,0.));
  A10len.resize(2,std::vector <double> (2,0.));
  B2len.resize(2,std::vector <double> (2,0.));
  B5len.resize(2,std::vector <double> (2,0.));
  B7len.resize(2,std::vector <double> (2,0.));
  B10len.resize(2,std::vector <double> (2,0.));

  binLowerEdge.resize(numBins,0.);
  binUpperEdge.resize(numBins,0.);
  
  for (unsigned int i=0; i<numBins; i++) {
    binLowerEdge[i] = (double)(i)*energyBinWidth;
    binUpperEdge[i] = (double)(i+1)*energyBinWidth;
  }  

  readOctetFile();
 
};

void AsymmetryBase::readOctetFile() {
  char fileName[200];
  sprintf(fileName,"%s/All_Octets/octet_list_%i.dat",getenv("OCTET_LIST"),octet);
  std::ifstream infile(fileName);
  std::string runTypeHold;
  int runNumberHold;
  int numRuns = 0;
  // Populate the map with the runType and runNumber from this octet
  while (infile >> runTypeHold >> runNumberHold) {
    runType[runTypeHold] = runNumberHold;
    numRuns++;
  }
  infile.close();
  runsInOctet = numRuns;
  std::cout << "Read in octet file for octet " << octet << " with " << runsInOctet << " runs\n";
};

bool AsymmetryBase::isFullOctet() {
  bool A2=false, B2=false, A5=false, B5=false, A7=false, B7=false, A10=false, B10=false;
  std::map <std::string,int>::iterator it;
  for ( it = runType.begin();it!=runType.end(); it++) {
    if (it->first=="A2") {A2=true; continue;}
    if (it->first=="B2") {B2=true; continue;}
    if (it->first=="A5") {A5=true; continue;}
    if (it->first=="B5") {B5=true; continue;}
    if (it->first=="A7") {A7=true; continue;}
    if (it->first=="B7") {B7=true; continue;}
    if (it->first=="A10") {A10=true; continue;}
    if (it->first=="B10") {B10=true; continue;}
  }
  if (A2 && B2 && A5 && B5 && A7 && B7 && A10 && B10) return true;
  else return false;
};

void AsymmetryBase::loadRates() {
  std::map<std::string,int>::iterator it = runType.begin();
  BGSubtractedRate *bg;
  while (it!=runType.end()) {
    if (checkIfBetaRun(it->first)) {
      bg = new BGSubtractedRate(it->second,energyBinWidth,fiducialCut,UKdata,Simulation);

      std::cout << "initialized BGStubtractedRate for run " << it->second << " (BG run " 
		<< bg->getBackgroundRun(it->second) << ") \n";

      bg->calcBGSubtRates();

      for (int type=0;type<4;type++) {
	for (int side=0; side<2; side++) {  
      
	  if (it->first=="A2") {
	    A2[type][side] = bg->returnBGSubtRate(side,type);
	    A2err[type][side] = bg->returnBGSubtRateError(side,type);
	    A2len[0] = bg->returnRunLengths(true);
	    A2len[1] = bg->returnRunLengths(false);
	  }
	  else if (it->first=="B2") {
	    B2[type][side] = bg->returnBGSubtRate(side,type);
	    B2err[type][side] = bg->returnBGSubtRateError(side,type);
	    B2len[0] = bg->returnRunLengths(true);
	    B2len[1] = bg->returnRunLengths(false);
	  }
	  else if (it->first=="A5") {
	    A5[type][side] = bg->returnBGSubtRate(side,type);
	    A5err[type][side] = bg->returnBGSubtRateError(side,type);
	    A5len[0] = bg->returnRunLengths(true);
	    A5len[1] = bg->returnRunLengths(false);
	  }
	  else if (it->first=="B5") {
	    B5[type][side] = bg->returnBGSubtRate(side,type);
	    B5err[type][side] = bg->returnBGSubtRateError(side,type);
	    B5len[0] = bg->returnRunLengths(true);
	    B5len[1] = bg->returnRunLengths(false);
	  }
	  else if (it->first=="A7") {
	    A7[type][side] = bg->returnBGSubtRate(side,type);
	    A7err[type][side] = bg->returnBGSubtRateError(side,type);
	    A7len[0] = bg->returnRunLengths(true);
	    A7len[1] = bg->returnRunLengths(false);
	  }
	  else if (it->first=="B7") {
	    B7[type][side] = bg->returnBGSubtRate(side,type);
	    B7err[type][side] = bg->returnBGSubtRateError(side,type);
	    B7len[0] = bg->returnRunLengths(true);
	    B7len[1] = bg->returnRunLengths(false);
	  }
	  else if (it->first=="A10") {
	    A10[type][side] = bg->returnBGSubtRate(side,type);
	    A10err[type][side] = bg->returnBGSubtRateError(side,type);
	    A10len[0] = bg->returnRunLengths(true);
	    A10len[1] = bg->returnRunLengths(false);
	  }
	  else if (it->first=="B10") {
	    B10[type][side] = bg->returnBGSubtRate(side,type);
	    B10err[type][side] = bg->returnBGSubtRateError(side,type);
	    B10len[0] = bg->returnRunLengths(true);
	    B10len[1] = bg->returnRunLengths(false);
	  }
	  else throw "Run misidentified in loadRates";
	}
      }
	  
      delete bg;
    }
    it++;
  }
};

bool AsymmetryBase::checkIfBetaRun(std::string type) {
  if (type=="A2" || type=="A5" || type=="A7" || type=="A10" || type=="B2" || type=="B5" || type=="B7" || type=="B10") return true;
  else return false;
};

void AsymmetryBase::calcBGsubtractedEvts() {
  std::map<std::string,int>::iterator it = runType.begin();
  BGSubtractedRate *bg;
  while (it!=runType.end()) {
    if (checkIfBetaRun(it->first)) {
      bg = new BGSubtractedRate(it->second,energyBinWidth,fiducialCut,UKdata,Simulation);

      std::cout << "initialized BGStubtractedRate for run " << it->second << " (BG run " 
		<< bg->getBackgroundRun(it->second) << ") \n";

      bg->calcBGSubtRates();
      
      if (!Simulation) {
	std::vector<double> runLengthBeta = bg->returnRunLengths(true);
	std::vector<double> runLengthBG = bg->returnRunLengths(false);
	std::cout << "RunLength: \tE\tW\n\t\t" 
		  << runLengthBeta[0] << "\t" << runLengthBeta[1] << std::endl
		  << "\t\t" << runLengthBG[0] << "\t" << runLengthBG[1] << std::endl;
	
	for (int type=0;type<4;type++) {
	  std::vector <double> evecbg = bg->returnBGSubtRate(0,type);
	  std::vector <double> wvecbg = bg->returnBGSubtRate(1,type);
	  for (unsigned int bin=0; bin<evecbg.size(); bin++) {
	    numEvtsEastByTypeByBin[type][bin]+=runLengthBeta[0]*evecbg[bin];
	    numEvtsWestByTypeByBin[type][bin]+=runLengthBeta[1]*wvecbg[bin];
	  }
	}
      }
      else {
	for (int type=0;type<4;type++) {
	  std::vector <double> evecbg = bg->returnBGSubtRate(0,type);
	  std::vector <double> wvecbg = bg->returnBGSubtRate(1,type);
	  for (unsigned int bin=0; bin<evecbg.size(); bin++) {
	    numEvtsEastByTypeByBin[type][bin]+=evecbg[bin];
	    numEvtsWestByTypeByBin[type][bin]+=wvecbg[bin];
	  }
	}
      }
      delete bg;
    }
    it++;
  }
};

std::vector <double> AsymmetryBase::getNumBGsubtrEvts(double enWinLow, double enWinHigh, int evtType) {
  unsigned int binLow = (unsigned int)(enWinLow/energyBinWidth);
  unsigned int binHigh = (unsigned int)(enWinHigh/energyBinWidth)-1;

  std::vector <double> numEvts(2,0.);

  std::cout << "Number of Type " << evtType << " events for energy window " << binLowerEdge[binLow] << " - " << binUpperEdge[binHigh] << ": ";
 
  for (unsigned int i=binLow; i<=binHigh; i++) {
    numEvts[0]+=numEvtsEastByTypeByBin[evtType][i];
    numEvts[1]+=numEvtsWestByTypeByBin[evtType][i];
  }
  std::cout << "E->" << numEvts[0] << "  W->" << numEvts[1] << std::endl;
  return numEvts;
};

void AsymmetryBase::makeAnalysisChoiceRateVectors(int anaChoice) {
  if (isAnaChoiceRateVectors()) {
    unsigned int numBins = (unsigned int)(1200./energyBinWidth);
    anaChoice_A2.resize(2,std::vector <double> (numBins,0.));
    anaChoice_B2.resize(2,std::vector <double> (numBins,0.));
    anaChoice_A5.resize(2,std::vector <double> (numBins,0.));
    anaChoice_B5.resize(2,std::vector <double> (numBins,0.));
    anaChoice_A7.resize(2,std::vector <double> (numBins,0.));
    anaChoice_B7.resize(2,std::vector <double> (numBins,0.));
    anaChoice_A10.resize(2,std::vector <double> (numBins,0.));
    anaChoice_B10.resize(2,std::vector <double> (numBins,0.));
    anaChoice_A2err.resize(2,std::vector <double> (numBins,0.));
    anaChoice_B2err.resize(2,std::vector <double> (numBins,0.));
    anaChoice_A5err.resize(2,std::vector <double> (numBins,0.));
    anaChoice_B5err.resize(2,std::vector <double> (numBins,0.));
    anaChoice_A7err.resize(2,std::vector <double> (numBins,0.));
    anaChoice_B7err.resize(2,std::vector <double> (numBins,0.));
    anaChoice_A10err.resize(2,std::vector <double> (numBins,0.));
    anaChoice_B10err.resize(2,std::vector <double> (numBins,0.));
  }
  analysisChoice = anaChoice;
  unsigned int numBins = (unsigned int)(1200./energyBinWidth);
  unsigned int type_low=0, type_high=0;
  if (anaChoice==1 || anaChoice==3 || anaChoice==5) { type_low=0; type_high=3;}
  else if (anaChoice==2) { type_low=0; type_high=1;}
  else if (anaChoice==4) { type_low=0; type_high=0;}
  else if (anaChoice==6) { type_low=1; type_high=1;}
  else if (anaChoice==7 || anaChoice==8 || anaChoice==9) { type_low=2; type_high=3;}
  else throw "Bad Analysis Choice";
  
  for (unsigned int side=0; side<2; side++) {
    for (unsigned int bin=0; bin<numBins; bin++) {
      for (unsigned int type=type_low; type<=type_high; type++) {	 
	anaChoice_A2[side][bin]+=A2[type][side][bin];
	anaChoice_A5[side][bin]+=A5[type][side][bin];
	anaChoice_A7[side][bin]+=A7[type][side][bin];
	anaChoice_A10[side][bin]+=A10[type][side][bin];
	anaChoice_B2[side][bin]+=B2[type][side][bin];
	anaChoice_B5[side][bin]+=B5[type][side][bin];
	anaChoice_B7[side][bin]+=B7[type][side][bin];
	anaChoice_B10[side][bin]+=B10[type][side][bin];
	
	anaChoice_A2err[side][bin]+=power(A2err[type][side][bin],2);
	anaChoice_A5err[side][bin]+=power(A5err[type][side][bin],2);
	anaChoice_A7err[side][bin]+=power(A7err[type][side][bin],2);
	anaChoice_A10err[side][bin]+=power(A10err[type][side][bin],2);
	anaChoice_B2err[side][bin]+=power(B2err[type][side][bin],2);
	anaChoice_B5err[side][bin]+=power(B5err[type][side][bin],2);
	anaChoice_B7err[side][bin]+=power(B7err[type][side][bin],2);
	anaChoice_B10err[side][bin]+=power(B10err[type][side][bin],2);

	//std::cout << anaChoice_A2[side][bin] << " " << anaChoice_A2err[side][bin] << std::endl;
      }
    
      anaChoice_A2err[side][bin]=sqrt(anaChoice_A2err[side][bin]);
      anaChoice_A5err[side][bin]=sqrt(anaChoice_A5err[side][bin]);
      anaChoice_A7err[side][bin]=sqrt(anaChoice_A7err[side][bin]);
      anaChoice_A10err[side][bin]=sqrt(anaChoice_A10err[side][bin]);
      anaChoice_B2err[side][bin]=sqrt(anaChoice_B2err[side][bin]);
      anaChoice_B5err[side][bin]=sqrt(anaChoice_B5err[side][bin]);
      anaChoice_B7err[side][bin]=sqrt(anaChoice_B7err[side][bin]);
      anaChoice_B10err[side][bin]=sqrt(anaChoice_B10err[side][bin]);
    }
  }
  boolAnaChRtVecs = true;
  std::cout << "Constructed rate vectors for analysis choice " << anaChoice << ".\n";	  
};


std::vector < std::vector < std::vector<double> > > AsymmetryBase::returnBGsubtractedRate(std::string runType) {
  if (runType=="A2") return A2;
  else if (runType=="A5") return A5;
  else if (runType=="A7") return A7;
  else if (runType=="A10") return A10;
  else if (runType=="B2") return B2;
  else if (runType=="B5") return B5;
  else if (runType=="B7") return B7;
  else if (runType=="B10") return B10;
  else throw "Bad Run Type given to returnBGsubtractedRate(string runType)";

};

std::vector < std::vector < std::vector<double> > > AsymmetryBase::returnBGsubtractedRateError(std::string runType) {
  if (runType=="A2") return A2err;
  else if (runType=="A5") return A5err;
  else if (runType=="A7") return A7err;
  else if (runType=="A10") return A10err;
  else if (runType=="B2") return B2err;
  else if (runType=="B5") return B5err;
  else if (runType=="B7") return B7err;
  else if (runType=="B10") return B10err;
  else throw "Bad Run Type given to returnBGsubtractedRateError(string runType)";

};

///////////////////////////////////////////////////////////////////////////////////////////

OctetAsymmetry::OctetAsymmetry(int oct, double enBinWidth, double fidCut, bool ukdata, bool simulation) : AsymmetryBase(oct,enBinWidth,fidCut,ukdata,simulation), totalAsymmetry(0.), totalAsymmetryError(0.) {
  if (isFullOctet()) {
    unsigned int numBins = (unsigned int)(1200./energyBinWidth);
    asymmetry.resize(numBins,0.);
    asymmetryError.resize(numBins,0.);
    superSum.resize(numBins,0.);
    superSumError.resize(numBins,0.);
    loadRates(); // load the rates in the rate vectors for each run
    std::cout <<"//////////////////////////////////////////////////////////////////\n"
	      <<"Initialized OctetAsymmetry for octet " << octet << std::endl;
  }
  else throw "OCTET NOT A COMPLETE OCTET";
};

void OctetAsymmetry::calcAsymmetryBinByBin(int anaChoice) {
  if (!isAnaChoiceRateVectors() || getCurrentAnaChoice()!=anaChoice) {
    makeAnalysisChoiceRateVectors(anaChoice);
  }

  double sfON[2]={0.};
  double sfOFF[2]={0.};
  double sfON_err[2]={0.};
  double sfOFF_err[2]={0.};

  for (unsigned int bin=0; bin<asymmetry.size(); bin++) {
    double R = 0.;
    double deltaR = 0.;
    for (unsigned int side=0; side<2; side++) {
      double weightsum=0.;

      sfOFF[side] = (anaChoice_A2[side][bin]>0.?power(1./anaChoice_A2err[side][bin],2)*anaChoice_A2[side][bin]:0.
		     + anaChoice_A10[side][bin]>0.?power(1./anaChoice_A10err[side][bin],2)*anaChoice_A10[side][bin]:0.
		     + anaChoice_B5[side][bin]>0.?power(1./anaChoice_B5err[side][bin],2)*anaChoice_B5[side][bin]:0.
		     + anaChoice_B7[side][bin]>0.?power(1./anaChoice_B7err[side][bin],2)*anaChoice_B7[side][bin]:0.);
      weightsum = (anaChoice_A2[side][bin]>0.?power(1./anaChoice_A2err[side][bin],2):0.
		   + anaChoice_A10[side][bin]>0.?power(1./anaChoice_A10err[side][bin],2):0.
		   + anaChoice_B5[side][bin]>0.?power(1./anaChoice_B5err[side][bin],2):0.
		   + anaChoice_B7[side][bin]>0.?power(1./anaChoice_B7err[side][bin],2):0.);

      sfOFF[side] = weightSum>0. ? sfOFF[side]/weightsum : 0.;
      sfOFF_err[side] = weightSum>0. ? 1./sqrt(weightsum) : 0.;

      weightsum=0.;
      sfON[side] = (anaChoice_A5[side][bin]>0.?power(1./anaChoice_A5err[side][bin],2)*anaChoice_A5[side][bin]:0.
		     + anaChoice_A7[side][bin]>0.?power(1./anaChoice_A7err[side][bin],2)*anaChoice_A7[side][bin]:0.
		     + anaChoice_B2[side][bin]>0.?power(1./anaChoice_B2err[side][bin],2)*anaChoice_B2[side][bin]:0.
		     + anaChoice_B10[side][bin]>0.?power(1./anaChoice_B10err[side][bin],2)*anaChoice_B10[side][bin]:0.);
      weightsum = (anaChoice_A5[side][bin]>0.?power(1./anaChoice_A5err[side][bin],2):0.
		   + anaChoice_A7[side][bin]>0.?power(1./anaChoice_A7err[side][bin],2):0.
		   + anaChoice_B2[side][bin]>0.?power(1./anaChoice_B2err[side][bin],2):0.
		   + anaChoice_B10[side][bin]>0.?power(1./anaChoice_B10err[side][bin],2):0.);
      
      //sfON[side] = (power(1./anaChoice_A5err[side][bin],2)*anaChoice_A5[side][bin]+power(1./anaChoice_A7err[side][bin],2)*anaChoice_A7[side][bin]
      //	    +power(1./anaChoice_B2err[side][bin],2)*anaChoice_B2[side][bin]+power(1./anaChoice_B10err[side][bin],2)*anaChoice_B10[side][bin]);
      //weightsum = power(1./anaChoice_A5err[side][bin],2)+power(1./anaChoice_A7err[side][bin],2)
      //+power(1./anaChoice_B2err[side][bin],2)+power(1./anaChoice_B10err[side][bin],2);

      sfON[side] = weightSum>0. ? sfON[side]/weightsum : 0.;
      sfON_err[side] = weightSum>0. ? 1./sqrt(weightsum) : 0.;
      
      //if (bin==73 || bin==74) std::cout << sfOFF[side] << " " << sfON[side] << std::endl;
      //if (side==1) {
      //std::cout << anaChoice_A10[side][bin] << " " << anaChoice_A10err[side][bin] << std::endl;
      //}
    }
    if (sfOFF[0]>0. && sfOFF[1]>0. && sfON[0]>0. && sfON[1]>0.) { 
      R = sfOFF[0]*sfON[1]/(sfON[0]*sfOFF[1]);
      deltaR = sqrt(R*R*(power(sfOFF_err[0]/sfOFF[0],2)+power(sfON_err[1]/sfON[1],2)+power(sfOFF_err[1]/sfOFF[1],2)+power(sfON_err[0]/sfON[0],2)));
      asymmetry[bin] = (1.-sqrt(R))/(1+sqrt(R));
      asymmetryError[bin] = (deltaR)/(sqrt(std::abs(R))*power((sqrt(std::abs(R))+1.),2));
    }
    else {
      asymmetry[bin] = 0.;
      asymmetryError[bin] = 0.;
    }

  } 
 
  writeAsymToFile(anaChoice);
};

void OctetAsymmetry::calcTotalAsymmetry(double enWinLow, double enWinHigh, int anaChoice) {
  if (!isAnaChoiceRateVectors() || getCurrentAnaChoice()!=anaChoice) {
    makeAnalysisChoiceRateVectors(anaChoice);
  }
  unsigned int binLow = (unsigned int)(enWinLow/energyBinWidth);
  unsigned int binHigh = (unsigned int)(enWinHigh/energyBinWidth)-1;

  double sf_ON[2]={0.}, sf_OFF[2]={0.};
  double sf_ON_err[2]={0.}, sf_OFF_err[2]={0.};
  double sumA2[2]={0.}, sumA5[2]={0.}, sumA7[2]={0.}, sumA10[2]={0.}, sumB2[2]={0.}, sumB5[2]={0.}, sumB7[2]={0.}, sumB10[2]={0.};
  double sumA2_err[2]={0.}, sumA5_err[2]={0.}, sumA7_err[2]={0.}, sumA10_err[2]={0.}, sumB2_err[2]={0.}, sumB5_err[2]={0.}, sumB7_err[2]={0.}, sumB10_err[2]={0.};

  double R = 0.;
  double deltaR = 0.;

  for (unsigned int side=0; side<2; side++) {
    for (unsigned int bin=binLow; bin<=binHigh; bin++) {
      sumA2[side]+=anaChoice_A2[side][bin];
      sumA2_err[side]+=power(anaChoice_A2err[side][bin],2);
      sumA5[side]+=anaChoice_A5[side][bin];
      sumA5_err[side]+=power(anaChoice_A5err[side][bin],2);
      sumA7[side]+=anaChoice_A7[side][bin];
      sumA7_err[side]+=power(anaChoice_A7err[side][bin],2);
      sumA10[side]+=anaChoice_A10[side][bin];
      sumA10_err[side]+=power(anaChoice_A10err[side][bin],2);
      sumB2[side]+=anaChoice_B2[side][bin];
      sumB2_err[side]+=power(anaChoice_B2err[side][bin],2);
      sumB5[side]+=anaChoice_B5[side][bin];
      sumB5_err[side]+=power(anaChoice_B5err[side][bin],2);
      sumB7[side]+=anaChoice_B7[side][bin];
      sumB7_err[side]+=power(anaChoice_B7err[side][bin],2);
      sumB10[side]+=anaChoice_B10[side][bin];
      sumB10_err[side]+=power(anaChoice_B10err[side][bin],2);

      /*if (side==1) {
	std::cout << anaChoice_A10[side][bin] << " " << anaChoice_A10err[side][bin] << std::endl;
	}*/ 
    }
    sumA2_err[side] = sqrt(sumA2_err[side]);
    sumA5_err[side] = sqrt(sumA5_err[side]);
    sumA7_err[side] = sqrt(sumA7_err[side]);
    sumA10_err[side] = sqrt(sumA10_err[side]);
    sumB2_err[side] = sqrt(sumB2_err[side]);
    sumB5_err[side] = sqrt(sumB5_err[side]);
    sumB7_err[side] = sqrt(sumB7_err[side]);
    sumB10_err[side] = sqrt(sumB10_err[side]);
    
    /*if (side==1) {
      std::cout << sumA2[side] << " " << sumA2_err[side] << std::endl;
      //std::cout << sumA5[side] << " " << sumA5_err[side] << std::endl;
      //std::cout << sumA7[side] << " " << sumA7_err[side] << std::endl;
      std::cout << sumA10[side] << " " << sumA10_err[side] << std::endl;
      //std::cout << sumB2[side] << " " << sumB2_err[side] << std::endl;
      std::cout << sumB5[side] << " " << sumB5_err[side] << std::endl;
      std::cout << sumB7[side] << " " << sumB7_err[side] << std::endl;
      //   std::cout << sumB10[side] << " " << sumB10_err[side] << std::endl << std::endl;
      }*/
  }
    		   
  for (unsigned int side=0; side<2; side++) {
  
    double weightsum=0.;
    sf_OFF[side] = (power(1./sumA2_err[side],2)*sumA2[side]+power(1./sumA10_err[side],2)*sumA10[side]
		   +power(1./sumB5_err[side],2)*sumB5[side]+power(1./sumB7_err[side],2)*sumB7[side]);
    weightsum = power(1./sumA2_err[side],2)+power(1./sumA10_err[side],2)
      +power(1./sumB5_err[side],2)+power(1./sumB7_err[side],2);
    
    sf_OFF[side] = sf_OFF[side]/weightsum;
    sf_OFF_err[side] = 1./sqrt(weightsum);
    
    weightsum=0.;
    sf_ON[side] = (power(1./sumA5_err[side],2)*sumA5[side]+power(1./sumA7_err[side],2)*sumA7[side]
		  +power(1./sumB2_err[side],2)*sumB2[side]+power(1./sumB10_err[side],2)*sumB10[side]);
    weightsum = power(1./sumA5_err[side],2)+power(1./sumA7_err[side],2)
      +power(1./sumB2_err[side],2)+power(1./sumB10_err[side],2);
    
    sf_ON[side] = sf_ON[side]/weightsum;
    sf_ON_err[side] = 1./sqrt(weightsum);
    std::cout << sf_OFF[side] << " " << sf_ON[side] << std::endl;
  }
  
  R = sf_OFF[0]*sf_ON[1]/(sf_ON[0]*sf_OFF[1]);
  deltaR = sqrt(R*R*(power(sf_OFF_err[0]/sf_OFF[0],2)+power(sf_ON_err[1]/sf_ON[1],2)+power(sf_OFF_err[1]/sf_OFF[1],2)+power(sf_ON_err[0]/sf_ON[0],2)));
  totalAsymmetry = (1.-sqrt(R))/(1+sqrt(R));
  totalAsymmetryError = (deltaR)/(sqrt((R))*power((sqrt((R))+1.),2));

  std::cout << std::endl << R << " " << deltaR << "\n";
  
};


void OctetAsymmetry::calcSuperSum(int anaChoice) {
  if (!isAnaChoiceRateVectors() || getCurrentAnaChoice()!=anaChoice) {
    makeAnalysisChoiceRateVectors(anaChoice);
  }

  double sfON[2]={0.};
  double sfOFF[2]={0.};
  double sfON_err[2]={0.};
  double sfOFF_err[2]={0.};

  for (unsigned int bin=0; bin<superSum.size(); bin++) {
 
    for (unsigned int side=0; side<2; side++) {
      double weightsum=0.;

      sfOFF[side] = (anaChoice_A2[side][bin]>0.?power(1./anaChoice_A2err[side][bin],2)*anaChoice_A2[side][bin]:0.
		     + anaChoice_A10[side][bin]>0.?power(1./anaChoice_A10err[side][bin],2)*anaChoice_A10[side][bin]:0.
		     + anaChoice_B5[side][bin]>0.?power(1./anaChoice_B5err[side][bin],2)*anaChoice_B5[side][bin]:0.
		     + anaChoice_B7[side][bin]>0.?power(1./anaChoice_B7err[side][bin],2)*anaChoice_B7[side][bin]:0.);
      weightsum = (anaChoice_A2[side][bin]>0.?power(1./anaChoice_A2err[side][bin],2):0.
		   + anaChoice_A10[side][bin]>0.?power(1./anaChoice_A10err[side][bin],2):0.
		   + anaChoice_B5[side][bin]>0.?power(1./anaChoice_B5err[side][bin],2):0.
		   + anaChoice_B7[side][bin]>0.?power(1./anaChoice_B7err[side][bin],2):0.);

      sfOFF[side] = weightSum>0. ? sfON[side]/weightsum : 0.;
      sfOFF_err[side] = weightSum>0. ? 1./sqrt(weightsum) : 0.;
      
      weightsum=0.;
      sfON[side] = (anaChoice_A5[side][bin]>0.?power(1./anaChoice_A5err[side][bin],2)*anaChoice_A5[side][bin]:0.
		     + anaChoice_A7[side][bin]>0.?power(1./anaChoice_A7err[side][bin],2)*anaChoice_A7[side][bin]:0.
		     + anaChoice_B2[side][bin]>0.?power(1./anaChoice_B2err[side][bin],2)*anaChoice_B2[side][bin]:0.
		     + anaChoice_B10[side][bin]>0.?power(1./anaChoice_B10err[side][bin],2)*anaChoice_B10[side][bin]:0.);
      weightsum = (anaChoice_A5[side][bin]>0.?power(1./anaChoice_A5err[side][bin],2):0.
		   + anaChoice_A7[side][bin]>0.?power(1./anaChoice_A7err[side][bin],2):0.
		   + anaChoice_B2[side][bin]>0.?power(1./anaChoice_B2err[side][bin],2):0.
		   + anaChoice_B10[side][bin]>0.?power(1./anaChoice_B10err[side][bin],2):0.);


      sfON[side] = weightSum>0. ? sfON[side]/weightsum : 0.;
      sfON_err[side] = weightSum>0. ? 1./sqrt(weightsum) : 0.;

      /*sfOFF[side] = (power(1./anaChoice_A2err[side][bin],2)*anaChoice_A2[side][bin]+power(1./anaChoice_A10err[side][bin],2)*anaChoice_A10[side][bin]
		     +power(1./anaChoice_B5err[side][bin],2)*anaChoice_B5[side][bin]+power(1./anaChoice_B7err[side][bin],2)*anaChoice_B7[side][bin]);
      weightsum = power(1./anaChoice_A2err[side][bin],2)+power(1./anaChoice_A10err[side][bin],2)
      +power(1./anaChoice_B5err[side][bin],2)+power(1./anaChoice_B7err[side][bin],2);

      sfOFF[side] = sfOFF[side]/weightsum;
      sfOFF_err[side] = 1./sqrt(weightsum);

      weightsum=0.;
      sfON[side] = (power(1./anaChoice_A5err[side][bin],2)*anaChoice_A5[side][bin]+power(1./anaChoice_A7err[side][bin],2)*anaChoice_A7[side][bin]
		    +power(1./anaChoice_B2err[side][bin],2)*anaChoice_B2[side][bin]+power(1./anaChoice_B10err[side][bin],2)*anaChoice_B10[side][bin]);
      weightsum = power(1./anaChoice_A5err[side][bin],2)+power(1./anaChoice_A7err[side][bin],2)
	+power(1./anaChoice_B2err[side][bin],2)+power(1./anaChoice_B10err[side][bin],2);
      
      sfON[side] = sfON[side]/weightsum;
      sfON_err[side] = 1./sqrt(weightsum);*/
      
      //if (bin==73 || bin==74) std::cout << sfOFF[side] << " " << sfON[side] << std::endl;
      //if (side==1) {
      //std::cout << anaChoice_A10[side][bin] << " " << anaChoice_A10err[side][bin] << std::endl;
      //}
    }

    superSum[bin] = 0.5*sqrt(sfOFF[0]*sfON[1])+0.5*sqrt(sfON[0]*sfOFF[1]);
    superSumError[bin] = sqrt(0.5*
			      (sfON[0]>0. && sfOFF[1]>0.)?(sfOFF[1]*sfON_err[0]+sfON[0]*sfOFF_err[1])/sqrt(sfON[0]*sfOFF[1]):0.
			      + (sfON[1]>0. && sfOFF[0]>0.)?(sfON[1]*sfOFF_err[0]+sfOFF[0]*sfON_err[1])/sqrt(sfON[1]*sfOFF[0]):0.);
    
  } 
 
  //writeAsymToFile(anaChoice);
};

void OctetAsymmetry::writeAsymToFile(int anaChoice) {
  std::string outpath = std::string(getenv("ANALYSIS_RESULTS")) + "Octet_" + itos(octet) + "/OctetAsymmetry/rawAsymmetry_Octet" + itos(octet)+".dat";
  std::ofstream outfile(outpath.c_str());
  
  for (unsigned int i=0; i<asymmetry.size(); i++) {
    outfile << binLowerEdge[i] << " " << asymmetry[i] << " " << asymmetryError[i] << std::endl;
  }
  outfile.close(); 
};


void makePlots() {

};
    

  
