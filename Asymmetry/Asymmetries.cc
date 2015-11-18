/*
Defines different types of Asymmetries to be calculated within an Octet.
ppair, quartet, octet
*/

#include "Asymmetries.hh"

#include <fstream>

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
///////////////////////////////////////////////////////////////////////////////////////////

OctetAsymmetry::OctetAsymmetry(int oct, double enBinWidth, double fidCut) : AsymmetryBase(oct),
									    energyBinWidth(enBinWidth), fiducialCut(fidCut) {
  unsigned int numBins = (unsigned int)(1200./energyBinWidth);
  numEvtsByTypeByBin.resize(4,std::vector<double>(numBins,0.));
  binLowerEdge.resize(numBins,0.);
  binUpperEdge.resize(numBins,0.);
  
  for (unsigned int i=0; i<numBins; i++) {
    binLowerEdge[i] = (double)(i)*energyBinWidth;
    binUpperEdge[i] = (double)(i+1)*energyBinWidth;
  }  
  std::cout <<"//////////////////////////////////////////////////////////////////\n"
	    <<"Initialized OctetAsymmetry for octet " << octet << std::endl;
};

void OctetAsymmetry::calcBGsubtractedEvts() {
  std::map<std::string,int>::iterator it = runType.begin();
  BGSubtractedRate *bg;
  while (it!=runType.end()) {
    if (checkIfBetaRun(it->first)) {
      bg = new BGSubtractedRate(it->second,energyBinWidth,fiducialCut,true,false);

      std::cout << "initialized BGStubtractedRate for run " << it->second << " (BG run " 
		<< bg->getBackgroundRun(it->second) << ") \n";

      bg->calcBGSubtRates();
      std::vector<double> runLengthBeta = bg->returnRunLengths(true);
      std::vector<double> runLengthBG = bg->returnRunLengths(false);
      std::cout << "RunLength: \tE\tW\n\t\t" 
      << runLengthBeta[0] << "\t" << runLengthBeta[1] << std::endl
      << "\t\t" << runLengthBG[0] << "\t" << runLengthBG[1] << std::endl;

      for (int type=0;type<4;type++) {
	std::vector <double> evecbg = bg->returnBGSubtRate(0,type);
	std::vector <double> wvecbg = bg->returnBGSubtRate(1,type);
	for (unsigned int bin=0; bin<evecbg.size(); bin++) {
	  numEvtsByTypeByBin[type][bin]+=runLengthBeta[0]*evecbg[bin]+runLengthBeta[1]*wvecbg[bin];
	  }
      }
      delete bg;
    }
    it++;
  }
};

double OctetAsymmetry::getNumBGsubtrEvts(double enWinLow, double enWinHigh, int evtType) {
  unsigned int binLow = (unsigned int)(enWinLow/energyBinWidth);
  unsigned int binHigh = (unsigned int)(enWinHigh/energyBinWidth)-1;

  std::cout << "Number of Type " << evtType << " events for energy window " << binLowerEdge[binLow] << " - " << binUpperEdge[binHigh] << ": ";
  double numEvts = 0.;

  for (unsigned int i=binLow; i<=binHigh; i++) {
    numEvts+=numEvtsByTypeByBin[evtType][i];
  }
  std::cout << numEvts << std::endl;
  return numEvts;
};


bool OctetAsymmetry::checkIfBetaRun(std::string type) {
  if (type=="A2" || type=="A5" || type=="A7" || type=="A10" || type=="B2" || type=="B5" || type=="B7" || type=="B10") return true;
  else return false;
};



void makePlots() {

};
    

  
