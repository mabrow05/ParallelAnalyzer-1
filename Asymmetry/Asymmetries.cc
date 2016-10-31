/*
Defines different types of Asymmetries to be calculated within an Octet.
ppair, quartet, octet
*/

#include "Asymmetries.hh"
#include "MBUtils.hh"

#include <fstream>
#include <cstdlib>
#include <cmath>

std::map<int,int> BGrunReplace = {{23017,23006},{21175,21164},{21176,21187}};
// 21175, 21176 beam out in runlog
// 23017 very low statistics

AsymmetryBase::AsymmetryBase(int oct, double enBinWidth, double fidCut, bool ukdata, bool simulation, bool applyAsym, bool unblind) : UKdata(ukdata), Simulation(simulation), applyAsymmetry(applyAsym), UNBLIND(unblind), octet(oct), energyBinWidth(enBinWidth), numEnBins(0), fiducialCut(fidCut), boolAnaChRtVecs(false), runsInOctet(0), analysisChoice(1) {
  numEnBins = (unsigned int)(1200./energyBinWidth);

  if (UNBLIND) {
    std::cout << "***************************************************************\n"
	      << "***************************************************************\n\n"
	      << "                UNBLINDING OCTET " << " " << octet << "\n\n"
	      << "***************************************************************\n"
	      << "***************************************************************\n\n";
  }


  A2.resize(4,std::vector < std::vector <double> > (2,std::vector <double> (numEnBins,0.)));
  B2.resize(4,std::vector < std::vector <double> > (2,std::vector <double> (numEnBins,0.)));
  A5.resize(4,std::vector < std::vector <double> > (2,std::vector <double> (numEnBins,0.)));
  B5.resize(4,std::vector < std::vector <double> > (2,std::vector <double> (numEnBins,0.)));
  A7.resize(4,std::vector < std::vector <double> > (2,std::vector <double> (numEnBins,0.)));
  B7.resize(4,std::vector < std::vector <double> > (2,std::vector <double> (numEnBins,0.)));
  A10.resize(4,std::vector < std::vector <double> > (2,std::vector <double> (numEnBins,0.)));
  B10.resize(4,std::vector < std::vector <double> > (2,std::vector <double> (numEnBins,0.)));
  A2err.resize(4,std::vector < std::vector <double> > (2,std::vector <double> (numEnBins,0.)));
  B2err.resize(4,std::vector < std::vector <double> > (2,std::vector <double> (numEnBins,0.)));
  A5err.resize(4,std::vector < std::vector <double> > (2,std::vector <double> (numEnBins,0.)));
  B5err.resize(4,std::vector < std::vector <double> > (2,std::vector <double> (numEnBins,0.)));
  A7err.resize(4,std::vector < std::vector <double> > (2,std::vector <double> (numEnBins,0.)));
  B7err.resize(4,std::vector < std::vector <double> > (2,std::vector <double> (numEnBins,0.)));
  A10err.resize(4,std::vector < std::vector <double> > (2,std::vector <double> (numEnBins,0.)));
  B10err.resize(4,std::vector < std::vector <double> > (2,std::vector <double> (numEnBins,0.)));
  anaChoice_A2.resize(2,std::vector <double> (numEnBins,0.));
  anaChoice_B2.resize(2,std::vector <double> (numEnBins,0.));
  anaChoice_A5.resize(2,std::vector <double> (numEnBins,0.));
  anaChoice_B5.resize(2,std::vector <double> (numEnBins,0.));
  anaChoice_A7.resize(2,std::vector <double> (numEnBins,0.));
  anaChoice_B7.resize(2,std::vector <double> (numEnBins,0.));
  anaChoice_A10.resize(2,std::vector <double> (numEnBins,0.));
  anaChoice_B10.resize(2,std::vector <double> (numEnBins,0.));
  anaChoice_A2err.resize(2,std::vector <double> (numEnBins,0.));
  anaChoice_B2err.resize(2,std::vector <double> (numEnBins,0.));
  anaChoice_A5err.resize(2,std::vector <double> (numEnBins,0.));
  anaChoice_B5err.resize(2,std::vector <double> (numEnBins,0.));
  anaChoice_A7err.resize(2,std::vector <double> (numEnBins,0.));
  anaChoice_B7err.resize(2,std::vector <double> (numEnBins,0.));
  anaChoice_A10err.resize(2,std::vector <double> (numEnBins,0.));
  anaChoice_B10err.resize(2,std::vector <double> (numEnBins,0.));
  numEvtsEastByTypeByBin.resize(4,std::vector<double>(numEnBins,0.));
  numEvtsWestByTypeByBin.resize(4,std::vector<double>(numEnBins,0.));
  A2len.resize(2,std::vector <double> (2,0.));
  A5len.resize(2,std::vector <double> (2,0.));
  A7len.resize(2,std::vector <double> (2,0.));
  A10len.resize(2,std::vector <double> (2,0.));
  B2len.resize(2,std::vector <double> (2,0.));
  B5len.resize(2,std::vector <double> (2,0.));
  B7len.resize(2,std::vector <double> (2,0.));
  B10len.resize(2,std::vector <double> (2,0.));

  binLowerEdge.resize(numEnBins,0.);
  binUpperEdge.resize(numEnBins,0.);
  
  for (unsigned int i=0; i<numEnBins; i++) {
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
    runType[runNumberHold] = runTypeHold;
    numRuns++;
  }
  infile.close();
  runsInOctet = numRuns;
  std::cout << "Read in octet file for octet " << octet << " with " << runsInOctet << " runs\n";
};

int AsymmetryBase::getBGrun(int run) {
  std::string type = runType.at(run);
  std::string bgType = type=="A2"?"A1" : ( type=="A5"?"A4" : ( type=="A7"?"A9" : ( type=="A10"?"A12" : ( type=="B2"?"B1" : ( type=="B5"?"B4" : ( type=="B7"?"B9" : ( type=="B10"?"B12" : "BAD" )))))));
  
  if (bgType=="BAD") { throw "PASSED RUN WITH BAD TYPE TO getBGrun() IN ASYMMETRY BASE!"; }
  //std::map <int, std::string>::iterator it;
  for ( auto& it: runType ) { //= runType.begin();it!=runType.end(); it++) {

    if ( it.second==bgType ) {
      if ( BGrunReplace.find(it.first)==BGrunReplace.end() ) return it.first;
      else return BGrunReplace[it.first];
    }
  }

  std::cout<< "NO BACKGROUND RUN IN BG SUBTRACTED RATE FOR RUN " << run << "!\n";
  throw "FAILED IN getBGrun(int run)";
};

bool AsymmetryBase::isFullOctet() {
  bool A2=false, B2=false, A5=false, B5=false, A7=false, B7=false, A10=false, B10=false;
  bool A2bg=false, B2bg=false, A5bg=false, B5bg=false, A7bg=false, B7bg=false, A10bg=false, B10bg=false;
  std::map <int, std::string>::iterator it;
  for ( it = runType.begin();it!=runType.end(); it++) {
    if (it->second=="A2") {A2=true; continue;}
    if (it->second=="B2") {B2=true; continue;}
    if (it->second=="A5") {A5=true; continue;}
    if (it->second=="B5") {B5=true; continue;}
    if (it->second=="A7") {A7=true; continue;}
    if (it->second=="B7") {B7=true; continue;}
    if (it->second=="A10") {A10=true; continue;}
    if (it->second=="B10") {B10=true; continue;}
    if (it->second=="A1") {A2bg=true; continue;}
    if (it->second=="B1") {B2bg=true; continue;}
    if (it->second=="A4") {A5bg=true; continue;}
    if (it->second=="B4") {B5bg=true; continue;}
    if (it->second=="A9") {A7bg=true; continue;}
    if (it->second=="B9") {B7bg=true; continue;}
    if (it->second=="A12") {A10bg=true; continue;}
    if (it->second=="B12") {B10bg=true; continue;}
  }
  return (A2 && B2 && A5 && B5 && A7 && B7 && A10 && B10 && A2bg && B2bg && A5bg && B5bg && A7bg && B7bg && A10bg && B10bg);
};

bool AsymmetryBase::isPair(int quartNum, int pairNum) {
  bool A2=false, B2=false, A5=false, B5=false, A7=false, B7=false, A10=false, B10=false;
  bool A2bg=false, B2bg=false, A5bg=false, B5bg=false, A7bg=false, B7bg=false, A10bg=false, B10bg=false;
  std::map<int,std::string>::iterator it;
  for ( it = runType.begin();it!=runType.end(); it++) {
    if (it->second=="A2") {A2=true; continue;}
    if (it->second=="B2") {B2=true; continue;}
    if (it->second=="A5") {A5=true; continue;}
    if (it->second=="B5") {B5=true; continue;}
    if (it->second=="A7") {A7=true; continue;}
    if (it->second=="B7") {B7=true; continue;}
    if (it->second=="A10") {A10=true; continue;}
    if (it->second=="B10") {B10=true; continue;}
    if (it->second=="A1") {A2bg=true; continue;}
    if (it->second=="B1") {B2bg=true; continue;}
    if (it->second=="A4") {A5bg=true; continue;}
    if (it->second=="B4") {B5bg=true; continue;}
    if (it->second=="A9") {A7bg=true; continue;}
    if (it->second=="B9") {B7bg=true; continue;}
    if (it->second=="A12") {A10bg=true; continue;}
    if (it->second=="B12") {B10bg=true; continue;}
  }
  if (quartNum==0) {
    if (pairNum==0) return (A2 && A5 && A2bg && A5bg);
    if (pairNum==1) return (A7 && A10 && A7bg && A10bg);
    else throw "BAD PAIR NUMBER GIVEN in AsymmetryBase::isPair(int QuartNum=0,1, int pairNum=0,1)";
  }
  else if (quartNum==1) {
    if (pairNum==0) return (B2 && B5 && B2bg && B5bg);
    if (pairNum==1) return (B7 && B10 && B7bg && B10bg);
    else throw "BAD PAIR NUMBER GIVEN in AsymmetryBase::isPair(int QuartNum=0,1, int pairNum=0,1)";
  }
  else throw "BAD QUARTET NUMBER GIVEN in AsymmetryBase::isPair(int QuartNum=0,1, int pairNum=0,1)";
  
};

bool AsymmetryBase::isFullQuartet(int quartNum) {
  bool A2=false, B2=false, A5=false, B5=false, A7=false, B7=false, A10=false, B10=false;
  bool A2bg=false, B2bg=false, A5bg=false, B5bg=false, A7bg=false, B7bg=false, A10bg=false, B10bg=false;
  std::map<int,std::string>::iterator it;
  for ( it = runType.begin();it!=runType.end(); it++) {
    if (it->second=="A2") {A2=true; continue;}
    if (it->second=="B2") {B2=true; continue;}
    if (it->second=="A5") {A5=true; continue;}
    if (it->second=="B5") {B5=true; continue;}
    if (it->second=="A7") {A7=true; continue;}
    if (it->second=="B7") {B7=true; continue;}
    if (it->second=="A10") {A10=true; continue;}
    if (it->second=="B10") {B10=true; continue;}
    if (it->second=="A1") {A2bg=true; continue;}
    if (it->second=="B1") {B2bg=true; continue;}
    if (it->second=="A4") {A5bg=true; continue;}
    if (it->second=="B4") {B5bg=true; continue;}
    if (it->second=="A9") {A7bg=true; continue;}
    if (it->second=="B9") {B7bg=true; continue;}
    if (it->second=="A12") {A10bg=true; continue;}
    if (it->second=="B12") {B10bg=true; continue;}
  }
  if (quartNum==0) return (A2 && A5 && A7 && A10 && A2bg && A5bg && A7bg && A10bg);
  else if (quartNum==1) return (B2 && B5 && B7 && B10 && B2bg && B5bg && B7bg && B10bg);
  else throw "BAD QUARTET NUMBER GIVEN in AsymemtryBase::isQuartet(int QuartNum)";
  
};


void AsymmetryBase::loadRates() {
  std::map<int,std::string>::iterator it = runType.begin();
  BGSubtractedRate *bgSubtr;
  while (it!=runType.end()) {
        
    if (checkIfBetaRun(it->second)) {
      bgSubtr = new BGSubtractedRate(it->first,getBGrun(it->first),energyBinWidth,fiducialCut,UKdata,Simulation,applyAsymmetry, UNBLIND); //UNBLIND defaults to false

      std::cout << "initialized BGStubtractedRate for run " << it->first << " (BG run " 
		<< getBGrun(it->first) << ") \n";

      bgSubtr->calcBGSubtRates();

      for (int type=0;type<4;type++) {
	for (int side=0; side<2; side++) {  

	  // vectors to hold the rates bin by bin so that we can do addition of multiple beta runs of the same type.. dumb
	  std::vector < double > rate (numEnBins, 0.);
	  rate = bgSubtr->returnBGSubtRate(side,type);
	  std::vector < double > rateErr (numEnBins, 0.);
	  rateErr = bgSubtr->returnBGSubtRateError(side,type);	  

	  // Note that the len vectors have the issue that they take on the length of the last beta run of that type when 
	  // there are multiple (some octets have split runs, say 2 consecutive A2 runs...)
	  // These lengths aren't currently used for anything though...

	  for (unsigned int bin=0 ; bin<numEnBins ; bin++) {
      
	    if (it->second=="A2") {
	      A2[type][side][bin] += rate[bin];
	      A2err[type][side][bin] += power( rateErr[bin], 2 );
	      A2len[0] = bgSubtr->returnRunLengths(true);
	      A2len[1] = bgSubtr->returnRunLengths(false);
	    }
	    else if (it->second=="B2") {
	      B2[type][side][bin] += rate[bin];
	      B2err[type][side][bin] += power( rateErr[bin], 2 );
	      B2len[0] = bgSubtr->returnRunLengths(true);
	      B2len[1] = bgSubtr->returnRunLengths(false);
	    }
	    else if (it->second=="A5") {
	      A5[type][side][bin] += rate[bin];
	      A5err[type][side][bin] += power( rateErr[bin], 2 );
	      A5len[0] = bgSubtr->returnRunLengths(true);
	      A5len[1] = bgSubtr->returnRunLengths(false);
	    }
	    else if (it->second=="B5") {
	      B5[type][side][bin] += rate[bin];
	      B5err[type][side][bin] += power( rateErr[bin], 2 );
	      B5len[0] = bgSubtr->returnRunLengths(true);
	      B5len[1] = bgSubtr->returnRunLengths(false);
	    }
	    else if (it->second=="A7") {
	      A7[type][side][bin] += rate[bin];
	      A7err[type][side][bin] += power( rateErr[bin], 2 );
	      A7len[0] = bgSubtr->returnRunLengths(true);
	      A7len[1] = bgSubtr->returnRunLengths(false);
	    }
	    else if (it->second=="B7") {
	      B7[type][side][bin] += rate[bin];
	      B7err[type][side][bin] += power( rateErr[bin], 2 );
	      B7len[0] = bgSubtr->returnRunLengths(true);
	      B7len[1] = bgSubtr->returnRunLengths(false);
	    }
	    else if (it->second=="A10") {
	      A10[type][side][bin] += rate[bin];
	      A10err[type][side][bin] += power( rateErr[bin], 2 );
	      A10len[0] = bgSubtr->returnRunLengths(true);
	      A10len[1] = bgSubtr->returnRunLengths(false);
	    }
	    else if (it->second=="B10") {
	      B10[type][side][bin] += rate[bin];
	      B10err[type][side][bin] += power( rateErr[bin], 2 );
	      B10len[0] = bgSubtr->returnRunLengths(true);
	      B10len[1] = bgSubtr->returnRunLengths(false);
	    }
	    else throw "Run misidentified in loadRates";

	  
	  }
	}
      }
	  
      delete bgSubtr;
    }
    it++;
  }
  for (int type=0;type<4;type++) {
    for (int side=0; side<2; side++) {  
      for (unsigned int bin=0 ; bin<numEnBins ; bin++) {
	A2err[type][side][bin] = sqrt(A2err[type][side][bin]);
	B2err[type][side][bin] = sqrt(B2err[type][side][bin]);
	A5err[type][side][bin] = sqrt(A5err[type][side][bin]);
	B5err[type][side][bin] = sqrt(B5err[type][side][bin]);
	A7err[type][side][bin] = sqrt(A7err[type][side][bin]);
	B7err[type][side][bin] = sqrt(B7err[type][side][bin]);
	A10err[type][side][bin] = sqrt(A10err[type][side][bin]);
	B10err[type][side][bin] = sqrt(B10err[type][side][bin]);
      }
    }
  }
};

bool AsymmetryBase::checkIfBetaRun(std::string type) {
  if (type=="A2" || type=="A5" || type=="A7" || type=="A10" || type=="B2" || type=="B5" || type=="B7" || type=="B10") return true;
  else return false;
};

void AsymmetryBase::calcBGsubtractedEvts() {
  std::map<int,std::string>::iterator it = runType.begin();
  BGSubtractedRate *bg;
  while (it!=runType.end()) {
    if (checkIfBetaRun(it->second)) {
      bg = new BGSubtractedRate(it->first,getBGrun(it->first),energyBinWidth,fiducialCut,UKdata,Simulation,applyAsymmetry);

      std::cout << "initialized BGStubtractedRate for run " << it->first << " (BG run " 
		<< getBGrun(it->first) << ") \n";

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
    for (auto &elem : anaChoice_A2) std::fill(elem.begin(), elem.end(), 0.);
    for (auto &elem : anaChoice_B2) std::fill(elem.begin(), elem.end(), 0.);
    for (auto &elem : anaChoice_A5) std::fill(elem.begin(), elem.end(), 0.);
    for (auto &elem : anaChoice_B5) std::fill(elem.begin(), elem.end(), 0.);
    for (auto &elem : anaChoice_A7) std::fill(elem.begin(), elem.end(), 0.);
    for (auto &elem : anaChoice_B7) std::fill(elem.begin(), elem.end(), 0.);
    for (auto &elem : anaChoice_A10) std::fill(elem.begin(), elem.end(), 0.);
    for (auto &elem : anaChoice_B10) std::fill(elem.begin(), elem.end(), 0.);
    for (auto &elem : anaChoice_A2err) std::fill(elem.begin(), elem.end(), 0.);
    for (auto &elem : anaChoice_B2err) std::fill(elem.begin(), elem.end(), 0.);
    for (auto &elem : anaChoice_A5err) std::fill(elem.begin(), elem.end(), 0.);
    for (auto &elem : anaChoice_B5err) std::fill(elem.begin(), elem.end(), 0.);
    for (auto &elem : anaChoice_A7err) std::fill(elem.begin(), elem.end(), 0.);
    for (auto &elem : anaChoice_B7err) std::fill(elem.begin(), elem.end(), 0.);
    for (auto &elem : anaChoice_A10err) std::fill(elem.begin(), elem.end(), 0.);
    for (auto &elem : anaChoice_B10err) std::fill(elem.begin(), elem.end(), 0.);    
  }
  analysisChoice = anaChoice;
  //unsigned int numBins = (unsigned int)(1200./energyBinWidth);
  unsigned int type_low=0, type_high=0;
  if (anaChoice==1) { type_low=0; type_high=2; }
  else if (anaChoice==3 || anaChoice==5) { type_low=0; type_high=3;}
  else if (anaChoice==2) { type_low=0; type_high=1;}
  else if (anaChoice==4) { type_low=0; type_high=0;}
  else if (anaChoice==6) { type_low=1; type_high=1;}
  else if (anaChoice==7 || anaChoice==8 || anaChoice==9) { type_low=2; type_high=3;}
  else throw "Bad Analysis Choice";
  
  for (unsigned int side=0; side<2; side++) {
    for (unsigned int bin=0; bin<numEnBins; bin++) {
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
      }
    
      anaChoice_A2err[side][bin]=sqrt(anaChoice_A2err[side][bin]);
      anaChoice_A5err[side][bin]=sqrt(anaChoice_A5err[side][bin]);
      anaChoice_A7err[side][bin]=sqrt(anaChoice_A7err[side][bin]);
      anaChoice_A10err[side][bin]=sqrt(anaChoice_A10err[side][bin]);
      anaChoice_B2err[side][bin]=sqrt(anaChoice_B2err[side][bin]);
      anaChoice_B5err[side][bin]=sqrt(anaChoice_B5err[side][bin]);
      anaChoice_B7err[side][bin]=sqrt(anaChoice_B7err[side][bin]);
      anaChoice_B10err[side][bin]=sqrt(anaChoice_B10err[side][bin]);
      //std::cout << anaChoice_A2[side][bin] << " " << anaChoice_A2err[side][bin] << std::endl;
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

OctetAsymmetry::OctetAsymmetry(int oct, double enBinWidth, double fidCut, bool ukdata, bool simulation, bool applyAsym, bool unblind) : AsymmetryBase(oct,enBinWidth,fidCut,ukdata,simulation,applyAsym,unblind), totalAsymmetry(0.), totalAsymmetryError(0.), boolSuperSum(false), boolAsymmetry(false) {
  if (isFullOctet()) {
    //unsigned int numBins = (unsigned int)(1200./energyBinWidth);
    asymmetry.resize(numEnBins,0.);
    asymmetryError.resize(numEnBins,0.);
    superSum.resize(numEnBins,0.);
    superSumError.resize(numEnBins,0.);
    loadRates(); // load the rates in the rate vectors for each run
    std::cout <<"//////////////////////////////////////////////////////////////////\n"
	      <<"OctetAsymmetry for octet " << octet << std::endl;
  }
  else { 
    boolAsymmetry=boolSuperSum=true; //Writing "bad octet" to all files so that it isn't used in the future
    for (int anach=1; anach<10; anach++) {
      analysisChoice = anach;
      writeAsymToFile();
      writeSuperSumToFile();
    }
    throw "OCTET NOT A COMPLETE OCTET";
  }
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

      if (anaChoice_A2[side][bin]>0. && anaChoice_A10[side][bin]>0. && anaChoice_B5[side][bin]>0. && anaChoice_B7[side][bin]>0. ) {

	sfOFF[side] = ( anaChoice_A2[side][bin]/power(anaChoice_A2err[side][bin],2) +
			anaChoice_A10[side][bin]/power(anaChoice_A10err[side][bin],2) +
			anaChoice_B5[side][bin]/power(anaChoice_B5err[side][bin],2) +
			anaChoice_B7[side][bin]/power(anaChoice_B7err[side][bin],2) );

	weightsum = ( 1./power(anaChoice_A2err[side][bin],2) + 
		      1./power(anaChoice_A10err[side][bin],2) +
		      1./power(anaChoice_B5err[side][bin],2) +
		      1./power(anaChoice_B7err[side][bin],2) );

      }
      else weightsum = 0.;

      sfOFF[side] = weightsum>0. ? sfOFF[side] / weightsum : 0.;
      sfOFF_err[side] = weightsum>0. ? 1./sqrt(weightsum) : 0.;


      if (anaChoice_A5[side][bin]>0. && anaChoice_A7[side][bin]>0. && anaChoice_B2[side][bin]>0. && anaChoice_B10[side][bin]>0. ) {

	sfON[side] = ( anaChoice_A5[side][bin]/power(anaChoice_A5err[side][bin],2) +
			anaChoice_A7[side][bin]/power(anaChoice_A7err[side][bin],2) +
			anaChoice_B2[side][bin]/power(anaChoice_B2err[side][bin],2) +
			anaChoice_B10[side][bin]/power(anaChoice_B10err[side][bin],2) );

	weightsum = ( 1./power(anaChoice_A5err[side][bin],2) + 
		      1./power(anaChoice_A7err[side][bin],2) +
		      1./power(anaChoice_B2err[side][bin],2) +
		      1./power(anaChoice_B10err[side][bin],2) );


      }
      else weightsum = 0.;
      //if (anaCh==1 && bin<80) std::cout << weightsum << std::endl;                                                                         
      sfON[side] = weightsum>0. ? sfON[side] / weightsum : 0.;
      sfON_err[side] = weightsum>0. ? 1./sqrt(weightsum) : 0.;

      
    }
    if (sfOFF[0]>0. && sfOFF[1]>0. && sfON[0]>0. && sfON[1]>0.) { 
      R = sfOFF[1]*sfON[0]/(sfON[1]*sfOFF[0]);
      deltaR = sqrt(R*R*(power(sfOFF_err[0]/sfOFF[0],2)+power(sfON_err[1]/sfON[1],2)+power(sfOFF_err[1]/sfOFF[1],2)+power(sfON_err[0]/sfON[0],2)));
      asymmetry[bin] = (1.-sqrt(R))/(1+sqrt(R));
      asymmetryError[bin] = (deltaR)/(sqrt(std::abs(R))*power((sqrt(std::abs(R))+1.),2));
    }
    else {
      asymmetry[bin] = 0.;
      asymmetryError[bin] = 0.;
    }

  } 
 
  boolAsymmetry = true;
};

void OctetAsymmetry::calcTotalAsymmetry(double enWinLow, double enWinHigh, int anaChoice) {
  if (!isAnaChoiceRateVectors() || getCurrentAnaChoice()!=anaChoice) {
    makeAnalysisChoiceRateVectors(anaChoice);
  }
  unsigned int binLow = (unsigned int)(enWinLow/energyBinWidth);
  unsigned int binHigh = (unsigned int)(enWinHigh/energyBinWidth)-1;

  double sfON[2]={0.}, sfOFF[2]={0.};
  double sfON_err[2]={0.}, sfOFF_err[2]={0.};
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
    sfOFF[side] = (power(1./sumA2_err[side],2)*sumA2[side]+power(1./sumA10_err[side],2)*sumA10[side] + power(1./sumB5_err[side],2)*sumB5[side]+power(1./sumB7_err[side],2)*sumB7[side]);
    weightsum = power(1./sumA2_err[side],2)+power(1./sumA10_err[side],2) + power(1./sumB5_err[side],2)+power(1./sumB7_err[side],2);
    
    sfOFF[side] = sfOFF[side]/weightsum;
    sfOFF_err[side] = 1./sqrt(weightsum);
    
    weightsum=0.;
    sfON[side] = (power(1./sumA5_err[side],2)*sumA5[side]+power(1./sumA7_err[side],2)*sumA7[side] + power(1./sumB2_err[side],2)*sumB2[side]+power(1./sumB10_err[side],2)*sumB10[side]);
    weightsum = power(1./sumA5_err[side],2)+power(1./sumA7_err[side],2) + power(1./sumB2_err[side],2)+power(1./sumB10_err[side],2);
    
    sfON[side] = sfON[side]/weightsum;
    sfON_err[side] = 1./sqrt(weightsum);
    std::cout << sfOFF[side] << " " << sfON[side] << std::endl;
  }
  
  if (sfOFF[0]>0. && sfOFF[1]>0. && sfON[0]>0. && sfON[1]>0.) { 
    R = sfOFF[1]*sfON[0]/(sfON[1]*sfOFF[0]);
    deltaR = sqrt(R*R*(power(sfOFF_err[0]/sfOFF[0],2)+power(sfON_err[1]/sfON[1],2)+power(sfOFF_err[1]/sfOFF[1],2)+power(sfON_err[0]/sfON[0],2)));
    totalAsymmetry= (1.-sqrt(R))/(1+sqrt(R));
    totalAsymmetryError = (deltaR)/(sqrt(std::abs(R))*power((sqrt(std::abs(R))+1.),2));
  }
  else {
    totalAsymmetry = 0.;
    totalAsymmetryError= 0.;
  }

  std::cout << std::endl << R << " " << deltaR << "\n";
  
};


void OctetAsymmetry::calcSuperSum(int anaChoice) {
  if (!isAnaChoiceRateVectors() || getCurrentAnaChoice()!=anaChoice) {
    makeAnalysisChoiceRateVectors(anaChoice);
  }
  //unsigned int numBins = (unsigned int)(1200./energyBinWidth);
  superSum.resize(numEnBins,0.);
  superSumError.resize(numEnBins,0.);

  double sfON[2]={0.};
  double sfOFF[2]={0.};
  double sfON_err[2]={0.};
  double sfOFF_err[2]={0.};

  for (unsigned int bin=0; bin<superSum.size(); bin++) {
 
    for (unsigned int side=0; side<2; side++) {
      double weightsum=0.;

      sfOFF[side] = (anaChoice_A2err[side][bin]>0.?power(1./anaChoice_A2err[side][bin],2)*anaChoice_A2[side][bin]:0.) + (anaChoice_A10err[side][bin]>0.?power(1./anaChoice_A10err[side][bin],2)*anaChoice_A10[side][bin]:0.) + (anaChoice_B5err[side][bin]>0.?power(1./anaChoice_B5err[side][bin],2)*anaChoice_B5[side][bin]:0.) + (anaChoice_B7err[side][bin]>0.?power(1./anaChoice_B7err[side][bin],2)*anaChoice_B7[side][bin]:0.);
      weightsum = (anaChoice_A2err[side][bin]>0.?power(1./anaChoice_A2err[side][bin],2):0.) + (anaChoice_A10err[side][bin]>0.?power(1./anaChoice_A10err[side][bin],2):0.) + (anaChoice_B5err[side][bin]>0.?power(1./anaChoice_B5err[side][bin],2):0.) + (anaChoice_B7err[side][bin]>0.?power(1./anaChoice_B7err[side][bin],2):0.);

      sfOFF[side] = weightsum>0. ? sfOFF[side]/weightsum : 0.;
      sfOFF_err[side] = weightsum>0. ? 1./sqrt(weightsum) : 0.;

      weightsum=0.;
      sfON[side] = (anaChoice_A5err[side][bin]>0.?power(1./anaChoice_A5err[side][bin],2)*anaChoice_A5[side][bin]:0.) + (anaChoice_A7err[side][bin]>0.?power(1./anaChoice_A7err[side][bin],2)*anaChoice_A7[side][bin]:0.) + (anaChoice_B2err[side][bin]>0.?power(1./anaChoice_B2err[side][bin],2)*anaChoice_B2[side][bin]:0.) + (anaChoice_B10err[side][bin]>0.?power(1./anaChoice_B10err[side][bin],2)*anaChoice_B10[side][bin]:0.);
      weightsum = (anaChoice_A5err[side][bin]>0.?power(1./anaChoice_A5err[side][bin],2):0.) + (anaChoice_A7err[side][bin]>0.?power(1./anaChoice_A7err[side][bin],2):0.) + (anaChoice_B2err[side][bin]>0.?power(1./anaChoice_B2err[side][bin],2):0.) + (anaChoice_B10err[side][bin]>0.?power(1./anaChoice_B10err[side][bin],2):0.);
      
      
      sfON[side] = weightsum>0. ? sfON[side]/weightsum : 0.;
      sfON_err[side] = weightsum>0. ? 1./sqrt(weightsum) : 0.;

    }

    //The geometric mean as defined by http://www.arpapress.com/Volumes/Vol11Issue3/IJRRAS_11_3_08.pdf
    double R1 = (sfON[1]>0. && sfOFF[0]>0.) ? sqrt(sfOFF[0]*sfON[1]) : (sfON[1]<0. && sfOFF[0]<0.) ? -sqrt(sfOFF[0]*sfON[1]) : 0.5*(sfOFF[0]+sfON[1]); 
    double R2 = (sfON[0]>0. && sfOFF[1]>0.) ? sqrt(sfOFF[1]*sfON[0]) : (sfON[0]<0. && sfOFF[1]<0.) ? -sqrt(sfOFF[1]*sfON[0]) : 0.5*(sfOFF[1]+sfON[0]);
    double deltaR1 = ((sfON[1]>0. && sfOFF[0]>0.) || (sfON[1]<0. && sfOFF[0]<0.)) ? 0.5*sqrt((power(sfOFF[0]*sfON_err[1],2)+power(sfON[1]*sfOFF_err[0],2))/(sfON[1]*sfOFF[0])) : 
      0.5*sqrt(power(sfON_err[1],2) + power(sfOFF_err[0],2));
    double deltaR2 = ((sfON[0]>0. && sfOFF[1]>0.) || (sfON[0]<0. && sfOFF[1]<0.)) ? 0.5*sqrt((power(sfOFF[1]*sfON_err[0],2)+power(sfON[0]*sfOFF_err[1],2))/(sfON[0]*sfOFF[1])) : 
      0.5*sqrt(power(sfON_err[0],2) + power(sfOFF_err[1],2));

    superSum[bin] = R1 + R2;
    //std::cout << sfOFF[0] << " " << sfOFF_err[0] << std::endl;
    superSumError[bin] = sqrt(power(deltaR1,2) + power(deltaR2,2));
    
  } 
  boolSuperSum = true;
};

void OctetAsymmetry::writeAsymToFile() {
  if (!isAsymmetry()) {
    std::cout << "No Asymmetry has been calculated. Not writing to file!\n\n";
    return;
  }

  int anaChoice = getCurrentAnaChoice();
  std::string outpath;  
  
  if (Simulation) outpath = std::string(getenv("SIM_ANALYSIS_RESULTS")) + "Octet_" + itos(octet) + "/OctetAsymmetry/"+ (UNBLIND?"UNBLINDED_":"")+"rawAsymmetry_Octet" + itos(octet)+"_AnaCh"+itos(anaChoice)+".dat";
  else if (UKdata) outpath = std::string(getenv("ANALYSIS_RESULTS")) + "Octet_" + itos(octet) + "/OctetAsymmetry/"+ (UNBLIND?"UNBLINDED_":"")+"rawAsymmetry_Octet" + itos(octet)+"_AnaCh"+itos(anaChoice)+".dat";
  else outpath = std::string(getenv("MPM_ANALYSIS_RESULTS")) + "Octet_" + itos(octet) + "/OctetAsymmetry/"+ (UNBLIND?"UNBLINDED_":"")+"rawAsymmetry_Octet" + itos(octet)+"_AnaCh"+itos(anaChoice)+".dat";
  std::ofstream outfile(outpath.c_str());
  
  if (isFullOctet()) {
    for (unsigned int i=0; i<asymmetry.size(); i++) {
      outfile << binLowerEdge[i] << " " << asymmetry[i] << " " << asymmetryError[i] << std::endl;
      //std::cout << binLowerEdge[i] << " " << asymmetry[i] << " " << asymmetryError[i] << std::endl;
    }
  }
  else {
    outfile << "BAD OCTET";
  }
  outfile.close();
  std::cout << "Wrote Asymmetry to file for " << anaChoice << " in " << outpath << "\n";
};

void OctetAsymmetry::writeSuperSumToFile() {
  if (!isSuperSum()) {
    std::cout << "No Super-Sum has been calculated. Not writing to file!\n\n";
    return;
  }

  int anaChoice = getCurrentAnaChoice();
  std::string outpath;

  if (Simulation) outpath = std::string(getenv("SIM_ANALYSIS_RESULTS")) + "Octet_" + itos(octet) + "/OctetAsymmetry/"+ (UNBLIND?"UNBLINDED_":"")+"superSum_Octet" + itos(octet)+"_AnaCh"+itos(anaChoice)+".dat";
  else if (UKdata) outpath = std::string(getenv("ANALYSIS_RESULTS")) + "Octet_" + itos(octet) + "/OctetAsymmetry/"+ (UNBLIND?"UNBLINDED_":"")+"superSum_Octet" + itos(octet)+"_AnaCh"+itos(anaChoice)+".dat";
  else outpath = std::string(getenv("MPM_ANALYSIS_RESULTS")) + "Octet_" + itos(octet) + "/OctetAsymmetry/superSum_Octet"+ (UNBLIND?"UNBLINDED_":"")+"" + itos(octet)+"_AnaCh"+itos(anaChoice)+".dat";

  
  std::ofstream outfile(outpath.c_str());
  
  if (isFullOctet()) {
    for (unsigned int i=0; i<superSum.size(); i++) {
      outfile << binLowerEdge[i] << " " << superSum[i] << " " << superSumError[i] << std::endl;
      //std::cout << binLowerEdge[i] << " " << superSum[i] << " " << superSumError[i] << std::endl;
    }
  }
  else {
    outfile << "BAD OCTET";
  }

  outfile.close(); 
};

    

///////////////////////////////////////////////////////////////////////////////////////////

QuartetAsymmetry::QuartetAsymmetry(int oct, double enBinWidth, double fidCut, bool ukdata, bool simulation, bool applyAsym, bool unblind) : AsymmetryBase(oct,enBinWidth,fidCut,ukdata,simulation,applyAsym,unblind), totalAsymmetryA(0.), totalAsymmetryErrorA(0.), totalAsymmetryB(0.), totalAsymmetryErrorB(0.), boolSuperSum(false), boolAsymmetry(false) {

  isGoodQuartet.push_back(isFullQuartet(0));
  isGoodQuartet.push_back(isFullQuartet(1));

  if (isGoodQuartet[0] || isGoodQuartet[1]) {
    //unsigned int numBins = (unsigned int)(1200./energyBinWidth);
    asymmetry.resize(2, std::vector<double> (numEnBins,0.));
    asymmetryError.resize(2, std::vector<double> (numEnBins,0.));
    superSum.resize(2, std::vector<double> (numEnBins,0.));
    superSumError.resize(2, std::vector<double> (numEnBins,0.));
    loadRates(); // load the rates in the rate vectors for each run
    std::cout <<"//////////////////////////////////////////////////////////////////\n"
	      <<"QuartetAsymmetry for octet " << octet << std::endl
	      <<"\nQuartet Status (0=bad, 1=good): First: " << isGoodQuartet[0] << " " 
	      <<"Second: " << isGoodQuartet[1] << std::endl;
  }
  else {
    boolAsymmetry=boolSuperSum=true; //Writing "bad quartet" to all files so that it isn't used in the future
    for (int anach=1; anach<10; anach++) {
      analysisChoice = anach;
      writeAsymToFile();
      writeSuperSumToFile();
    }
    throw "NO QUARTETS IN THIS OCTET";
  }
};

void QuartetAsymmetry::calcAsymmetryBinByBin(int anaChoice) {
  if (!isAnaChoiceRateVectors() || getCurrentAnaChoice()!=anaChoice) {
    makeAnalysisChoiceRateVectors(anaChoice);
  }

  double sfON[2]={0.};
  double sfOFF[2]={0.};
  double sfON_err[2]={0.};
  double sfOFF_err[2]={0.};

  // A type runs
  if (isGoodQuartet[0]) {
    for (unsigned int bin=0; bin<asymmetry[0].size(); bin++) {
      double R = 0.;
      double deltaR = 0.;
      
      for (unsigned int side=0; side<2; side++) {
	double weightsum=0.;
	
	// AFP Off
	sfOFF[side] = (anaChoice_A2err[side][bin]>0.?power(1./anaChoice_A2err[side][bin],2)*anaChoice_A2[side][bin]:0.) + (anaChoice_A10err[side][bin]>0.?power(1./anaChoice_A10err[side][bin],2)*anaChoice_A10[side][bin]:0.);
	weightsum = (anaChoice_A2err[side][bin]>0.?power(1./anaChoice_A2err[side][bin],2):0.) + (anaChoice_A10err[side][bin]>0.?power(1./anaChoice_A10err[side][bin],2):0.);
	
	sfOFF[side] = weightsum>0. ? sfOFF[side]/weightsum : 0.;
	sfOFF_err[side] = weightsum>0. ? 1./sqrt(weightsum) : 0.;
	
	weightsum=0.;
	
	// AFP ON
	sfON[side] = (anaChoice_A5err[side][bin]>0.?power(1./anaChoice_A5err[side][bin],2)*anaChoice_A5[side][bin]:0.) + (anaChoice_A7err[side][bin]>0.?power(1./anaChoice_A7err[side][bin],2)*anaChoice_A7[side][bin]:0.);
	weightsum = (anaChoice_A5err[side][bin]>0.?power(1./anaChoice_A5err[side][bin],2):0.) + (anaChoice_A7err[side][bin]>0.?power(1./anaChoice_A7err[side][bin],2):0.);
	
	sfON[side] = weightsum>0. ? sfON[side]/weightsum : 0.;
	sfON_err[side] = weightsum>0. ? 1./sqrt(weightsum) : 0.;
	
      }
      if (sfOFF[0]>0. && sfOFF[1]>0. && sfON[0]>0. && sfON[1]>0.) { 
	R = sfOFF[1]*sfON[0]/(sfON[1]*sfOFF[0]);
	deltaR = sqrt(R*R*(power(sfOFF_err[0]/sfOFF[0],2)+power(sfON_err[1]/sfON[1],2)+power(sfOFF_err[1]/sfOFF[1],2)+power(sfON_err[0]/sfON[0],2)));
	asymmetry[0][bin] = (1.-sqrt(R))/(1+sqrt(R));
	asymmetryError[0][bin] = (deltaR)/(sqrt(std::abs(R))*power((sqrt(std::abs(R))+1.),2));
      }
      else {
	asymmetry[0][bin] = 0.;
	asymmetryError[0][bin] = 0.;
      }
    }
  }

  // B type runs
  if (isGoodQuartet[1]) {
    for (unsigned int bin=0; bin<asymmetry[1].size(); bin++) {
      double R = 0.;
      double deltaR = 0.;
      
      for (unsigned int side=0; side<2; side++) {
	double weightsum=0.;
	
	sfOFF[side]=0.;
	sfON[side]=0.;
	sfOFF_err[side]=0.;
	sfON_err[side]=0.;

	// AFP Off
	sfOFF[side] = (anaChoice_B5err[side][bin]>0.?power(1./anaChoice_B5err[side][bin],2)*anaChoice_B5[side][bin]:0.) + (anaChoice_B7err[side][bin]>0.?power(1./anaChoice_B7err[side][bin],2)*anaChoice_B7[side][bin]:0.);
	weightsum = (anaChoice_B5err[side][bin]>0.?power(1./anaChoice_B5err[side][bin],2):0.) + (anaChoice_B7err[side][bin]>0.?power(1./anaChoice_B7err[side][bin],2):0.);
	
	sfOFF[side] = weightsum>0. ? sfOFF[side]/weightsum : 0.;
	sfOFF_err[side] = weightsum>0. ? 1./sqrt(weightsum) : 0.;
	
	weightsum=0.;
	
	// AFP ON
	sfON[side] = (anaChoice_B2err[side][bin]>0.?power(1./anaChoice_B2err[side][bin],2)*anaChoice_B2[side][bin]:0.) + (anaChoice_B10err[side][bin]>0.?power(1./anaChoice_B10err[side][bin],2)*anaChoice_B10[side][bin]:0.);

	weightsum = (anaChoice_B2err[side][bin]>0.?power(1./anaChoice_B2err[side][bin],2):0.) + (anaChoice_B10err[side][bin]>0.?power(1./anaChoice_B10err[side][bin],2):0.);
	
	sfON[side] = weightsum>0. ? sfON[side]/weightsum : 0.;
	sfON_err[side] = weightsum>0. ? 1./sqrt(weightsum) : 0.;
	
      }
      if (sfOFF[0]>0. && sfOFF[1]>0. && sfON[0]>0. && sfON[1]>0.) { 
	R = sfOFF[1]*sfON[0]/(sfON[1]*sfOFF[0]);
	deltaR = sqrt(R*R*(power(sfOFF_err[0]/sfOFF[0],2)+power(sfON_err[1]/sfON[1],2)+power(sfOFF_err[1]/sfOFF[1],2)+power(sfON_err[0]/sfON[0],2)));
	asymmetry[1][bin] = (1.-sqrt(R))/(1+sqrt(R));
	asymmetryError[1][bin] = (deltaR)/(sqrt(std::abs(R))*power((sqrt(std::abs(R))+1.),2));
      }
      else {
	asymmetry[1][bin] = 0.;
	asymmetryError[1][bin] = 0.;
      }
    }
  }
  boolAsymmetry = true;
};

void QuartetAsymmetry::calcTotalAsymmetry(double enWinLow, double enWinHigh, int anaChoice) {
  if (!isAnaChoiceRateVectors() || getCurrentAnaChoice()!=anaChoice) {
    makeAnalysisChoiceRateVectors(anaChoice);
  }
  unsigned int binLow = (unsigned int)(enWinLow/energyBinWidth);
  unsigned int binHigh = (unsigned int)(enWinHigh/energyBinWidth)-1;

  double sfON[2]={0.}, sfOFF[2]={0.};
  double sfON_err[2]={0.}, sfOFF_err[2]={0.};
  double sumA2[2]={0.}, sumA5[2]={0.}, sumA7[2]={0.}, sumA10[2]={0.}, sumB2[2]={0.}, sumB5[2]={0.}, sumB7[2]={0.}, sumB10[2]={0.};
  double sumA2_err[2]={0.}, sumA5_err[2]={0.}, sumA7_err[2]={0.}, sumA10_err[2]={0.}, sumB2_err[2]={0.}, sumB5_err[2]={0.}, sumB7_err[2]={0.}, sumB10_err[2]={0.};

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

    }
    sumA2_err[side] = sqrt(sumA2_err[side]);
    sumA5_err[side] = sqrt(sumA5_err[side]);
    sumA7_err[side] = sqrt(sumA7_err[side]);
    sumA10_err[side] = sqrt(sumA10_err[side]);
    sumB2_err[side] = sqrt(sumB2_err[side]);
    sumB5_err[side] = sqrt(sumB5_err[side]);
    sumB7_err[side] = sqrt(sumB7_err[side]);
    sumB10_err[side] = sqrt(sumB10_err[side]);
    
  }
    
  double R = 0.;
  double deltaR = 0.; 

  if (isGoodQuartet[0]) {
    for (unsigned int side=0; side<2; side++) {
      
      double weightsum=0.;
      sfOFF[side] = (power(1./sumA2_err[side],2)*sumA2[side]+power(1./sumA10_err[side],2)*sumA10[side]);
      weightsum = power(1./sumA2_err[side],2)+power(1./sumA10_err[side],2);
      
      sfOFF[side] = sfOFF[side]/weightsum;
      sfOFF_err[side] = 1./sqrt(weightsum);
      
      weightsum=0.;
      sfON[side] = (power(1./sumA5_err[side],2)*sumA5[side]+power(1./sumA7_err[side],2)*sumA7[side]);
      weightsum = power(1./sumA5_err[side],2)+power(1./sumA7_err[side],2);
      
      sfON[side] = sfON[side]/weightsum;
      sfON_err[side] = 1./sqrt(weightsum);
      //std::cout << sfOFF[side] << " " << sfON[side] << std::endl;
    }
    
    if (sfOFF[0]>0. && sfOFF[1]>0. && sfON[0]>0. && sfON[1]>0.) { 
      R = sfOFF[1]*sfON[0]/(sfON[1]*sfOFF[0]);
      deltaR = sqrt(R*R*(power(sfOFF_err[0]/sfOFF[0],2)+power(sfON_err[1]/sfON[1],2)+power(sfOFF_err[1]/sfOFF[1],2)+power(sfON_err[0]/sfON[0],2)));
      totalAsymmetryA = (1.-sqrt(R))/(1+sqrt(R));
      totalAsymmetryErrorA = (deltaR)/(sqrt(std::abs(R))*power((sqrt(std::abs(R))+1.),2));
    }
    else {
      totalAsymmetryA = 0.;
      totalAsymmetryErrorA= 0.;
    }
    
    //std::cout << std::endl << R << " " << deltaR << "\n";
  }

  R=0.;
  deltaR=0.;

  if (isGoodQuartet[1]) {
    for (unsigned int side=0; side<2; side++) {

      sfOFF[side]=0.;
      sfON[side]=0.;
      sfOFF_err[side]=0.;
      sfON_err[side]=0.;
      
      double weightsum=0.;
      sfOFF[side] = (power(1./sumB5_err[side],2)*sumB5[side]+power(1./sumB7_err[side],2)*sumB7[side]);
      weightsum = power(1./sumB5_err[side],2)+power(1./sumB7_err[side],2);
      
      sfOFF[side] = sfOFF[side]/weightsum;
      sfOFF_err[side] = 1./sqrt(weightsum);
      
      weightsum=0.;
      sfON[side] = (power(1./sumB2_err[side],2)*sumB2[side]+power(1./sumB10_err[side],2)*sumB10[side]);
      weightsum = power(1./sumB2_err[side],2)+power(1./sumB10_err[side],2);
      
      sfON[side] = sfON[side]/weightsum;
      sfON_err[side] = 1./sqrt(weightsum);
      //std::cout << sfOFF[side] << " " << sfON[side] << std::endl;
    }
    
    if (sfOFF[0]>0. && sfOFF[1]>0. && sfON[0]>0. && sfON[1]>0.) { 
      R = sfOFF[1]*sfON[0]/(sfON[1]*sfOFF[0]);
      deltaR = sqrt(R*R*(power(sfOFF_err[0]/sfOFF[0],2)+power(sfON_err[1]/sfON[1],2)+power(sfOFF_err[1]/sfOFF[1],2)+power(sfON_err[0]/sfON[0],2)));
      totalAsymmetryB= (1.-sqrt(R))/(1+sqrt(R));
      totalAsymmetryErrorB = (deltaR)/(sqrt(std::abs(R))*power((sqrt(std::abs(R))+1.),2));
    }
    else {
      totalAsymmetryB = 0.;
      totalAsymmetryErrorB= 0.;
    }
    //std::cout << std::endl << R << " " << deltaR << "\n";
  }

  
};


void QuartetAsymmetry::calcSuperSum(int anaChoice) {
  if (!isAnaChoiceRateVectors() || getCurrentAnaChoice()!=anaChoice) {
    makeAnalysisChoiceRateVectors(anaChoice);
  }
  //unsigned int numBins = (unsigned int)(1200./energyBinWidth);
  superSum.resize(2, std::vector<double> (numEnBins,0.));
  superSumError.resize(2, std::vector<double> (numEnBins,0.));

  double sfON[2]={0.};
  double sfOFF[2]={0.};
  double sfON_err[2]={0.};
  double sfOFF_err[2]={0.};

  if (isGoodQuartet[0]) {
    for (unsigned int bin=0; bin<superSum[0].size(); bin++) {
      
      for (unsigned int side=0; side<2; side++) {
	double weightsum=0.;

	// AFP Off
	sfOFF[side] = (anaChoice_A2err[side][bin]>0.?power(1./anaChoice_A2err[side][bin],2)*anaChoice_A2[side][bin]:0.) + (anaChoice_A10err[side][bin]>0.?power(1./anaChoice_A10err[side][bin],2)*anaChoice_A10[side][bin]:0.);
	weightsum = (anaChoice_A2err[side][bin]>0.?power(1./anaChoice_A2err[side][bin],2):0.) + (anaChoice_A10err[side][bin]>0.?power(1./anaChoice_A10err[side][bin],2):0.);
	
	sfOFF[side] = weightsum>0. ? sfOFF[side]/weightsum : 0.;
	sfOFF_err[side] = weightsum>0. ? 1./sqrt(weightsum) : 0.;
	
	weightsum=0.;
	
	// AFP ON
	sfON[side] = (anaChoice_A5err[side][bin]>0.?power(1./anaChoice_A5err[side][bin],2)*anaChoice_A5[side][bin]:0.) + (anaChoice_A7err[side][bin]>0.?power(1./anaChoice_A7err[side][bin],2)*anaChoice_A7[side][bin]:0.);
	weightsum = (anaChoice_A5err[side][bin]>0.?power(1./anaChoice_A5err[side][bin],2):0.) + (anaChoice_A7err[side][bin]>0.?power(1./anaChoice_A7err[side][bin],2):0.);
	
	sfON[side] = weightsum>0. ? sfON[side]/weightsum : 0.;
	sfON_err[side] = weightsum>0. ? 1./sqrt(weightsum) : 0.;
	
      }
      
      //The geometric mean as defined by http://www.arpapress.com/Volumes/Vol11Issue3/IJRRAS_11_3_08.pdf
      double R1 = (sfON[1]>0. && sfOFF[0]>0.) ? sqrt(sfOFF[0]*sfON[1]) : (sfON[1]<0. && sfOFF[0]<0.) ? -sqrt(sfOFF[0]*sfON[1]) : 0.5*(sfOFF[0]+sfON[1]); 
      double R2 = (sfON[0]>0. && sfOFF[1]>0.) ? sqrt(sfOFF[1]*sfON[0]) : (sfON[0]<0. && sfOFF[1]<0.) ? -sqrt(sfOFF[1]*sfON[0]) : 0.5*(sfOFF[1]+sfON[0]);
      double deltaR1 = ((sfON[1]>0. && sfOFF[0]>0.) || (sfON[1]<0. && sfOFF[0]<0.)) ? 0.5*sqrt((power(sfOFF[0]*sfON_err[1],2)+power(sfON[1]*sfOFF_err[0],2))/(sfON[1]*sfOFF[0])) : 
	0.5*sqrt(power(sfON_err[1],2) + power(sfOFF_err[0],2));
      double deltaR2 = ((sfON[0]>0. && sfOFF[1]>0.) || (sfON[0]<0. && sfOFF[1]<0.)) ? 0.5*sqrt((power(sfOFF[1]*sfON_err[0],2)+power(sfON[0]*sfOFF_err[1],2))/(sfON[0]*sfOFF[1])) : 
	0.5*sqrt(power(sfON_err[0],2) + power(sfOFF_err[1],2));
      
      superSum[0][bin] = R1 + R2;
      //std::cout << sfOFF[0] << " " << sfOFF_err[0] << std::endl;
      superSumError[0][bin] = sqrt(power(deltaR1,2) + power(deltaR2,2));
    } 
  }

  if (isGoodQuartet[1]) {
    for (unsigned int bin=0; bin<superSum[1].size(); bin++) {
      
      for (unsigned int side=0; side<2; side++) {

	sfOFF[side]=0.;
	sfON[side]=0.;
	sfOFF_err[side]=0.;
	sfON_err[side]=0.;
	double weightsum=0.;

	// AFP Off
	sfOFF[side] = (anaChoice_B5err[side][bin]>0.?power(1./anaChoice_B5err[side][bin],2)*anaChoice_B5[side][bin]:0.) + (anaChoice_B7err[side][bin]>0.?power(1./anaChoice_B7err[side][bin],2)*anaChoice_B7[side][bin]:0.);
	weightsum = (anaChoice_B5err[side][bin]>0.?power(1./anaChoice_B5err[side][bin],2):0.) + (anaChoice_B7err[side][bin]>0.?power(1./anaChoice_B7err[side][bin],2):0.);
	
	sfOFF[side] = weightsum>0. ? sfOFF[side]/weightsum : 0.;
	sfOFF_err[side] = weightsum>0. ? 1./sqrt(weightsum) : 0.;
	
	weightsum=0.;
	
	// AFP ON
	sfON[side] = (anaChoice_B2err[side][bin]>0.?power(1./anaChoice_B2err[side][bin],2)*anaChoice_B2[side][bin]:0.) + (anaChoice_B10err[side][bin]>0.?power(1./anaChoice_B10err[side][bin],2)*anaChoice_B10[side][bin]:0.);

	weightsum = (anaChoice_B2err[side][bin]>0.?power(1./anaChoice_B2err[side][bin],2):0.) + (anaChoice_B10err[side][bin]>0.?power(1./anaChoice_B10err[side][bin],2):0.);
	
	sfON[side] = weightsum>0. ? sfON[side]/weightsum : 0.;
	sfON_err[side] = weightsum>0. ? 1./sqrt(weightsum) : 0.;
	
      }

      //The geometric mean as defined by http://www.arpapress.com/Volumes/Vol11Issue3/IJRRAS_11_3_08.pdf
      double R1 = (sfON[1]>0. && sfOFF[0]>0.) ? sqrt(sfOFF[0]*sfON[1]) : (sfON[1]<0. && sfOFF[0]<0.) ? -sqrt(sfOFF[0]*sfON[1]) : 0.5*(sfOFF[0]+sfON[1]); 
      double R2 = (sfON[0]>0. && sfOFF[1]>0.) ? sqrt(sfOFF[1]*sfON[0]) : (sfON[0]<0. && sfOFF[1]<0.) ? -sqrt(sfOFF[1]*sfON[0]) : 0.5*(sfOFF[1]+sfON[0]);
      double deltaR1 = ((sfON[1]>0. && sfOFF[0]>0.) || (sfON[1]<0. && sfOFF[0]<0.)) ? 0.5*sqrt((power(sfOFF[0]*sfON_err[1],2)+power(sfON[1]*sfOFF_err[0],2))/(sfON[1]*sfOFF[0])) : 
	0.5*sqrt(power(sfON_err[1],2) + power(sfOFF_err[0],2));
      double deltaR2 = ((sfON[0]>0. && sfOFF[1]>0.) || (sfON[0]<0. && sfOFF[1]<0.)) ? 0.5*sqrt((power(sfOFF[1]*sfON_err[0],2)+power(sfON[0]*sfOFF_err[1],2))/(sfON[0]*sfOFF[1])) : 
	0.5*sqrt(power(sfON_err[0],2) + power(sfOFF_err[1],2));
      
      superSum[1][bin] = R1 + R2;
      //std::cout << sfOFF[0] << " " << sfOFF_err[0] << std::endl;
      superSumError[1][bin] = sqrt(power(deltaR1,2) + power(deltaR2,2));
      
    } 
  }
 
  boolSuperSum = true;
};

void QuartetAsymmetry::writeAsymToFile() {
  if (!isAsymmetry()) {
    std::cout << "No Asymmetry has been calculated. Not writing to file!\n\n";
    return;
  }

  int anaChoice = getCurrentAnaChoice();
  //Setting paths to output files
  std::string outpathA = (Simulation ? std::string(getenv("SIM_ANALYSIS_RESULTS")) : UKdata ? std::string(getenv("ANALYSIS_RESULTS")) : std::string(getenv("MPM_ANALYSIS_RESULTS"))) + 
    "Octet_" + itos(octet) + "/QuartetAsymmetry/"+ (UNBLIND?"UNBLINDED_":"")+"rawAsymmetry_Octet" + itos(octet)+"_AnaCh"+itos(anaChoice)+"_Quartet_A.dat";
  std::string outpathB = (Simulation ? std::string(getenv("SIM_ANALYSIS_RESULTS")) : UKdata ? std::string(getenv("ANALYSIS_RESULTS")) : std::string(getenv("MPM_ANALYSIS_RESULTS"))) 
    + "Octet_" + itos(octet) + "/QuartetAsymmetry/"+ (UNBLIND?"UNBLINDED_":"")+"rawAsymmetry_Octet" + itos(octet)+"_AnaCh"+itos(anaChoice)+"_Quartet_B.dat";
  

  //Open and fill file for quartet A
  std::ofstream outfileA(outpathA.c_str()); 
  
  if (isGoodQuartet[0]) {
    for (unsigned int i=0; i<asymmetry[0].size(); i++) {
      outfileA << binLowerEdge[i] << " " << asymmetry[0][i] << " " << asymmetryError[0][i] << std::endl;
      //std::cout << binLowerEdge[i] << " " << asymmetry[i] << " " << asymmetryError[i] << std::endl;
    }
  }
  else outfileA << "BAD QUARTET";
  outfileA.close();

  //Open and fill file for quartet B
  std::ofstream outfileB(outpathB.c_str());

  if (isGoodQuartet[1]) {
    for (unsigned int i=0; i<asymmetry[1].size(); i++) {
      outfileB << binLowerEdge[i] << " " << asymmetry[1][i] << " " << asymmetryError[1][i] << std::endl;
      //std::cout << binLowerEdge[i] << " " << asymmetry[i] << " " << asymmetryError[i] << std::endl;
    }
  }
  else outfileB << "BAD QUARTET";
  outfileB.close();

  std::cout << "Wrote Asymmetry to file for anaChoice " << anaChoice << " in:\n" << outpathA << std::endl << outpathB << "\n";
};

void QuartetAsymmetry::writeSuperSumToFile() {
  if (!isSuperSum()) {
    std::cout << "No Super-sum has been calculated. Not writing to file!\n\n";
    return;
  }

  int anaChoice = getCurrentAnaChoice();
  //Setting paths to output files
  std::string outpathA = Simulation ? std::string(getenv("SIM_ANALYSIS_RESULTS")) : UKdata ? std::string(getenv("ANALYSIS_RESULTS")) : std::string(getenv("MPM_ANALYSIS_RESULTS")) + 
    "Octet_" + itos(octet) + "/QuartetAsymmetry/"+ (UNBLIND?"UNBLINDED_":"")+"superSum_Octet" + itos(octet)+"_AnaCh"+itos(anaChoice)+"_Quartet_A.dat";
  std::string outpathB = Simulation ? std::string(getenv("SIM_ANALYSIS_RESULTS")) : UKdata ? std::string(getenv("ANALYSIS_RESULTS")) : std::string(getenv("MPM_ANALYSIS_RESULTS")) 
    + "Octet_" + itos(octet) + "/QuartetAsymmetry/"+ (UNBLIND?"UNBLINDED_":"")+"superSum_Octet" + itos(octet)+"_AnaCh"+itos(anaChoice)+"_Quartet_B.dat";

  //Open and fill file for quartet A
  std::ofstream outfileA(outpathA.c_str()); 
  
  if (isGoodQuartet[0]) {
    for (unsigned int i=0; i<superSum[0].size(); i++) {
      outfileA << binLowerEdge[i] << " " << superSum[0][i] << " " << superSumError[0][i] << std::endl;
      //std::cout << binLowerEdge[i] << " " << superSum[i] << " " << superSumError[i] << std::endl;
    }
  }
  else outfileA << "BAD QUARTET";
  outfileA.close();

  //Open and fill file for quartet B
  std::ofstream outfileB(outpathB.c_str());

  if (isGoodQuartet[1]) {
    for (unsigned int i=0; i<superSum[1].size(); i++) {
      outfileB << binLowerEdge[i] << " " << superSum[1][i] << " " << superSumError[1][i] << std::endl;
      //std::cout << binLowerEdge[i] << " " << superSum[i] << " " << superSumError[i] << std::endl;
    }
  }
  else outfileB << "BAD QUARTET";
  outfileB.close();

};


///////////////////////////////////////////////////////////////////////////////////////////

PairAsymmetry::PairAsymmetry(int oct, double enBinWidth, double fidCut, bool ukdata, bool simulation, bool applyAsym, bool unblind) : AsymmetryBase(oct,enBinWidth,fidCut,ukdata,simulation,applyAsym,unblind), totalAsymmetryA0(0.), totalAsymmetryErrorA0(0.), totalAsymmetryB0(0.), totalAsymmetryErrorB0(0.), totalAsymmetryA1(0.), totalAsymmetryErrorA1(0.), totalAsymmetryB1(0.), totalAsymmetryErrorB1(0.), boolSuperSum(false), boolAsymmetry(false) {

  isGoodPair.resize(2,std::vector <bool> (2));

  isGoodPair[0][0] = isPair(0,0);
  isGoodPair[0][1] = isPair(0,1);
  isGoodPair[1][0] = isPair(1,0);
  isGoodPair[1][1] = isPair(1,1);


  if (isGoodPair[0][0] || isGoodPair[1][0] || isGoodPair[0][1] || isGoodPair[1][1]) {
    //unsigned int numBins = (unsigned int)(1200./energyBinWidth);
    asymmetry.resize(2, std::vector < std::vector <double> > (2, std::vector <double>(numEnBins,0.)));
    asymmetryError.resize(2, std::vector < std::vector <double> > (2, std::vector <double>(numEnBins,0.)));
    superSum.resize(2, std::vector < std::vector <double> > (2, std::vector <double>(numEnBins,0.)));
    superSumError.resize(2, std::vector < std::vector <double> > (2, std::vector <double>(numEnBins,0.)));
    loadRates(); // load the rates in the rate vectors for each run
    std::cout <<"//////////////////////////////////////////////////////////////////\n"
	      <<"PairAsymmetry for octet " << octet << std::endl
	      <<"\nPair Status (0=bad, 1=good): First: " << isGoodPair[0][0] << " " 
	      <<"Second: " << isGoodPair[0][1]  << " " <<"Third: " << isGoodPair[1][0] << " "
	      <<"Fourth: " << isGoodPair[1][1] << std::endl;
  }
  else {
    boolAsymmetry=boolSuperSum=true; //Writing "bad pair" to all files so that it isn't used in the future
    for (int anach=1; anach<10; anach++) {
      analysisChoice = anach;
      writeAsymToFile();
      writeSuperSumToFile();
    }
    throw "NO PAIRS IN THIS OCTET";
  }
};

void PairAsymmetry::calcAsymmetryBinByBin(int anaChoice) {
  if (!isAnaChoiceRateVectors() || getCurrentAnaChoice()!=anaChoice) {
    makeAnalysisChoiceRateVectors(anaChoice);
  }

  double sfON[2]={0.};
  double sfOFF[2]={0.};
  double sfON_err[2]={0.};
  double sfOFF_err[2]={0.};

  // A type runs, pair 0
  if (isGoodPair[0][0]) {
    for (unsigned int bin=0; bin<asymmetry[0][0].size(); bin++) {
      double R = 0.;
      double deltaR = 0.;
      
      for (unsigned int side=0; side<2; side++) {
	double weightsum=0.;
	
	// AFP Off
	sfOFF[side] = anaChoice_A2[side][bin];
	weightsum = (anaChoice_A2err[side][bin]>0.?power(1./anaChoice_A2err[side][bin],2):0.);
	
	sfOFF_err[side] = weightsum>0. ? 1./sqrt(weightsum) : 0.;
	
	weightsum=0.;
	
	// AFP ON
	sfON[side] = anaChoice_A5[side][bin];
	weightsum = (anaChoice_A5err[side][bin]>0.?power(1./anaChoice_A5err[side][bin],2):0.);
	
	sfON_err[side] = weightsum>0. ? 1./sqrt(weightsum) : 0.;
	
      }
      if (sfOFF[0]>0. && sfOFF[1]>0. && sfON[0]>0. && sfON[1]>0.) { 
	R = sfOFF[1]*sfON[0]/(sfON[1]*sfOFF[0]);
	deltaR = sqrt(R*R*(power(sfOFF_err[0]/sfOFF[0],2)+power(sfON_err[1]/sfON[1],2)+power(sfOFF_err[1]/sfOFF[1],2)+power(sfON_err[0]/sfON[0],2)));
	asymmetry[0][0][bin] = (1.-sqrt(R))/(1+sqrt(R));
	asymmetryError[0][0][bin] = (deltaR)/(sqrt(std::abs(R))*power((sqrt(std::abs(R))+1.),2));
      }
      else {
	asymmetry[0][0][bin] = 0.;
	asymmetryError[0][0][bin] = 0.;
      }
    }
  }

  if (isGoodPair[0][1]) {
    for (unsigned int bin=0; bin<asymmetry[0][1].size(); bin++) {
      double R = 0.;
      double deltaR = 0.;
      
      for (unsigned int side=0; side<2; side++) {

	sfOFF[side]=0.;
	sfON[side]=0.;
	sfOFF_err[side]=0.;
	sfON_err[side]=0.;
	double weightsum=0.;
	
	// AFP Off
	sfOFF[side] = anaChoice_A10[side][bin];
	weightsum = (anaChoice_A10err[side][bin]>0.?power(1./anaChoice_A10err[side][bin],2):0.);
	
	sfOFF_err[side] = weightsum>0. ? 1./sqrt(weightsum) : 0.;
	
	weightsum=0.;
	
	// AFP ON
	sfON[side] = anaChoice_A7[side][bin];
	weightsum = (anaChoice_A7err[side][bin]>0.?power(1./anaChoice_A7err[side][bin],2):0.);
	
	sfON_err[side] = weightsum>0. ? 1./sqrt(weightsum) : 0.;
	
      }
      if (sfOFF[0]>0. && sfOFF[1]>0. && sfON[0]>0. && sfON[1]>0.) { 
	R = sfOFF[1]*sfON[0]/(sfON[1]*sfOFF[0]);
	deltaR = sqrt(R*R*(power(sfOFF_err[0]/sfOFF[0],2)+power(sfON_err[1]/sfON[1],2)+power(sfOFF_err[1]/sfOFF[1],2)+power(sfON_err[0]/sfON[0],2)));
	asymmetry[0][1][bin] = (1.-sqrt(R))/(1+sqrt(R));
	asymmetryError[0][1][bin] = (deltaR)/(sqrt(std::abs(R))*power((sqrt(std::abs(R))+1.),2));
      }
      else {
	asymmetry[0][1][bin] = 0.;
	asymmetryError[0][1][bin] = 0.;
      }
    }
  }

  if (isGoodPair[1][0]) {
    for (unsigned int bin=0; bin<asymmetry[1][0].size(); bin++) {
      double R = 0.;
      double deltaR = 0.;
      
      for (unsigned int side=0; side<2; side++) {
	sfOFF[side]=0.;
	sfON[side]=0.;
	sfOFF_err[side]=0.;
	sfON_err[side]=0.;
	double weightsum=0.;
	
	// AFP Off
	sfOFF[side] = anaChoice_B5[side][bin];
	weightsum = (anaChoice_B5err[side][bin]>0.?power(1./anaChoice_B5err[side][bin],2):0.);
	
	sfOFF_err[side] = weightsum>0. ? 1./sqrt(weightsum) : 0.;
	
	weightsum=0.;
	
	// AFP ON
	sfON[side] = anaChoice_B2[side][bin];
	weightsum = (anaChoice_B2err[side][bin]>0.?power(1./anaChoice_B2err[side][bin],2):0.);
	
	sfON_err[side] = weightsum>0. ? 1./sqrt(weightsum) : 0.;
	
      }
      if (sfOFF[0]>0. && sfOFF[1]>0. && sfON[0]>0. && sfON[1]>0.) { 
	R = sfOFF[1]*sfON[0]/(sfON[1]*sfOFF[0]);
	deltaR = sqrt(R*R*(power(sfOFF_err[0]/sfOFF[0],2)+power(sfON_err[1]/sfON[1],2)+power(sfOFF_err[1]/sfOFF[1],2)+power(sfON_err[0]/sfON[0],2)));
	asymmetry[1][0][bin] = (1.-sqrt(R))/(1+sqrt(R));
	asymmetryError[1][0][bin] = (deltaR)/(sqrt(std::abs(R))*power((sqrt(std::abs(R))+1.),2));
      }
      else {
	asymmetry[1][0][bin] = 0.;
	asymmetryError[1][0][bin] = 0.;
      }
    }
  }

  if (isGoodPair[1][1]) {
    for (unsigned int bin=0; bin<asymmetry[1][1].size(); bin++) {
      double R = 0.;
      double deltaR = 0.;
      
      for (unsigned int side=0; side<2; side++) {
	sfOFF[side]=0.;
	sfON[side]=0.;
	sfOFF_err[side]=0.;
	sfON_err[side]=0.;
	double weightsum=0.;
	
	// AFP Off
	sfOFF[side] = anaChoice_B7[side][bin];
	weightsum = (anaChoice_B7err[side][bin]>0.?power(1./anaChoice_B7err[side][bin],2):0.);
	
	sfOFF_err[side] = weightsum>0. ? 1./sqrt(weightsum) : 0.;
	
	weightsum=0.;
	
	// AFP ON
	sfON[side] = anaChoice_B10[side][bin];
	weightsum = (anaChoice_B10err[side][bin]>0.?power(1./anaChoice_B10err[side][bin],2):0.);
	
	sfON_err[side] = weightsum>0. ? 1./sqrt(weightsum) : 0.;
      }

      //std::cout << bin << " " << sfON[0] << " " << sfON[1] << " "  << sfOFF[0] << " " << sfOFF[1] << std::endl;

      if (sfOFF[0]>0. && sfOFF[1]>0. && sfON[0]>0. && sfON[1]>0.) { 
	R = sfOFF[1]*sfON[0]/(sfON[1]*sfOFF[0]);
	deltaR = sqrt(R*R*(power(sfOFF_err[0]/sfOFF[0],2)+power(sfON_err[1]/sfON[1],2)+power(sfOFF_err[1]/sfOFF[1],2)+power(sfON_err[0]/sfON[0],2)));
	asymmetry[1][1][bin] = (1.-sqrt(R))/(1+sqrt(R));
	asymmetryError[1][1][bin] = (deltaR)/(sqrt(std::abs(R))*power((sqrt(std::abs(R))+1.),2));
      }
      else {
	asymmetry[1][1][bin] = 0.;
	asymmetryError[1][1][bin] = 0.;
      }
    }
  }
  
  boolAsymmetry = true;
};

void PairAsymmetry::calcTotalAsymmetry(double enWinLow, double enWinHigh, int anaChoice) {
  if (!isAnaChoiceRateVectors() || getCurrentAnaChoice()!=anaChoice) {
    makeAnalysisChoiceRateVectors(anaChoice);
  }
  unsigned int binLow = (unsigned int)(enWinLow/energyBinWidth);
  unsigned int binHigh = (unsigned int)(enWinHigh/energyBinWidth)-1;

  double sfON[2]={0.}, sfOFF[2]={0.};
  double sfON_err[2]={0.}, sfOFF_err[2]={0.};
  double sumA2[2]={0.}, sumA5[2]={0.}, sumA7[2]={0.}, sumA10[2]={0.}, sumB2[2]={0.}, sumB5[2]={0.}, sumB7[2]={0.}, sumB10[2]={0.};
  double sumA2_err[2]={0.}, sumA5_err[2]={0.}, sumA7_err[2]={0.}, sumA10_err[2]={0.}, sumB2_err[2]={0.}, sumB5_err[2]={0.}, sumB7_err[2]={0.}, sumB10_err[2]={0.};

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

    }
    sumA2_err[side] = sqrt(sumA2_err[side]);
    sumA5_err[side] = sqrt(sumA5_err[side]);
    sumA7_err[side] = sqrt(sumA7_err[side]);
    sumA10_err[side] = sqrt(sumA10_err[side]);
    sumB2_err[side] = sqrt(sumB2_err[side]);
    sumB5_err[side] = sqrt(sumB5_err[side]);
    sumB7_err[side] = sqrt(sumB7_err[side]);
    sumB10_err[side] = sqrt(sumB10_err[side]);
    
  }
    
  
  // A type runs, pair 0
  if (isGoodPair[0][0]) {
    
    double R = 0.;
    double deltaR = 0.;
    
    for (unsigned int side=0; side<2; side++) {
      double weightsum=0.;
      
      // AFP Off
      sfOFF[side] = sumA2[side];
      weightsum = (sumA2_err[side]>0.?power(1./sumA2_err[side],2):0.);
      
      sfOFF_err[side] = weightsum>0. ? 1./sqrt(weightsum) : 0.;
      
      weightsum=0.;
      
      // AFP ON
      sfON[side] = sumA5[side];
	weightsum = (sumA5_err[side]>0.?power(1./sumA5_err[side],2):0.);
	
	sfON_err[side] = weightsum>0. ? 1./sqrt(weightsum) : 0.;
	
    }
    
    if (sfOFF[0]>0. && sfOFF[1]>0. && sfON[0]>0. && sfON[1]>0.) { 
      R = sfOFF[1]*sfON[0]/(sfON[1]*sfOFF[0]);
      deltaR = sqrt(R*R*(power(sfOFF_err[0]/sfOFF[0],2)+power(sfON_err[1]/sfON[1],2)+power(sfOFF_err[1]/sfOFF[1],2)+power(sfON_err[0]/sfON[0],2)));
      totalAsymmetryA0= (1.-sqrt(R))/(1+sqrt(R));
      totalAsymmetryErrorA0 = (deltaR)/(sqrt(std::abs(R))*power((sqrt(std::abs(R))+1.),2));
    }
    else {
      totalAsymmetryA0 = 0.;
      totalAsymmetryErrorA0= 0.;
    }
  }
  
  
  if (isGoodPair[0][1]) {
    double R = 0.;
    double deltaR = 0.;
    
    for (unsigned int side=0; side<2; side++) {
      sfOFF[side]=0.;
      sfON[side]=0.;
      sfOFF_err[side]=0.;
      sfON_err[side]=0.;
      double weightsum=0.;
      
      // AFP Off
      sfOFF[side] = sumA10[side];
	weightsum = (sumA10_err[side]>0.?power(1./sumA10_err[side],2):0.);
	
	sfOFF_err[side] = weightsum>0. ? 1./sqrt(weightsum) : 0.;
	
	weightsum=0.;
	
	// AFP ON
	sfON[side] = sumA7[side];
	weightsum = (sumA7_err[side]>0.?power(1./sumA7_err[side],2):0.);
	
	sfON_err[side] = weightsum>0. ? 1./sqrt(weightsum) : 0.;
	
    }
    if (sfOFF[0]>0. && sfOFF[1]>0. && sfON[0]>0. && sfON[1]>0.) { 
      R = sfOFF[1]*sfON[0]/(sfON[1]*sfOFF[0]);
      deltaR = sqrt(R*R*(power(sfOFF_err[0]/sfOFF[0],2)+power(sfON_err[1]/sfON[1],2)+power(sfOFF_err[1]/sfOFF[1],2)+power(sfON_err[0]/sfON[0],2)));
      totalAsymmetryA1= (1.-sqrt(R))/(1+sqrt(R));
      totalAsymmetryErrorA1 = (deltaR)/(sqrt(std::abs(R))*power((sqrt(std::abs(R))+1.),2));
    }
    else {
      totalAsymmetryA1 = 0.;
      totalAsymmetryErrorA1= 0.;
    }
  }

  
  if (isGoodPair[1][0]) {
    double R = 0.;
    double deltaR = 0.;
    
    for (unsigned int side=0; side<2; side++) {
      sfOFF[side]=0.;
      sfON[side]=0.;
      sfOFF_err[side]=0.;
      sfON_err[side]=0.;
      double weightsum=0.;
      
      // AFP Off
      sfOFF[side] = sumB5[side];
      weightsum = (sumB5_err[side]>0.?power(1./sumB5_err[side],2):0.);
      
      sfOFF_err[side] = weightsum>0. ? 1./sqrt(weightsum) : 0.;
      
      weightsum=0.;
      
      // AFP ON
      sfON[side] = sumB2[side];
      weightsum = (sumB2_err[side]>0.?power(1./sumB2_err[side],2):0.);
      
      sfON_err[side] = weightsum>0. ? 1./sqrt(weightsum) : 0.;
      
    }

    if (sfOFF[0]>0. && sfOFF[1]>0. && sfON[0]>0. && sfON[1]>0.) { 
      R = sfOFF[1]*sfON[0]/(sfON[1]*sfOFF[0]);
      deltaR = sqrt(R*R*(power(sfOFF_err[0]/sfOFF[0],2)+power(sfON_err[1]/sfON[1],2)+power(sfOFF_err[1]/sfOFF[1],2)+power(sfON_err[0]/sfON[0],2)));
      totalAsymmetryB0= (1.-sqrt(R))/(1+sqrt(R));
      totalAsymmetryErrorB0 = (deltaR)/(sqrt(std::abs(R))*power((sqrt(std::abs(R))+1.),2));
    }
    else {
      totalAsymmetryB0 = 0.;
      totalAsymmetryErrorB0= 0.;
    }

  }


  if (isGoodPair[1][1]) {
    double R = 0.;
    double deltaR = 0.;
    
    for (unsigned int side=0; side<2; side++) {
      sfOFF[side]=0.;
      sfON[side]=0.;
      sfOFF_err[side]=0.;
      sfON_err[side]=0.;
      double weightsum=0.;
      
      // AFP Off
      sfOFF[side] = sumB7[side];
      weightsum = (sumB7_err[side]>0.?power(1./sumB7_err[side],2):0.);
      
      sfOFF_err[side] = weightsum>0. ? 1./sqrt(weightsum) : 0.;
      
      weightsum=0.;
      
      // AFP ON
      sfON[side] = sumB10[side];
      weightsum = (sumB10_err[side]>0.?power(1./sumB10_err[side],2):0.);
      
      sfON_err[side] = weightsum>0. ? 1./sqrt(weightsum) : 0.;
      
    }

    if (sfOFF[0]>0. && sfOFF[1]>0. && sfON[0]>0. && sfON[1]>0.) { 
      R = sfOFF[1]*sfON[0]/(sfON[1]*sfOFF[0]);
      deltaR = sqrt(R*R*(power(sfOFF_err[0]/sfOFF[0],2)+power(sfON_err[1]/sfON[1],2)+power(sfOFF_err[1]/sfOFF[1],2)+power(sfON_err[0]/sfON[0],2)));
      totalAsymmetryB1= (1.-sqrt(R))/(1+sqrt(R));
      totalAsymmetryErrorB1 = (deltaR)/(sqrt(std::abs(R))*power((sqrt(std::abs(R))+1.),2));
    }
    else {
      totalAsymmetryB1 = 0.;
      totalAsymmetryErrorB1= 0.;
    }
    

  }
};


void PairAsymmetry::calcSuperSum(int anaChoice) {
  if (!isAnaChoiceRateVectors() || getCurrentAnaChoice()!=anaChoice) {
    makeAnalysisChoiceRateVectors(anaChoice);
  }
  //unsigned int numBins = (unsigned int)(1200./energyBinWidth);
  superSum.resize(2, std::vector < std::vector <double> > (2,std::vector <double> (numEnBins, 0.)));
  superSumError.resize(2, std::vector < std::vector <double> > (2,std::vector <double> (numEnBins, 0.)));;

  double sfON[2]={0.};
  double sfOFF[2]={0.};
  double sfON_err[2]={0.};
  double sfOFF_err[2]={0.};

  

  // A type runs, pair 0
  if (isGoodPair[0][0]) {
    for (unsigned int bin=0; bin<asymmetry[0][0].size(); bin++) {     
      for (unsigned int side=0; side<2; side++) {
	double weightsum=0.;
	
	// AFP Off
	sfOFF[side] = anaChoice_A2[side][bin];
	weightsum = (anaChoice_A2err[side][bin]>0.?power(1./anaChoice_A2err[side][bin],2):0.);
	
	sfOFF_err[side] = weightsum>0. ? 1./sqrt(weightsum) : 0.;
	
	weightsum=0.;
	
	// AFP ON
	sfON[side] = anaChoice_A5[side][bin];
	weightsum = (anaChoice_A5err[side][bin]>0.?power(1./anaChoice_A5err[side][bin],2):0.);
	
	sfON_err[side] = weightsum>0. ? 1./sqrt(weightsum) : 0.;
	
      }
      //The geometric mean as defined by http://www.arpapress.com/Volumes/Vol11Issue3/IJRRAS_11_3_08.pdf
      double R1 = (sfON[1]>0. && sfOFF[0]>0.) ? sqrt(sfOFF[0]*sfON[1]) : (sfON[1]<0. && sfOFF[0]<0.) ? -sqrt(sfOFF[0]*sfON[1]) : 0.5*(sfOFF[0]+sfON[1]); 
      double R2 = (sfON[0]>0. && sfOFF[1]>0.) ? sqrt(sfOFF[1]*sfON[0]) : (sfON[0]<0. && sfOFF[1]<0.) ? -sqrt(sfOFF[1]*sfON[0]) : 0.5*(sfOFF[1]+sfON[0]);
      double deltaR1 = ((sfON[1]>0. && sfOFF[0]>0.) || (sfON[1]<0. && sfOFF[0]<0.)) ? 0.5*sqrt((power(sfOFF[0]*sfON_err[1],2)+power(sfON[1]*sfOFF_err[0],2))/(sfON[1]*sfOFF[0])) : 
	0.5*sqrt(power(sfON_err[1],2) + power(sfOFF_err[0],2));
      double deltaR2 = ((sfON[0]>0. && sfOFF[1]>0.) || (sfON[0]<0. && sfOFF[1]<0.)) ? 0.5*sqrt((power(sfOFF[1]*sfON_err[0],2)+power(sfON[0]*sfOFF_err[1],2))/(sfON[0]*sfOFF[1])) : 
	0.5*sqrt(power(sfON_err[0],2) + power(sfOFF_err[1],2));
      
      superSum[0][0][bin] = R1 + R2;
      //std::cout << sfOFF[0] << " " << sfOFF_err[0] << std::endl;
      superSumError[0][0][bin] = sqrt(power(deltaR1,2) + power(deltaR2,2));
    }
  }

  if (isGoodPair[0][1]) {
    for (unsigned int bin=0; bin<asymmetry[0][1].size(); bin++) {
      for (unsigned int side=0; side<2; side++) {
	sfOFF[side]=0.;
	sfON[side]=0.;
	sfOFF_err[side]=0.;
	sfON_err[side]=0.;
	double weightsum=0.;
	
	// AFP Off
	sfOFF[side] = anaChoice_A10[side][bin];
	weightsum = (anaChoice_A10err[side][bin]>0.?power(1./anaChoice_A10err[side][bin],2):0.);
	
	sfOFF_err[side] = weightsum>0. ? 1./sqrt(weightsum) : 0.;
	
	weightsum=0.;
	
	// AFP ON
	sfON[side] = anaChoice_A7[side][bin];
	weightsum = (anaChoice_A7err[side][bin]>0.?power(1./anaChoice_A7err[side][bin],2):0.);
	
	sfON_err[side] = weightsum>0. ? 1./sqrt(weightsum) : 0.;
	
      }
      //The geometric mean as defined by http://www.arpapress.com/Volumes/Vol11Issue3/IJRRAS_11_3_08.pdf
      double R1 = (sfON[1]>0. && sfOFF[0]>0.) ? sqrt(sfOFF[0]*sfON[1]) : (sfON[1]<0. && sfOFF[0]<0.) ? -sqrt(sfOFF[0]*sfON[1]) : 0.5*(sfOFF[0]+sfON[1]); 
      double R2 = (sfON[0]>0. && sfOFF[1]>0.) ? sqrt(sfOFF[1]*sfON[0]) : (sfON[0]<0. && sfOFF[1]<0.) ? -sqrt(sfOFF[1]*sfON[0]) : 0.5*(sfOFF[1]+sfON[0]);
      double deltaR1 = ((sfON[1]>0. && sfOFF[0]>0.) || (sfON[1]<0. && sfOFF[0]<0.)) ? 0.5*sqrt((power(sfOFF[0]*sfON_err[1],2)+power(sfON[1]*sfOFF_err[0],2))/(sfON[1]*sfOFF[0])) : 
	0.5*sqrt(power(sfON_err[1],2) + power(sfOFF_err[0],2));
      double deltaR2 = ((sfON[0]>0. && sfOFF[1]>0.) || (sfON[0]<0. && sfOFF[1]<0.)) ? 0.5*sqrt((power(sfOFF[1]*sfON_err[0],2)+power(sfON[0]*sfOFF_err[1],2))/(sfON[0]*sfOFF[1])) : 
	0.5*sqrt(power(sfON_err[0],2) + power(sfOFF_err[1],2));
      
      superSum[0][1][bin] = R1 + R2;
      //std::cout << sfOFF[0] << " " << sfOFF_err[0] << std::endl;
      superSumError[0][1][bin] = sqrt(power(deltaR1,2) + power(deltaR2,2));
      
    }
  }

  if (isGoodPair[1][0]) {
    for (unsigned int bin=0; bin<asymmetry[1][0].size(); bin++) {
      for (unsigned int side=0; side<2; side++) {
	sfOFF[side]=0.;
	sfON[side]=0.;
	sfOFF_err[side]=0.;
	sfON_err[side]=0.;
	double weightsum=0.;
	
	// AFP Off
	sfOFF[side] = anaChoice_B5[side][bin];
	weightsum = (anaChoice_B5err[side][bin]>0.?power(1./anaChoice_B5err[side][bin],2):0.);
	
	sfOFF_err[side] = weightsum>0. ? 1./sqrt(weightsum) : 0.;
	
	weightsum=0.;
	
	// AFP ON
	sfON[side] = anaChoice_B2[side][bin];
	weightsum = (anaChoice_B2err[side][bin]>0.?power(1./anaChoice_B2err[side][bin],2):0.);
	
	sfON_err[side] = weightsum>0. ? 1./sqrt(weightsum) : 0.;
	
      }
      //The geometric mean as defined by http://www.arpapress.com/Volumes/Vol11Issue3/IJRRAS_11_3_08.pdf
      double R1 = (sfON[1]>0. && sfOFF[0]>0.) ? sqrt(sfOFF[0]*sfON[1]) : (sfON[1]<0. && sfOFF[0]<0.) ? -sqrt(sfOFF[0]*sfON[1]) : 0.5*(sfOFF[0]+sfON[1]); 
      double R2 = (sfON[0]>0. && sfOFF[1]>0.) ? sqrt(sfOFF[1]*sfON[0]) : (sfON[0]<0. && sfOFF[1]<0.) ? -sqrt(sfOFF[1]*sfON[0]) : 0.5*(sfOFF[1]+sfON[0]);
      double deltaR1 = ((sfON[1]>0. && sfOFF[0]>0.) || (sfON[1]<0. && sfOFF[0]<0.)) ? 0.5*sqrt((power(sfOFF[0]*sfON_err[1],2)+power(sfON[1]*sfOFF_err[0],2))/(sfON[1]*sfOFF[0])) : 
	0.5*sqrt(power(sfON_err[1],2) + power(sfOFF_err[0],2));
      double deltaR2 = ((sfON[0]>0. && sfOFF[1]>0.) || (sfON[0]<0. && sfOFF[1]<0.)) ? 0.5*sqrt((power(sfOFF[1]*sfON_err[0],2)+power(sfON[0]*sfOFF_err[1],2))/(sfON[0]*sfOFF[1])) : 
	0.5*sqrt(power(sfON_err[0],2) + power(sfOFF_err[1],2));
      
      superSum[1][0][bin] = R1 + R2;
      //std::cout << sfOFF[0] << " " << sfOFF_err[0] << std::endl;
      superSumError[1][0][bin] = sqrt(power(deltaR1,2) + power(deltaR2,2));
    }
  }

  if (isGoodPair[1][1]) {
    for (unsigned int bin=0; bin<asymmetry[1][1].size(); bin++) {
      
      for (unsigned int side=0; side<2; side++) {
	sfOFF[side]=0.;
	sfON[side]=0.;
	sfOFF_err[side]=0.;
	sfON_err[side]=0.;
	double weightsum=0.;
	
	// AFP Off
	sfOFF[side] = anaChoice_B7[side][bin];
	weightsum = (anaChoice_B7err[side][bin]>0.?power(1./anaChoice_B7err[side][bin],2):0.);
	
	sfOFF_err[side] = weightsum>0. ? 1./sqrt(weightsum) : 0.;
	
	weightsum=0.;
	
	// AFP ON
	sfON[side] = anaChoice_B10[side][bin];
	weightsum = (anaChoice_B10err[side][bin]>0.?power(1./anaChoice_B10err[side][bin],2):0.);
	
	sfON_err[side] = weightsum>0. ? 1./sqrt(weightsum) : 0.;
	
      }
      //The geometric mean as defined by http://www.arpapress.com/Volumes/Vol11Issue3/IJRRAS_11_3_08.pdf
      double R1 = (sfON[1]>0. && sfOFF[0]>0.) ? sqrt(sfOFF[0]*sfON[1]) : (sfON[1]<0. && sfOFF[0]<0.) ? -sqrt(sfOFF[0]*sfON[1]) : 0.5*(sfOFF[0]+sfON[1]); 
      double R2 = (sfON[0]>0. && sfOFF[1]>0.) ? sqrt(sfOFF[1]*sfON[0]) : (sfON[0]<0. && sfOFF[1]<0.) ? -sqrt(sfOFF[1]*sfON[0]) : 0.5*(sfOFF[1]+sfON[0]);
      double deltaR1 = ((sfON[1]>0. && sfOFF[0]>0.) || (sfON[1]<0. && sfOFF[0]<0.)) ? 0.5*sqrt((power(sfOFF[0]*sfON_err[1],2)+power(sfON[1]*sfOFF_err[0],2))/(sfON[1]*sfOFF[0])) : 
	0.5*sqrt(power(sfON_err[1],2) + power(sfOFF_err[0],2));
      double deltaR2 = ((sfON[0]>0. && sfOFF[1]>0.) || (sfON[0]<0. && sfOFF[1]<0.)) ? 0.5*sqrt((power(sfOFF[1]*sfON_err[0],2)+power(sfON[0]*sfOFF_err[1],2))/(sfON[0]*sfOFF[1])) : 
	0.5*sqrt(power(sfON_err[0],2) + power(sfOFF_err[1],2));
      
      superSum[1][1][bin] = R1 + R2;
      //std::cout << sfOFF[0] << " " << sfOFF_err[0] << std::endl;
      superSumError[1][1][bin] = sqrt(power(deltaR1,2) + power(deltaR2,2));
    }
  }  


  boolSuperSum = true;
};

void PairAsymmetry::writeAsymToFile() {
  if (!isAsymmetry()) {
    std::cout << "No Asymmetry has been calculated. Not writing to file!\n\n";
    return;
  }

  int anaChoice = getCurrentAnaChoice();
  //Setting paths to output files
  std::string outpathA0 = (Simulation ? std::string(getenv("SIM_ANALYSIS_RESULTS")) : UKdata ? std::string(getenv("ANALYSIS_RESULTS")) : std::string(getenv("MPM_ANALYSIS_RESULTS"))) + 
    "Octet_" + itos(octet) + "/PairAsymmetry/"+ (UNBLIND?"UNBLINDED_":"")+"rawAsymmetry_Octet" + itos(octet)+"_AnaCh"+itos(anaChoice)+"_Pair_A0.dat";
  std::string outpathA1 = (Simulation ? std::string(getenv("SIM_ANALYSIS_RESULTS")) : UKdata ? std::string(getenv("ANALYSIS_RESULTS")) : std::string(getenv("MPM_ANALYSIS_RESULTS"))) + 
    "Octet_" + itos(octet) + "/PairAsymmetry/"+ (UNBLIND?"UNBLINDED_":"")+"rawAsymmetry_Octet" + itos(octet)+"_AnaCh"+itos(anaChoice)+"_Pair_A1.dat";
  std::string outpathB0 = (Simulation ? std::string(getenv("SIM_ANALYSIS_RESULTS")) : UKdata ? std::string(getenv("ANALYSIS_RESULTS")) : std::string(getenv("MPM_ANALYSIS_RESULTS"))) +
    "Octet_" + itos(octet) + "/PairAsymmetry/"+ (UNBLIND?"UNBLINDED_":"")+"rawAsymmetry_Octet" + itos(octet)+"_AnaCh"+itos(anaChoice)+"_Pair_B0.dat";
  std::string outpathB1 = (Simulation ? std::string(getenv("SIM_ANALYSIS_RESULTS")) : UKdata ? std::string(getenv("ANALYSIS_RESULTS")) : std::string(getenv("MPM_ANALYSIS_RESULTS"))) +
    "Octet_" + itos(octet) + "/PairAsymmetry/"+ (UNBLIND?"UNBLINDED_":"")+"rawAsymmetry_Octet" + itos(octet)+"_AnaCh"+itos(anaChoice)+"_Pair_B1.dat";

  //Open and fill file for quartet A
  std::ofstream outfileA0(outpathA0.c_str()); 
  
  if (isGoodPair[0][0]) {
    for (unsigned int i=0; i<asymmetry[0][0].size(); i++) {
      outfileA0 << binLowerEdge[i] << " " << asymmetry[0][0][i] << " " << asymmetryError[0][0][i] << std::endl;
      //std::cout << binLowerEdge[i] << " " << asymmetry[i] << " " << asymmetryError[i] << std::endl;
    }
  }
  else outfileA0 << "BAD PAIR";
  outfileA0.close();

  std::ofstream outfileA1(outpathA1.c_str()); 
  
  if (isGoodPair[0][1]) {
    for (unsigned int i=0; i<asymmetry[0][1].size(); i++) {
      outfileA1 << binLowerEdge[i] << " " << asymmetry[0][1][i] << " " << asymmetryError[0][1][i] << std::endl;
      //std::cout << binLowerEdge[i] << " " << asymmetry[i] << " " << asymmetryError[i] << std::endl;
    }
  }
  else outfileA1 << "BAD PAIR";
  outfileA1.close();

  std::ofstream outfileB0(outpathB0.c_str()); 

  if (isGoodPair[1][0]) {
    for (unsigned int i=0; i<asymmetry[1][0].size(); i++) {
      outfileB0 << binLowerEdge[i] << " " << asymmetry[1][0][i] << " " << asymmetryError[1][0][i] << std::endl;
      //std::cout << binLowerEdge[i] << " " << asymmetry[i] << " " << asymmetryError[i] << std::endl;
    }
  }
  else outfileB0 << "BAD PAIR";
  outfileB0.close();

  std::ofstream outfileB1(outpathB1.c_str()); 
  
  if (isGoodPair[1][1]) {
    for (unsigned int i=0; i<asymmetry[1][1].size(); i++) {
      outfileB1 << binLowerEdge[i] << " " << asymmetry[1][1][i] << " " << asymmetryError[1][1][i] << std::endl;
      //std::cout << binLowerEdge[i] << " " << asymmetry[1][1][i] << " " << asymmetryError[1][1][i] << std::endl;
    }
  }
  else outfileB1 << "BAD PAIR";
  outfileB1.close();

  
  std::cout << "Wrote Asymmetry to file for anaChoice " << anaChoice << "\n";
};

void PairAsymmetry::writeSuperSumToFile() {
  if (!isSuperSum()) {
    std::cout << "No Super-Sum has been calculated. Not writing to file!\n\n";
    return;
  }

  int anaChoice = getCurrentAnaChoice();
  //Setting paths to output files
  std::string outpathA0 = (Simulation ? std::string(getenv("SIM_ANALYSIS_RESULTS")) : UKdata ? std::string(getenv("ANALYSIS_RESULTS")) : std::string(getenv("MPM_ANALYSIS_RESULTS"))) + 
    "Octet_" + itos(octet) + "/PairAsymmetry/"+ (UNBLIND?"UNBLINDED_":"")+"superSum_Octet" + itos(octet)+"_AnaCh"+itos(anaChoice)+"_Pair_A0.dat";
  std::string outpathA1 = (Simulation ? std::string(getenv("SIM_ANALYSIS_RESULTS")) : UKdata ? std::string(getenv("ANALYSIS_RESULTS")) : std::string(getenv("MPM_ANALYSIS_RESULTS"))) + 
    "Octet_" + itos(octet) + "/PairAsymmetry/"+ (UNBLIND?"UNBLINDED_":"")+"superSum_Octet" + itos(octet)+"_AnaCh"+itos(anaChoice)+"_Pair_A1.dat";
  std::string outpathB0 = (Simulation ? std::string(getenv("SIM_ANALYSIS_RESULTS")) : UKdata ? std::string(getenv("ANALYSIS_RESULTS")) : std::string(getenv("MPM_ANALYSIS_RESULTS"))) +
    "Octet_" + itos(octet) + "/PairAsymmetry/"+ (UNBLIND?"UNBLINDED_":"")+"superSum_Octet" + itos(octet)+"_AnaCh"+itos(anaChoice)+"_Pair_B0.dat";
  std::string outpathB1 = (Simulation ? std::string(getenv("SIM_ANALYSIS_RESULTS")) : UKdata ? std::string(getenv("ANALYSIS_RESULTS")) : std::string(getenv("MPM_ANALYSIS_RESULTS"))) +
    "Octet_" + itos(octet) + "/PairAsymmetry/"+ (UNBLIND?"UNBLINDED_":"")+"superSum_Octet" + itos(octet)+"_AnaCh"+itos(anaChoice)+"_Pair_B1.dat";
  
  //Open and fill file for quartet A
  std::ofstream outfileA0(outpathA0.c_str()); 
  
  if (isGoodPair[0][0]) {
    for (unsigned int i=0; i<superSum[0][0].size(); i++) {
      outfileA0 << binLowerEdge[i] << " " << superSum[0][0][i] << " " << superSumError[0][0][i] << std::endl;
      //std::cout << binLowerEdge[i] << " " << superSum[0][0][i] << " " << superSumError[0][0][i] << std::endl;
    }
  }
  else outfileA0 << "BAD PAIR";
  outfileA0.close();
  
  std::ofstream outfileA1(outpathA1.c_str()); 
  
  if (isGoodPair[0][1]) {
    for (unsigned int i=0; i<superSum[0][1].size(); i++) {
      outfileA1 << binLowerEdge[i] << " " << superSum[0][1][i] << " " << superSumError[0][1][i] << std::endl;
    }
  }
  else outfileA1 << "BAD PAIR";
  outfileA1.close();
  
  std::ofstream outfileB0(outpathB0.c_str()); 
  
  if (isGoodPair[1][0]) {
    for (unsigned int i=0; i<superSum[1][0].size(); i++) {
      outfileB0 << binLowerEdge[i] << " " << superSum[1][0][i] << " " << superSumError[1][0][i] << std::endl;
    }
  }
  else outfileB0 << "BAD PAIR";
  outfileB0.close();
  
  std::ofstream outfileB1(outpathB1.c_str()); 
  
  if (isGoodPair[1][1]) {
    for (unsigned int i=0; i<superSum[1][1].size(); i++) {
      outfileB1 << binLowerEdge[i] << " " << superSum[1][1][i] << " " << superSumError[1][1][i] << std::endl;
    }
  }
  else outfileB1 << "BAD PAIR";
  outfileB1.close();
  
  
  std::cout << "Wrote Super-Sum to file for anaChoice " << anaChoice << "\n";
  //std::cout << "Wrote to " << outpathB1 << "\n";
};

