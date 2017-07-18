/*
Compilable code which houses the heart of the analysis. 
Here are classes which handle calculating asymmetries of either data 
or simulation. There will be flags to be passed dictating what type of 
asymmetry, what type of data, and what octet/runs to use. 

The code should be able to do BG subtraction, construct asymmetries, 
or do both depending on the flags passed.

*/


#include "MBUtils.hh"
#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <cmath>
#include <fstream>
#include <algorithm>
#include <map>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <TStyle.h>

#include "BetaSpectrum.hh"

TString anaChoices[10] = {"A","B","C","D","E","F","G","H","J","K"};


//Types of Corrections to apply
std::string corr ("AllCorr");//{"UnCorr","DeltaExpOnly","DeltaTheoryOnly","AllCorr"};
                             

Double_t POL_minus = 0.9981;
Double_t POL_plus = 0.9937;
Double_t delta_POL = POL_plus-POL_minus;
Double_t POL_ave = (POL_plus+POL_minus) / 2.;

bool withPOL = false; //Set this to true to correct DATA for the polarimetry measurement


std::vector <Int_t> badOct {7,60,61,62,63,64,65,66}; // 9,59 used to be part of this, but not totally sure why... they are just not complete octets
                                  // Octet 7 had W anode dead for part of run
                                  // Either need to discard, or apply the 
                                  // Charge cloud method to determine a 
                                  // coincidence
                                  // 60,61,62,63,64,65,66 have bad TDC spectra
                                  // 70 and 92 just have low statistics so they had bad corrections. Put them back in when using high statistics corrections

std::vector <Double_t> LoadFieldDipCorrections(std::vector <Double_t> enBinMidpoint,Int_t oct);

//Returns the values of the systematic corrections and statistical error on this for now
std::vector < std::vector <Double_t> >  LoadOctetSystematics(Int_t octet, std::string anaChoice, std::vector <Double_t> enBinMidpoint);

//Returns a vector containing all the theory corrections to A0 for a particular bin
std::vector <Double_t> LoadTheoryCorrections(std::vector <Double_t> enBinMidpoint);

// Collects all the asymmetries, bin-by-bin, and produces a final asymmetry plot, both A_SR and 2*A/Beta, for whatever grouping provided ("Octet", "Quartet", "Pair"). Also makes a plot of the 
// integrated asymmetry vs octet/quartet/pair number
void PlotFinalAsymmetries(std::string groupType, Int_t octBegin, Int_t octEnd, std::string anaChoice, Double_t Elow=220., Double_t Ehigh=680., Double_t enBinWidth=10., bool UKdata=true, bool simulation=false, bool UNBLIND=false, int key=0);


// Function to return beta when given the kinetic energy of an electron
Double_t returnBeta(Double_t En) { 
  Double_t me = 510.998928; //rest mass energy of electron in keV
  return sqrt(En*En+2.*En*me)/(En+me);
};



int main(int argc, char* argv[])
{
  
 
  if (argc!=8) {
    std::cout << "USAGE: ./MC_SpectralCorrection.exe [anaCh] [octmin] [octmax] [analysis window min] [analysis window max] [corr] [sim]  \n";
    exit(0);
  }
    
  std::string analysisChoice = argv[1];
  Int_t octBegin = atoi(argv[2]);
  Int_t octEnd = atoi(argv[3]);
  Double_t enBinWidth = 10.;
  Double_t Elow = atoi(argv[4]);
  Double_t Ehigh = atoi(argv[5]);//680
  std::string sim = argv[7];
  corr = std::string(argv[6]);
  int key = 0;
 
  bool UKdata = true;//true;
  bool simulation = sim==std::string("true") ? true : false;

  if (simulation) withPOL=false;

  //I should keep track of the raw asymmetry, Experimental systematic corrected asymmetry, and theoretical 
  // Systematic asymmetry in a file each time I run to see the percent correction of each...

  //NEED TO ADD IN ABILITY TO READ IN SYSTEMATICS FOR QUARTET AND PAIR

  //****************************************************************
  //****************************************************************
  // ONLY TURN THIS ON WHEN READY TO UNBLIND. IT WILL USE THE TRUE 
  // CLOCK TIMES.
  //****************************************************************
  //****************************************************************
  bool UNBLIND = false;//unblind==std::string("true") ? true : false;

  
  try {
    PlotFinalAsymmetries("Octet",octBegin, octEnd, analysisChoice, Elow, Ehigh, enBinWidth, UKdata, simulation, UNBLIND, key);
  }
  catch(const char* ex){
    std::cerr << "Error: " << ex << std::endl;
  }
  
  
  return 0;
}





void PlotFinalAsymmetries(std::string groupType, Int_t octBegin, Int_t octEnd, std::string anaChoice, Double_t Elow, Double_t Ehigh, 
			  Double_t enBinWidth, bool UKdata, bool simulation, bool UNBLIND, int key) {

  if (groupType!="Quartet" && groupType!="Octet" && groupType!="Pair") throw "Bad group type given to PlotFinalAsymmetries. Options are \"Octet\", \"Quartet\", or \"Pair\""; 

  key = ( groupType==std::string("Octet") ? key : 0 );

  Int_t numBins = 1200/enBinWidth;

  std::vector < Double_t > enBinMedian; //Holds the center of the energy bins
  for (Int_t i=0; i<numBins; i++) {
    Double_t En = i*enBinWidth+enBinWidth/2.;
    enBinMedian.push_back(En);
  }

  //Loading theory systematics... These aren't dependent on grouping, only energy of bin
  std::vector <Double_t> theoryCorr = LoadTheoryCorrections(enBinMedian);
  
  std::vector < std::vector <Double_t> > groupRawAsymByBin(3,std::vector <Double_t> (numBins,0.));
  std::vector < std::vector <Double_t> > groupAsymByBin(3,std::vector <Double_t> (numBins,0.));

  
  //bool isOctet = groupType=="Octet" ? true : false;
  bool isQuartet = groupType=="Quartet" ? true : false;
  bool isPair = groupType=="Pair" ? true : false;
  std::map < std::string, Int_t > numQuartetLoops = {{"Octet",1},{"Quartet",2},{"Pair",2}};
  std::map < std::string, Int_t > numPairLoops = {{"Octet",1},{"Quartet",1},{"Pair",2}};
  std::string quartetName[2] = {"A","B"};
  
  std::string basePath = simulation ? getenv("SIM_ANALYSIS_RESULTS") : UKdata ? getenv("ANALYSIS_RESULTS") : getenv("MPM_ANALYSIS_RESULTS");
  
  std::string txt;
  Double_t binEdge, Asym, AsymError;

  Int_t pair = 0, quartet = 0;

  for (Int_t octet=octBegin; octet<=octEnd; octet++) {

    if (std::find(badOct.begin(), badOct.end(),octet) != badOct.end()) { pair+=4; quartet+=2; continue; } //Checking if octet should be ignored for data quality reasons
    
    for (Int_t q = 0; q<numQuartetLoops[groupType]; q++) {

      for (Int_t p=0; p<numPairLoops[groupType]; p++) {

	std::string path = ( basePath + "Octet_" + itos(octet) + "/" + groupType + "Asymmetry/" + (UNBLIND?"UNBLINDED_":"") + "rawAsymmetry_Octet" + 
			     itos(octet) + "_AnaCh" + anaChoice  + (isQuartet?("_Quartet_"+quartetName[q]):isPair?("_Pair_"+quartetName[q]+itos(p)):"")+
			     ( key/10==1 ? "_Quadrant"+itos(key-10) : ( key/10==2 ? "_RadialRing"+itos(key-20) :"" ) ) + ".dat" );
      

	std::ifstream infile(path.c_str());  
   
	if (infile.is_open()) {	  
	  
	  std::vector < std::vector <Double_t> > deltaSys(enBinMedian.size(),std::vector<Double_t>(2,1.));
	  //if (groupType=="Octet") deltaSys = LoadOctetSystematics(octet,anaChoice,enBinMedian);	  
	  std::vector <Double_t> fieldCorr = LoadFieldDipCorrections(enBinMedian,octet);

	  Int_t i = 0;
	  while (infile >> binEdge >> Asym >> AsymError) {
	    groupRawAsymByBin[0][i] = binEdge;
	    groupRawAsymByBin[1][i] += AsymError>0. ? 1./power(AsymError*fieldCorr[i]*deltaSys[i][0]/theoryCorr[i],2)*Asym*fieldCorr[i]*deltaSys[i][0]/theoryCorr[i] : 0.; //Applying Delta_exp here.. Should this affect AsymError???
	    groupRawAsymByBin[2][i] += AsymError>0. ? 1/power(AsymError*fieldCorr[i]*deltaSys[i][0]/theoryCorr[i],2) : 0.;         // Since Asym is multiplied by delta, AsymError would be as well...
	    //std::cout << binEdge << " " << groupRawAsymByBin[1][i] << " " << groupRawAsymByBin[2][i] << std::endl;
	    i++;
	  }
	  infile.close();
	}
	pair++;
      }
      quartet++;
    }
  }
  //Do final calculations of the rates in each bin and their associated errors
  for (unsigned int i=0; i<groupRawAsymByBin[1].size(); i++) {
    groupRawAsymByBin[1][i] = groupRawAsymByBin[2][i]>0. ? groupRawAsymByBin[1][i] / groupRawAsymByBin[2][i] : 0.; // This is sum of weights*Asym / sum of weights 
    groupRawAsymByBin[2][i] = groupRawAsymByBin[2][i]>0. ? 1./sqrt(groupRawAsymByBin[2][i]) : 0.; // Sqrt of sum of weights
    
    groupAsymByBin[0][i] = groupRawAsymByBin[0][i];
    groupAsymByBin[1][i] = 2*groupRawAsymByBin[1][i]/returnBeta(enBinMedian[i]); //Divide out the energy dependence...
    groupAsymByBin[2][i] = 2*groupRawAsymByBin[2][i]/returnBeta(enBinMedian[i]);

    // Apply polarimetry correction if necessary
    if (withPOL) {
      groupRawAsymByBin[1][i] = groupRawAsymByBin[1][i] / ( POL_ave ); 
      groupRawAsymByBin[2][i] = groupRawAsymByBin[2][i] / ( POL_ave );
      groupAsymByBin[1][i] = groupAsymByBin[1][i] / ( POL_ave ); 
      groupAsymByBin[2][i] = groupAsymByBin[2][i] / ( POL_ave );
    }
    
    std::cout << enBinMedian[i] << " " << groupAsymByBin[1][i] << " " << groupAsymByBin[2][i] << std::endl;
  }
  
  // Plotting stuff
  std::string outFile = std::string("asymm_hold.txt");
  
  std::ofstream asymFile(outFile.c_str());
	
  TGraphErrors *gBeta2 = new TGraphErrors(enBinMedian.size(), &enBinMedian[0], &groupAsymByBin[1][0], 0, &groupAsymByBin[2][0]);


  
  TF1 *fitBeta2 = new TF1("fitBeta2","[0]",Elow, Ehigh);
  fitBeta2->SetParameter(0,-0.12);
  
  gBeta2->Fit("fitBeta2","R");

  asymFile << "Asymm\t" << fitBeta2->GetParameter(0) << "\t" << fitBeta2->GetParError(0);
 
  asymFile.close();

  delete gBeta2; delete fitBeta2;

};



  
std::vector < std::vector <Double_t> > LoadOctetSystematics(Int_t octet, std::string anaChoice, std::vector <Double_t> enBinMidpoint) {

  Int_t iAnaChoice;

  for (UInt_t i = 0; i<10; i++) {
    
    if (anaChoices[i]==anaChoice) iAnaChoice = i+1;

  }
  
  //TString filename = TString::Format("%s/Octet_%i/OctetAsymmetry/Systematics/ThOverProc_Octet-%i_Analysis-%i.txt",getenv("ANALYSIS_RESULTS"),octet,octet,iAnaChoice);
  TString filename = TString::Format("%s/systematics/MC_Corrections/DeltaExp_OctetByOctetCorrections/ThOverProc_Octet-%i_Analysis-%i.txt",getenv("ANALYSIS_CODE"),octet,iAnaChoice);
  //std::cout << filename.Data() << std::endl;
  std::vector < std::vector <Double_t> > syst(enBinMidpoint.size(), std::vector<Double_t>(2,1.));

  if ( corr!=std::string("DeltaExpOnly") && corr!=std::string("AllCorr") ) return syst;
  
  std::ifstream infile(filename.Data());

  if (!infile.is_open()) throw "Couldn't open file in LoadOctetSystematics!";

  //Read in the header crap
  //  std::string firstline[8];
  //for (int i=0; i<8; i++) {
  //  infile >> firstline[i];
  //  std::cout << firstline[i] << " " ;
  // }
  std::cout << "\n";

  //Read in the systematics
  Double_t mid; //midpoint of bin
  Double_t ratio;  //systematics correction (Apure/Aproc)
  Double_t midErr; //BinWidth/2
  Double_t ratioErr; //Systematic error on ratio

  Int_t it = 0;

  while (infile >> mid >> ratio >> ratioErr) {
    if (mid==enBinMidpoint[it]) {
      syst[it][0] = ratio!=0. ? ratio : 1.;
      syst[it][1] = ratioErr;
      std::cout << mid << " " << ratio << " " << ratioErr << "\n";
      it++;
    }
  }

  std::cout << "Loaded Systematic Errors for Octet " << octet << " ...\n";

  return syst;
};

std::vector <Double_t> LoadTheoryCorrections(std::vector <Double_t> enBinMidpoint) {

  std::vector <Double_t> syst(enBinMidpoint.size(), 1.);

  if ( corr!=std::string("DeltaTheoryOnly") && corr!=std::string("AllCorr") ) return syst;

  for (UInt_t i=0; i<syst.size(); i++) {
   
    syst[i] = asymmetryCorrectionFactor(enBinMidpoint[i]); //As defined in BetaSpectrum.hh by MPM

  }

  return syst;

};

std::vector <Double_t> LoadFieldDipCorrections(std::vector <Double_t> enBinMidpoint,Int_t oct) {
  std::vector <Double_t> syst(enBinMidpoint.size(), 1.);
  if ( corr!=std::string("DeltaFieldDip") && corr!=std::string("AllCorr") ) return syst;
  
  TString filename = TString::Format("FieldDipCorr_%s.txt",oct<60?"2011-2012":"2012-2013");
  //std::cout << filename.Data() << std::endl;
  
  std::ifstream infile(filename.Data());

  if (!infile.is_open()) throw "Couldn't open file in LoadFieldDipCorrections!";

  std::cout << "\n";

  //Read in the systematics
  Double_t mid; //midpoint of bin
  Double_t ratio;  //systematics correction (Apure/Aproc)
  
  Int_t it = 0;

  while (infile >> mid >> ratio) {
    if (mid==enBinMidpoint[it]) {
      syst[it] = ratio!=0. ? ratio : 1.;
      std::cout << mid << " " << ratio << "\n";
      it++;
    }
  }

  std::cout << "Loaded Field Dip Errors for Octet " << oct << " ...\n";

  return syst;


};
