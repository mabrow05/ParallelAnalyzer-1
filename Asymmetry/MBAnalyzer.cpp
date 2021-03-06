/*
Compilable code which houses the heart of the analysis. 
Here are classes which handle calculating asymmetries of either data 
or simulation. There will be flags to be passed dictating what type of 
asymmetry, what type of data, and what octet/runs to use. 

The code should be able to do BG subtraction, construct asymmetries, 
or do both depending on the flags passed.

*/

#include "EvtRateHandler.hh"
#include "Asymmetries.hh"
#include "SystematicCorrections.hh"
#include "MBUtils.hh"
#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <cmath>
#include <fstream>
#include <algorithm>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <TStyle.h>
#include <TMath.h>

#include <iomanip>

#include "BetaSpectrum.hh"

TString anaChoices[10] = {"A","B","C","D","E","F","G","H","J","K"};


//Types of Corrections to apply
std::string corr ("UnCorr");//{"UnCorr","DeltaExpOnly","DeltaTheoryOnly","AllCorr"};
                             

bool withPOL = true; //Set this to true to correct DATA for the polarimetry measurement
bool dualPol = true; //Whether or not to use both pol measurements per Albert's Rx, or use ave.

bool realOutput = true; // Whether or not the real output should be recorded. If false, just a temporary file is made in systematics/

bool unblind = false;

Double_t POL_minus2011 = 0.997; Double_t POL_minus2011err = 0.003*POL_minus2011;
Double_t POL_plus2011 = 0.9939; Double_t POL_plus2011err = 0.0025*POL_plus2011;
Double_t delta_POL2011 = (POL_plus2011-POL_minus2011)/2.;
Double_t POL_ave2011 = (POL_plus2011+POL_minus2011) / 2.;

Double_t POL_minus2012 = 0.9979; Double_t POL_minus2012err = 0.0015*POL_minus2012;
Double_t POL_plus2012 = 0.9952;  Double_t POL_plus2012err = 0.0020*POL_plus2012;
Double_t delta_POL2012 = (POL_plus2012-POL_minus2012)/2.;
Double_t POL_ave2012 = (POL_plus2012+POL_minus2012) / 2.;


std::vector <Int_t> badOct {7,9,59,60,61,62,63,64,65,66,67,91,93,101,107,121}; // This includes all
                                  // Partial octets
                                  // Octet 7 had W anode dead for part of run
                                  // Either need to discard, or apply the 
                                  // Charge cloud method to determine a 
                                  // coincidence
                                  // 60,61,62,63,64,65,66 have bad TDC spectra
                                  // 70 and 92 just have low statistics so they had bad corrections. Put them back in when using high statistics corrections

// The Process functions will calculate all the raw asymmetries on a bin-by-bin basis and write them to file
void ProcessOctets(Int_t octBegin, Int_t octEnd, std::string anaChoice, Double_t enBinWidth=10., bool UKdata=true, bool simulation=false, bool UNBLIND=false);
void ProcessQuartets(Int_t octBegin, Int_t octEnd, std::string anaChoice, Double_t enBinWidth=10., bool UKdata=true, bool simulation=false, bool UNBLIND=false);
void ProcessPairs(Int_t octBegin, Int_t octEnd, std::string anaChoice, Double_t enBinWidth=10., bool UKdata=true, bool simulation=false, bool UNBLIND=false); 

// makes plots of 2*A/Beta for each octet (pair, quartet, and octet as a whole) and fits over the range specified, for whatever grouping provided ("Octet", "Quartet", "Pair")
void PlotAsymmetriesByGrouping(std::string groupType, Int_t octBegin, Int_t octEnd, std::string anaChoice, Double_t Elow=220., Double_t Ehigh=680., Double_t enBinWidth=10., bool UKdata=true, bool simulation=false, bool UNBLIND=false, int key=0);

//Returns the values of the systematic corrections and statistical error on this for now
std::vector < std::vector <Double_t> >  LoadOctetSystematics(Int_t octet, std::string anaChoice, std::vector <Double_t> enBinMidpoint);

//Returns a vector containing all the theory corrections to A0 for a particular bin
std::vector <Double_t> LoadTheoryCorrections(std::vector <Double_t> enBinMidpoint);

std::vector <Double_t> LoadAngleCorrections(std::vector <Double_t> enBinMidpoint,Int_t octet,std::string anaCh);
std::vector <Double_t> LoadBackscCorrections(std::vector <Double_t> enBinMidpoint,Int_t octet,std::string anaCh);


// Collects all the asymmetries, bin-by-bin, and produces a final asymmetry plot, both A_SR and 2*A/Beta, for whatever grouping provided ("Octet", "Quartet", "Pair"). Also makes a plot of the 
// integrated asymmetry vs octet/quartet/pair number
void PlotFinalAsymmetries(std::string groupType, Int_t octBegin, Int_t octEnd, std::string anaChoice, Double_t Elow=220., Double_t Ehigh=680., Double_t enBinWidth=10., bool UKdata=true, bool simulation=false, bool UNBLIND=false, int key=0);

// notes on key value...
// key==0 : Normal asymmetries
// key>=10 : Quadrant asymmtries, where 10 is quadrant0, 11 is quadrant 1, etc
// key>=20 : Radial Asymmetries, where 20 is radial cut 0, etc (in 10 mm incremements)
// key>=30 : Strip Asymmetries, where 30 is center strip, 31 is elsewhere(Strip is +/- 12 mm



// Function to return beta when given the kinetic energy of an electron
Double_t returnBeta(Double_t En) { 
  Double_t me = 510.998928; //rest mass energy of electron in keV
  return sqrt(En*En+2.*En*me)/(En+me);
};

Double_t backoutRfromA(Double_t A) {
  return TMath::Power(( (1.-A)/(1.+A) ),2.);
};

Double_t backoutDeltaRfromA(Double_t A,Double_t Aerr) {
  return TMath::Abs(4.*(A-1.)/TMath::Power((A+1.),3.)*Aerr);
};

Double_t dA_dP1(Double_t gamma, Double_t P1, Double_t P2) {
  return ( (2.*P1*P2-gamma*gamma*(P1*P2+P2*P2))/(2.*P2*P1*P1*TMath::Sqrt(gamma*gamma*(P1+P2)*(P1+P2)-4.*P1*P2)) 
	   + gamma/(2.*P1*P1) );
};
Double_t polErrorInEnergyBin(Double_t R, Double_t p1, Double_t p1err, Double_t p2, Double_t p2err) {
  Double_t gamma = (R+1.)/(R-1.);
  return TMath::Sqrt( TMath::Power(dA_dP1(gamma,p1,p2)*p1err,2.) + TMath::Power(dA_dP1(gamma,p2,p1)*p2err,2.) ) ;
};




int main(int argc, char* argv[])
{
  
  //if ( corr.size()>7 && corr.compare(corr.size()-8,8,"_withPOL") == 0 ) withPOL = true; 

  /*  std::cout << "Is this thing broken???\n";
  std::vector<double> vec{1.,3.,2.,7.,8.,4.,5.,6.,5.,10.};

    for ( auto i : vec ) {
      std::cout << i << " ";
    }
    std::cout << std::endl;
    
    std::vector<double> vecSorted = sortVecDouble(vec,true);

    for ( auto i : vecSorted ) {
      std::cout << i << " ";
    }

    std::cout << std::endl;

    std::vector <int> rn {17126};
  SimEvtRateHandler evt(rn,true,"C",10.,50.,false);
  evt.CalcRates();
  exit(0);
   
  

  */	
  

  if (argc==1) {
    std::cout << "USAGE: ./MBAnalyzer.exe [anaChoice] [octmin] [octmax] [analysis window min] [analysis window max] [corrections] [bool simulation] [key] [realOutput=true] [withPol=true] [UNBLIND=false] \n";
    exit(0);
  }
    
  std::string analysisChoice = argc>1 ? std::string(argv[1]) : "A";
  Int_t octBegin = argc>2 ? atoi(argv[2]) : 0;
  Int_t octEnd = argc>2 ? atoi(argv[3]) : 1;
  Double_t enBinWidth = 10.;
  Double_t Elow = argc>4 ? atoi(argv[4]) : 220.;
  Double_t Ehigh = argc>4 ? atoi(argv[5]) : 670.;
  if ( argc>6 ) corr = std::string(argv[6]);
  int key = 0;
  if ( argc>8 ) key = atoi(argv[8]);
  bool UKdata = true;//true;
  bool simulation = false;
  if ( argc>7 ) simulation = std::string(argv[7])==std::string("true") ? true: false;
  bool applyAsymm = false;

 
  if (argc>9) realOutput = TString(argv[9])==TString("true")?true:false;
  if (argc>10) withPOL = TString(argv[10])==TString("true")?true:false;
  if (argc>11) unblind = TString(argv[11])==TString("true")?true:false;

  if (simulation) withPOL=false;

  //I should keep track of the raw asymmetry, Experimental systematic corrected asymmetry, and theoretical 
  // Systematic asymmetry in a file each time I run to see the percent correction of each...

  //NEED TO ADD IN ABILITY TO READ IN SYSTEMATICS FOR QUARTET AND PAIR


  if (unblind && realOutput) {
    std::string decision;
    std::cout << "YOU ARE ABOUT TO DO SOME SORT OF UNBLINDING!!!! \n\n";
    std::cout << "To continue, type YES: ";
    std::cin >> decision;

    if (decision!=std::string("YES")) exit(0);

  }

  
  try {
    
    /*std::vector < Double_t > enBinMedian; //Holds the center of the energy bins
    for (Int_t i=0; i<120; i++) {
      Double_t En = i*enBinWidth+enBinWidth/2.;
      enBinMedian.push_back(En);
    }    
    std::vector <Double_t> theoryCorr = LoadTheoryCorrections(enBinMedian);    
    for (UInt_t i=0; i<theoryCorr.size(); i++) std::cout << enBinMedian[i] << " " << theoryCorr[i] << "\n";*/
    

    if (!realOutput) { //This is for systematic study in analysisWindow.py
      //PlotAsymmetriesByGrouping("Octet",octBegin, octEnd, analysisChoice, Elow, Ehigh, enBinWidth, UKdata, simulation, UNBLIND, key);
      PlotFinalAsymmetries("Octet",octBegin, octEnd, analysisChoice, Elow, Ehigh, enBinWidth, UKdata, simulation, unblind, key);
      exit(0);
    }
    

    std::vector<TString> aCh {};//{"A","B","C","D","F","G","H","J","K"};//{"C","D"};//;//{"J","K","C","D","F"};//{"A","B","C","D","F","G","H","J","K"};//{"G","H"};//{"J","K","G"};//{"F","A","H"};//{"A","B","G","H"};//{"C","J","K","H"};//"A","D"

    for (auto ach : aCh) {
      //ProcessOctets(octBegin, octEnd, std::string(ach.Data()), enBinWidth, UKdata, simulation, unblind);
      //ProcessPairs(octBegin, octEnd, std::string(ach.Data()), enBinWidth, UKdata, simulation, unblind);
      PlotAsymmetriesByGrouping("Octet",octBegin, octEnd, std::string(ach.Data()), Elow, Ehigh, enBinWidth, UKdata, simulation, unblind);
      PlotFinalAsymmetries("Octet",octBegin, octEnd, std::string(ach.Data()), Elow, Ehigh, enBinWidth, UKdata, simulation, unblind);
    }
    
    // Loop over keys
    /*int keys[] {0,10,11,12,13,20,21,22,23,24,25,30,31};

    for ( auto& k : keys ) {
      PlotAsymmetriesByGrouping("Octet",octBegin, octEnd, analysisChoice, Elow, Ehigh, enBinWidth, UKdata, simulation, unblind, k);
      PlotFinalAsymmetries("Octet",octBegin, octEnd, analysisChoice, Elow, Ehigh, enBinWidth, UKdata, simulation, unblind, k);
      }*/

   
    //ProcessOctets(octBegin, octEnd, analysisChoice, enBinWidth, UKdata, simulation, unblind);
    PlotAsymmetriesByGrouping("Octet",octBegin, octEnd, analysisChoice, Elow, Ehigh, enBinWidth, UKdata, simulation, unblind, key);
    PlotFinalAsymmetries("Octet",octBegin, octEnd, analysisChoice, Elow, Ehigh, enBinWidth, UKdata, simulation, unblind, key);
    

    //ProcessQuartets(octBegin, octEnd, analysisChoice, enBinWidth, UKdata, simulation, unblind);
    //PlotAsymmetriesByGrouping("Quartet",octBegin, octEnd, analysisChoice, Elow, Ehigh, enBinWidth, UKdata, simulation, unblind, key);
    //PlotFinalAsymmetries("Quartet",octBegin, octEnd, analysisChoice, Elow, Ehigh, enBinWidth, UKdata, simulation, unblind, key);
     
    //ProcessPairs(octBegin, octEnd, analysisChoice, enBinWidth, UKdata, simulation, unblind);
    //PlotAsymmetriesByGrouping("Pair",octBegin, octEnd, analysisChoice, Elow, Ehigh, enBinWidth, UKdata, simulation, unblind, key);
    //PlotFinalAsymmetries("Pair",octBegin, octEnd, analysisChoice, Elow, Ehigh, enBinWidth, UKdata, simulation, unblind, key);
    
    
  }
  catch(const char* ex){
    std::cerr << "Error: " << ex << std::endl;
  }
  
  
  return 0;
}


void ProcessOctets(Int_t octBegin, Int_t octEnd, std::string anaChoice, Double_t enBinWidth, bool UKdata, bool simulation, bool UNBLIND) {

  //std::ofstream octval("octvalUK.dat");
  //octval << "oct" << "\t" 
  //	 << "Asymm" << "\t" 
  //	 << "Error" << "\t" 
  //	 << "Pull" << "\n";
  

  for ( Int_t octet=octBegin; octet<=octEnd; octet++ ) {
    if ( std::find(badOct.begin(), badOct.end(),octet) != badOct.end() ) continue;  //Checking if octet should be ignored for data quality reasons
    try {
      OctetAsymmetry oct(octet,anaChoice,enBinWidth, 50., UKdata, simulation, UNBLIND);
      //oct.calcTotalAsymmetry(180.,780.);
      //oct.calcAsymmetryBinByBin(); 
      //oct.calcNCSUSumAsymmetryBinByBin(); 
      oct.calcSuperSum();
      //oct.calcSuperSumNCSUstyle();
      //oct.writeAsymToFile();
      oct.writeSuperSumToFile(); 
      
      //  octval << octet << "\t" 
      //     << oct.returnTotalAsymmetry() << "\t" 
      //     << oct.returnTotalAsymmetryError() << "\t" 
      //     << 0. << "\n";
      
    }
    catch(const char* ex){
      std::cerr << "Error: " << ex << std::endl;
    }

  }

  //octval.close();

};

void ProcessQuartets(Int_t octBegin, Int_t octEnd, std::string anaChoice, Double_t enBinWidth, bool UKdata, bool simulation, bool UNBLIND) {
  
  for (Int_t octet=octBegin; octet<=octEnd; octet++) {
    if ( std::find(badOct.begin(), badOct.end(),octet) != badOct.end() ) continue;  //Checking if octet should be ignored for data quality reasons
    try {
      QuartetAsymmetry quart(octet,anaChoice,enBinWidth, 50., UKdata, simulation, UNBLIND);
      quart.calcAsymmetryBinByBin();
      //quart.calcSuperSum();
      quart.writeAsymToFile();
      //quart.writeSuperSumToFile();
    }
    catch(const char* ex){
      std::cerr << "Error: " << ex << std::endl;
    }
  }
  
};

void ProcessPairs(Int_t octBegin, Int_t octEnd, std::string anaChoice, Double_t enBinWidth, bool UKdata, bool simulation, bool UNBLIND) {
  
  for (Int_t octet=octBegin; octet<=octEnd; octet++) {
    if ( std::find(badOct.begin(), badOct.end(),octet) != badOct.end() ) continue;  //Checking if octet should be ignored for data quality reasons
    try {
      PairAsymmetry pair(octet, anaChoice, enBinWidth, 50., UKdata, simulation, UNBLIND);
      pair.calcAsymmetryBinByBin();
      //pair.calcSuperSum();
      pair.writeAsymToFile();
      //pair.writeSuperSumToFile();
    }
    catch(const char* ex){
      std::cerr << "Error: " << ex << std::endl;
    }
  }
  
};


void PlotFinalAsymmetries(std::string groupType, Int_t octBegin, Int_t octEnd, std::string anaChoice, Double_t Elow, Double_t Ehigh, Double_t enBinWidth, bool UKdata, bool simulation, bool UNBLIND, int key) {


  if (groupType!=std::string("Quartet") && groupType!=std::string("Octet") 
      && groupType!=std::string("Pair") )
    throw "Bad group type given to PlotFinalAsymmetries. Options are \"Octet\", \"Quartet\", or \"Pair\""; 

  key = ( groupType==std::string("Octet") ? key : 0 );

  Int_t numBins = 1200/enBinWidth;

  std::vector < Double_t > enBinMedian; //Holds the center of the energy bins
  for (Int_t i=0; i<numBins; i++) {
    Double_t En = i*enBinWidth+enBinWidth/2.;
    enBinMedian.push_back(En);
  }

  //Loading theory systematics... These aren't dependent on grouping, only energy of bin
  std::vector <Double_t> theoryCorr = LoadTheoryCorrections(enBinMedian);

  std::vector < std::vector <Double_t> > rawAsymByGroup(3,std::vector <Double_t> (0));
  std::vector < std::vector <Double_t> > AsymByGroup(3,std::vector <Double_t> (0));
  
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

	std::string path = ( basePath+"Octet_"+itos(octet)+"/"+groupType+"Asymmetry/"+
			     (UNBLIND?"UNBLINDED_":"") + corr + "_" + (withPOL?"withPOL_":"") +
			     "FittedAsymmetry_Octet"+itos(octet)+"_AnaCh"+anaChoice+"_"+
			     (isQuartet?("Quartet_"+quartetName[q]+"_"):isPair?
			      ("Pair_"+quartetName[q]+itos(p)+"_"):"")+itos((int)Elow)+
			     "-"+itos((int)Ehigh) + 
			     ( key/10==1 ? "_Quadrant"+itos(key-10) : 
			       ( key/10==2 ? "_RadialRing"+itos(key-20) :
				 ( key/10==3 ? "_Strip"+itos(key-30) :"" ) ) ) + ".dat");
	//std::cout << path << std::endl;

	std::ifstream infile(path.c_str());  
   
	if (infile.is_open() || !realOutput) {

	  if (realOutput) {
	    infile >> txt >> Asym >> AsymError;
	    rawAsymByGroup[0].push_back(isQuartet?quartet:isPair?pair:octet);
	    rawAsymByGroup[1].push_back(Asym); // Note the negative sign here to turn the raw asymmetry into purely a positive number (by definition, it is negative)
	    rawAsymByGroup[2].push_back(AsymError);
	    
	    //std::cout << octet << " " << Asym << " " << AsymError << std::endl;
	    
	    infile >> txt >> Asym >> AsymError;
	    AsymByGroup[0].push_back(isQuartet?quartet:isPair?pair:octet);
	    AsymByGroup[1].push_back(Asym);
	    AsymByGroup[2].push_back(AsymError);
	    
	  }
	  infile.close();
	  path = ( basePath + "Octet_" + itos(octet) + "/" + groupType + "Asymmetry/" +
		   (UNBLIND?"UNBLINDED_":"") + "rawAsymmetry_Octet" + itos(octet) + "_AnaCh" + anaChoice  +
		   (isQuartet?("_Quartet_"+quartetName[q]):isPair?("_Pair_"+quartetName[q]+itos(p)):"")+
		   ( key/10==1 ? "_Quadrant"+itos(key-10) : 
		     ( key/10==2 ? "_RadialRing"+itos(key-20) :
		       ( key/10==3 ? "_Strip"+itos(key-30) :"" ) ) ) + ".dat" );
	  
	  infile.open(path.c_str());
	  //std::cout << path << std::endl;
	  
	  
	  Int_t i = 0;
	  while (infile >> binEdge >> Asym >> AsymError) {

	    // Apply polarimetry correction if necessary
	    /*if (withPOL) {
	      if (dualPol) {
		double P_p =  POL_plus2011; 
		double P_m =  POL_minus2011;
		if (octet>59) P_p = POL_plus2012; P_m = POL_minus2012;
		
		Double_t R = backoutRfromA(Asym);
		Double_t delta_R = backoutDeltaRfromA(Asym,AsymError);
		Double_t gam = (R+1.)/(R-1.);
		Double_t delta_gam = TMath::Abs(2.*delta_R/TMath::Power(R-1.,2.));
		
	        Asym = (-gam*(P_p+P_m)+TMath::Sqrt(gam*gam*(P_p+P_m)*(P_p+P_m)-4.*P_p*P_m))/(2.*P_p*P_m);
		AsymError = TMath::Abs((-(P_p+P_m)/(2.*P_p*P_m) + gam*(P_p+P_m)*(P_p+P_m)/(2.*P_p*P_m*TMath::Sqrt(gam*gam*(P_p+P_m)*(P_p+P_m)-4.*P_p*P_m)))*delta_gam);	
	      } else {
		double POL_ave=POL_ave2011;
		if (octBegin>59 && octEnd>59)  POL_ave = POL_ave2012;
	        Asym /= ( POL_ave ); 
		AsymError /= ( POL_ave );
	      }
	      }*/

	    groupRawAsymByBin[0][i] = binEdge;
	    groupRawAsymByBin[1][i] += AsymError>0. ? 1./power(AsymError,2)*Asym: 0.; //Applying Delta_exp here.. Should this affect AsymError???
	    groupRawAsymByBin[2][i] += AsymError>0. ? 1/power(AsymError,2) : 0.;         // Since Asym is multiplied by delta, AsymError would be as well...
	    std::cout << binEdge << " " << groupRawAsymByBin[1][i] << " " << groupRawAsymByBin[2][i] << std::endl;

	    i++;
	  }
	  infile.close();
	}
	pair++;
      }
      quartet++;
    }
  }

  //Apply systematics if run is an octet
  //TODO: put in systematics for other types of runs...
  //std::vector < std::vector <Double_t> > deltaSys(enBinMedian.size(),std::vector<Double_t>(2,1.));
  
  //if (groupType==std::string("Octet")) deltaSys = LoadOctetSystematics(octet,anaChoice,enBinMedian);
  std::vector <Double_t> angleCorr = LoadAngleCorrections(enBinMedian,octBegin,anaChoice);
  std::vector <Double_t> bsCorr = LoadBackscCorrections(enBinMedian,octBegin,anaChoice);

  std::vector <Double_t> polErrors(enBinMedian.size(),0.);

  //Do final calculations of the rates in each bin and their associated errors
  for (unsigned int i=0; i<groupRawAsymByBin[1].size(); i++) {
    groupRawAsymByBin[1][i] = groupRawAsymByBin[2][i]>0. ? groupRawAsymByBin[1][i] / groupRawAsymByBin[2][i] : 0.; // This is sum of weights*Asym / sum of weights 
    groupRawAsymByBin[2][i] = groupRawAsymByBin[2][i]>0. ? 1./sqrt(groupRawAsymByBin[2][i]) : 0.; // Sqrt of sum of weights
    
    // Apply polarimetry correction if necessary
    if (withPOL) {
      if (dualPol) {
	double P_p =  POL_plus2011; 
	double P_m =  POL_minus2011;
	double P_p_err = POL_plus2011err;
	double P_m_err = POL_minus2011err;
	if (octBegin>59) {P_p = POL_plus2012; P_m = POL_minus2012; P_p_err = POL_plus2012err; P_m_err = POL_minus2012err;}
	
	Double_t R = backoutRfromA(groupRawAsymByBin[1][i]);
	Double_t delta_R = backoutDeltaRfromA(groupRawAsymByBin[1][i],groupRawAsymByBin[2][i]);
	Double_t gam = (R+1.)/(R-1.);
	Double_t delta_gam = TMath::Abs(2.*delta_R/TMath::Power(R-1.,2.));
	
	groupRawAsymByBin[1][i] = (-gam*(P_p+P_m)+TMath::Sqrt(gam*gam*(P_p+P_m)*(P_p+P_m)-4.*P_p*P_m))/(2.*P_p*P_m);
	groupRawAsymByBin[2][i] = TMath::Abs((-(P_p+P_m)/(2.*P_p*P_m) + gam*(P_p+P_m)*(P_p+P_m)/(2.*P_p*P_m*TMath::Sqrt(gam*gam*(P_p+P_m)*(P_p+P_m)-4.*P_p*P_m)))*delta_gam);
	polErrors[i] = 	TMath::Abs( groupRawAsymByBin[1][i]!=0.?polErrorInEnergyBin(R, P_p, P_p_err, P_m, P_m_err)/groupRawAsymByBin[1][i]:1. );
      } else {
	double POL_ave=POL_ave2011;
	if (octBegin>59)  POL_ave = POL_ave2012;
	groupRawAsymByBin[1][i] /= ( POL_ave ); 
	groupRawAsymByBin[2][i] /= ( POL_ave );
      }
    }
      
    //Apply systematic corrections
    groupRawAsymByBin[1][i] = groupRawAsymByBin[1][i]*angleCorr[i]*bsCorr[i]/theoryCorr[i];
    groupRawAsymByBin[2][i] = groupRawAsymByBin[2][i]*angleCorr[i]*bsCorr[i]/theoryCorr[i];
    
    groupAsymByBin[0][i] = groupRawAsymByBin[0][i];
    groupAsymByBin[1][i] = 2*groupRawAsymByBin[1][i]/returnBeta(enBinMedian[i]); //Divide out the energy dependence...
    groupAsymByBin[2][i] = 2*groupRawAsymByBin[2][i]/returnBeta(enBinMedian[i]);
    
    //std::cout << enBinMedian[i] << " " << groupAsymByBin[1][i] << " " << groupAsymByBin[2][i] << std::endl;
  }



  // Plotting stuff
  std::string outFile;
  if (realOutput) outFile = ( basePath + "Asymmetries/"+(UNBLIND?"UNBLINDED_":"") + corr + "_" +
			      (withPOL?"withPOL_":"") + groupType +"Asymmetries_AnaCh" + anaChoice 
			      + std::string("_") + itos((int)Elow) + std::string("-") +
			      itos((int)Ehigh) + "_Octets_" +itos(octBegin)+"-"+itos(octEnd)+
			      ( key/10==1 ? "_Quadrant"+itos(key-10) : 
				( key/10==2 ? "_RadialRing"+itos(key-20) :
				  ( key/10==3 ? "_Strip"+itos(key-30) :"" ) ) ) );
  
  else outFile = "../UNBLINDING/asymm_hold";

  std::string pdfFile = outFile+std::string(".pdf");
  std::string txtFile = outFile+std::string(".txt");

  std::ofstream asymFile(txtFile.c_str());
  asymFile << std::setprecision(10);
  
  gStyle->SetTitleSize(0.1,"t");
  gStyle->SetTitleSize(0.06,"x");
  gStyle->SetTitleSize(0.06,"y");
  gStyle->SetTitleOffset(1.2,"y");
  gStyle->SetTitleOffset(0.9,"x");
  gStyle->SetPadTopMargin(0.08);
  gStyle->SetPadBottomMargin(0.13);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadRightMargin(0.02);
  gStyle->SetLabelSize(0.04,"xyz");

  gStyle->SetOptFit(1111);//1111);
  gStyle->SetTitleX(0.5);
  gStyle->SetStatX(0.75);
  gStyle->SetStatY(0.9);
  gStyle->SetStatFontSize(0);
  gStyle->SetStatBorderSize(0);
  
  gStyle->SetFitFormat("0.5g");
  
  
   
  std::vector <double> errorX(1000,0.);
  std::string title;
  
  TGraphErrors *g;
  TF1 *fit;

  TCanvas *c1 = new TCanvas("c1", "c1", 1000, 400);
  if (realOutput) {
    
    gStyle->SetStatH(0.3);
    gStyle->SetStatW(0.27);
    gStyle->SetStatX(0.85);
    gPad->SetLeftMargin(0.12); 
    gPad->SetBottomMargin(0.15); 
    gStyle->SetLabelSize(0.06,"xyz");
   //gPad->SetRightMargin(0.02);
    
    g = new TGraphErrors(rawAsymByGroup[0].size(), &rawAsymByGroup[0][0],&rawAsymByGroup[1][0],&errorX[0], &rawAsymByGroup[2][0]);
    title = "Raw Measured Asymmetry " + itos((int)Elow) + std::string("-") +itos((int)Ehigh) + std::string(" keV Window");
    g->SetTitle("");//title.c_str());
    
    g->SetMarkerStyle(20);
    g->SetLineWidth(1);
    g->Draw("AP");
    g->GetXaxis()->SetLimits(rawAsymByGroup[0][0]-2., rawAsymByGroup[0][rawAsymByGroup[0].size()-1]+2.);
    g->GetXaxis()->SetTitle(TString::Format("%s Number",groupType.c_str()));
    g->GetYaxis()->SetTitle("Raw Asymmetry A_{SR}");
    g->GetYaxis()->SetTitleOffset(0.7);
    g->GetXaxis()->SetTitleOffset(0.8);
    g->GetYaxis()->SetTitleSize(0.08);
    g->GetXaxis()->SetTitleSize(0.08);
    g->GetXaxis()->CenterTitle();
    g->GetYaxis()->CenterTitle();
    
    fit = new TF1("fit","[0]",rawAsymByGroup[0][0]-1., rawAsymByGroup[0][rawAsymByGroup[0].size()-1]+1.);
    fit->SetLineColor(kRed);
    fit->SetLineWidth(2);
    fit->SetParameter(0,0.05);
    fit->SetParName(0,"A");
    
    g->Fit("fit","R");
    
    asymFile << "RawA_oct_by_oct\t" << fit->GetParameter(0) << "\t" << fit->GetParError(0) << std::endl;
  
    //Writing to file the raw asymmetries of each octet
    /*std::ofstream octval("octvalUK.dat");
      octval << "oct" << "\t" 
      << "Asymm" << "\t" 
      << "Error" << "\t" 
      << "Pull" << "\n";
      
      for (UInt_t i=0; i<rawAsymByGroup[0].size(); i++) {
      octval << rawAsymByGroup[0][i] << "\t" 
      << rawAsymByGroup[1][i] << "\t" 
      << rawAsymByGroup[2][i] << "\t" 
      << ( rawAsymByGroup[1][i] - fit->GetParameter(0) ) / rawAsymByGroup[2][i] << "\n";
      } 
      octval.close();*/
    
    g->Draw("AP");
    g->SetMinimum(-0.07);//fit->GetParameter(0)-0.02);//((simulation && !AsymmOn) ? -0.05 : 0.03);
    g->SetMaximum(0.0);//fit->GetParameter(0)+0.02);//(simulation && !AsymmOn) ? 0.05 : 0.07);
    c1->Update();
  }

  gStyle->SetStatH(0.2);
  gStyle->SetStatW(0.3);
  gStyle->SetStatX(0.85);
  gStyle->SetLabelSize(0.04,"xyz");
    
  TGraphErrors *gBeta;
  TF1 *fitBeta;
  std::string corrText = corr=="UnCorr" ? " 2A_{SR}/#beta " : ( corr=="DeltaExpOnly" ? " 2A_{SR}#upoint(1+#Delta_{exp})/#beta " : 
								( corr=="DeltaTheoryOnly" ? " 2A_{SR}#upoint(1+#Delta_{Th})/#beta " : 
								  ( " 2A_{SR}#upoint(1+#Delta_{exp})#upoint(1+#Delta_{Th})/#beta ")));
  
  TCanvas *c2 = new TCanvas("c2", "c2", 700, 500);

  if (realOutput) {
    
    TGraphErrors *gBeta = new TGraphErrors(AsymByGroup[0].size(), &AsymByGroup[0][0],&AsymByGroup[1][0],&errorX[0], &AsymByGroup[2][0]);
    title = groupType+"-By-"+groupType+corrText + itos((int)Elow) + std::string("-") +itos((int)Ehigh) + std::string(" keV Window");
    gBeta->SetTitle("");//(title.c_str());
    gBeta->SetMarkerStyle(20);
    gBeta->SetLineWidth(1);
    gBeta->GetXaxis()->SetLimits(rawAsymByGroup[0][0]-2., rawAsymByGroup[0][rawAsymByGroup[0].size()-1]+2.);
    gBeta->GetXaxis()->SetTitle("Number");
    gBeta->GetYaxis()->SetTitle("Asymmetry");
    gBeta->GetXaxis()->CenterTitle();
    gBeta->GetYaxis()->CenterTitle();
    
    TF1 *fitBeta = new TF1("fitBeta","[0]",rawAsymByGroup[0][0]-1., rawAsymByGroup[0][rawAsymByGroup[0].size()-1]+1.);
    fitBeta->SetLineColor(kRed);
    fitBeta->SetLineWidth(2);
    fitBeta->SetParameter(0,0.05);
    fitBeta->SetParName(0,"A");
	
    gBeta->Fit("fitBeta","R");
    
    asymFile << "BetaCorrA_oct_by_oct\t" << fitBeta->GetParameter(0) << "\t" << fitBeta->GetParError(0) << std::endl;
    
    gBeta->Draw("AP");
    gBeta->SetMinimum(fitBeta->GetParameter(0)-0.02);//(simulation && !AsymmOn) ? -0.5 : -0.14);
    gBeta->SetMaximum(fitBeta->GetParameter(0)+0.03);//(simulation && !AsymmOn) ? 0.5 : -0.09);
    c2->Update();
  }
  
  //c1->Print(pdfFile.c_str());

  std::string corrText2 = ( corr=="UnCorr" ? " A_{SR} " : 
			    ( corr=="DeltaExpOnly" ? " A_{SR}#upoint(1+#Delta_{exp}) " : 
			      ( corr=="DeltaTheoryOnly" ? " A_{SR}#upoint(1+#Delta_{Th}) " : 
				( " A_{SR}#upoint(1+#Delta_{exp})#upoint(1+#Delta_{Th}) "))));

  TCanvas *c3 = new TCanvas("c3", "c3", 700, 500);

  Int_t offset = 4;
  TGraphErrors *g2 = new TGraphErrors(enBinMedian.size()-offset, &enBinMedian[offset], &groupRawAsymByBin[1][offset], 0, &groupRawAsymByBin[2][offset]);
  title = groupType + corrText2 + itos((int)Elow) + std::string("-") +itos((int)Ehigh) + std::string(" keV Window");
  g2->SetTitle("");//(title.c_str());
  g2->SetMarkerStyle(20);
  g2->SetLineWidth(1);
  g2->GetXaxis()->SetLimits(0., 800.);
  g2->GetXaxis()->SetTitle("Energy (keV)");
  g2->GetYaxis()->SetTitle("Asymmetry");
  g2->GetXaxis()->CenterTitle();
  g2->GetYaxis()->CenterTitle();
  
  g2->Draw("AP");
  g2->SetMinimum(-0.1);//(simulation && !AsymmOn) ? -0.05 : -0.08);
  g2->SetMaximum(0.);//(simulation && !AsymmOn) ? 0.05 : -0.01);

  TF1 *f2 = new TF1("f2","(-0.1184/2.*sqrt(x*x+2.*x*510.99892)/(x+510.99892))",0., 1000.);
  f2->SetLineColor(kBlack);
  f2->SetLineStyle(2);

  f2->Draw("SAME");

  c3->Update();
	
  TCanvas *c4 = new TCanvas("c4", "c4", 700, 500);
  TGraphErrors *gBeta2 = new TGraphErrors(enBinMedian.size()-offset, &enBinMedian[offset], &groupAsymByBin[1][offset], 0, &groupAsymByBin[2][offset]);
  title = groupType + corrText + itos((int)Elow) + std::string("-") +itos((int)Ehigh) + std::string(" keV Window");
  gBeta2->SetTitle("");//(title.c_str());
  gBeta2->SetMarkerStyle(20);
  gBeta2->SetLineWidth(1);
  gBeta2->GetXaxis()->SetLimits(0., 800.);
  gBeta2->GetXaxis()->SetTitle("Energy (keV)");
  gBeta2->GetYaxis()->SetTitle("Asymmetry");
  gBeta2->GetXaxis()->CenterTitle();
  gBeta2->GetYaxis()->CenterTitle();
  
  TF1 *fitBeta2 = new TF1("fitBeta2","[0]",Elow, Ehigh);
  fitBeta2->SetLineColor(kRed);
  fitBeta2->SetLineWidth(2);
  fitBeta2->SetParameter(0,-0.12);
  fitBeta2->SetParName(0,"A_{0}");
  
  gBeta2->Fit("fitBeta2","R");

  if (realOutput) asymFile << "BetaCorrA_binSummed\t" << fitBeta2->GetParameter(0) << "\t" << fitBeta2->GetParError(0) << std::endl;
  else asymFile << "BetaCorrA_binSummed\t" << fitBeta2->GetParameter(0) << "\t" << fitBeta2->GetParError(0) << std::endl;
  

  
  gBeta2->Draw("AP");
  gBeta2->SetMinimum(fitBeta2->GetParameter(0)-0.04);//(simulation && !AsymmOn) ? -0.5 : -0.16);
  gBeta2->SetMaximum(fitBeta2->GetParameter(0)+0.05);//(simulation && !AsymmOn) ? 0.5 : -0.07);
  c4->Update();
  if (realOutput) { 
    c1->Print(TString::Format("%s(",pdfFile.c_str()));
    c2->Print(pdfFile.c_str());
    c3->Print(pdfFile.c_str());
    c4->Print(TString::Format("%s)",pdfFile.c_str()));
  }
  asymFile.close();

  //Write out corrected (but energy dependent) bin-by-bin asymmetry for use in calculating effect from doing corrections

  
  outFile = ( basePath + "Asymmetries/"+(UNBLIND?"UNBLINDED_":"") + corr + "_" + (withPOL?"withPOL_":"") + groupType +"Asymmetries_AnaCh" + anaChoice 
	      + std::string("_")  + "Octets_" +itos(octBegin)+"-"+itos(octEnd)+
	      ( key/10==1 ? "_Quadrant"+itos(key-10) : 
		( key/10==2 ? "_RadialRing"+itos(key-20) :
		  ( key/10==3 ? "_Strip"+itos(key-30) :"" ) ) ) );
  
  txtFile = outFile+std::string("_BinByBin_withEnergyDependence.txt");
  asymFile.open(txtFile.c_str());
  asymFile << std::setprecision(15);
  
  for (UInt_t n=0; n<enBinMedian.size(); n++) {
    asymFile << enBinMedian[n] << "\t" << groupRawAsymByBin[1][n] << "\t" << groupRawAsymByBin[2][n] << "\n";
  }
  asymFile.close();
  
  //Write out corrected bin-by-bin asymmetry for use in calculating effect from doing corrections
  txtFile = outFile+std::string("_BinByBin.txt");
  asymFile.open(txtFile.c_str());
  asymFile << std::setprecision(10);
  
  for (UInt_t n=0; n<enBinMedian.size(); n++) {
    asymFile << enBinMedian[n] << "\t" << groupAsymByBin[1][n] << "\t" << groupAsymByBin[2][n] << "\n";
    //std::cout << enBinMedian[n] << "\t" << groupAsymByBin[1][n] << "\t" << groupAsymByBin[2][n] << "\n";
  }
  asymFile.close();


  //Write out polarimetry errors if applicable
  if (withPOL && dualPol) {
    outFile = ( basePath + "Corrections/PolUncert"+ groupType +"Asymmetries_AnaCh" + anaChoice 
		+ std::string("_")  + "Octets_" +itos(octBegin)+"-"+itos(octEnd)+
		( key/10==1 ? "_Quadrant"+itos(key-10) : 
		  ( key/10==2 ? "_RadialRing"+itos(key-20) :
		    ( key/10==3 ? "_Strip"+itos(key-30) :"" ) ) ) +".txt");
    std::ofstream polFile(outFile);
    polFile << std::setprecision(15);
    for (unsigned i=0;i<enBinMedian.size();++i) {
      polFile << enBinMedian[i] << "\t" << polErrors[i] << std::endl;
    }
    polFile.close();
  }
  
  if (true) {
    gStyle->SetTitleSize(0.04,"x");
    gStyle->SetTitleSize(0.04,"y");
    gStyle->SetStatX(0.75);
    gStyle->SetStatY(0.9);
    gStyle->SetStatBorderSize(0);
    //gStyle->SetOptStat(0);
    //gStyle->SetTitleOffset(0.85,"y");
    gStyle->SetPadTopMargin(0.05);
    gStyle->SetPadBottomMargin(0.15);
    gStyle->SetPadLeftMargin(0.17);
    gStyle->SetLabelSize(0.03,"xyz");
    //gStyle->SetPadLeftMargin(0.12);
    //gStyle->SetPadBottomMargin(0.12);
    TCanvas *cBrad = new TCanvas("cBrad","cBrad",700,500);
    gBeta2->SetTitle("");
    gBeta2->SetMarkerStyle(20);
    gBeta2->SetLineWidth(1);
    gBeta2->GetXaxis()->SetLimits(0., 800.);
    gBeta2->GetXaxis()->SetTitle("Energy (keV)");
    gBeta2->GetXaxis()->SetTitleOffset(1.1);
    gBeta2->GetYaxis()->SetTitleOffset(1.3);
    gBeta2->GetYaxis()->SetTitle("Asymmetry");
    gBeta2->GetXaxis()->CenterTitle();
    gBeta2->GetYaxis()->CenterTitle();
  
    //TF1 *fitBeta2 = new TF1("fitBeta2","[0]",Elow, Ehigh);
    //fitBeta2->SetLineColor(kRed);
    //fitBeta2->SetLineWidth(3);
    //fitBeta2->SetParameter(0,-0.12);
    //gBeta2->Fit("fitBeta2","R");
  
    gBeta2->Draw("AP");
    gBeta2->SetMinimum(fitBeta2->GetParameter(0)-0.04);//(simulation && !AsymmOn) ? -0.5 : -0.16);
    gBeta2->SetMaximum(fitBeta2->GetParameter(0)+0.05);//(simulation && !AsymmOn) ? 0.5 : -0.07);
    cBrad->Print("asymmetryForBrad.pdf");
    delete cBrad;
  }

  delete c1, delete c2, delete c3, delete c4; 
  if (realOutput ) { delete g; delete fit; delete gBeta; delete fitBeta; }
  delete g2; delete gBeta2; delete fitBeta2;

};


//This will create asymmetry plots for each pair, quartet, and octet, and also write the raw integrated asymmetry and Beta corrected integrated asymmetry
// to file for each grouping to be used by PlotFinalAsymmetries

void PlotAsymmetriesByGrouping(std::string groupType, Int_t octBegin, Int_t octEnd, std::string anaChoice, Double_t Elow, Double_t Ehigh, Double_t enBinWidth, bool UKdata, bool simulation, bool UNBLIND,int key) {
  if (groupType!="Quartet" && groupType!="Octet" && groupType!="Pair") throw "Bad group type given to PlotAsymmetriesByGrouping. Options are \"Octet\", \"Quartet\", or \"Pair\""; 
  
  key = ( groupType==std::string("Octet") ? key : 0 );

  std::string basePath = simulation ? getenv("SIM_ANALYSIS_RESULTS") : UKdata ? getenv("ANALYSIS_RESULTS") : getenv("MPM_ANALYSIS_RESULTS");

  for (Int_t octet=octBegin; octet<=octEnd; octet++) {

    if (std::find(badOct.begin(), badOct.end(),octet) != badOct.end()) {  continue; } //Checking if octet should be ignored for data quality reasons

    std::ifstream infile;

    if (groupType=="Octet")
    {
      std::string basePath2 =basePath+ "Octet_"+itos(octet)+"/" + groupType + "Asymmetry/"; 
      std::string infilePath = ( basePath2 + (UNBLIND?"UNBLINDED_":"") + "rawAsymmetry_Octet" + 
				 itos(octet) + "_AnaCh" + anaChoice +
				 ( key/10==1 ? "_Quadrant"+itos(key-10) : 
				   ( key/10==2 ? "_RadialRing"+itos(key-20) :
				     ( key/10==3 ? "_Strip"+itos(key-30) :"" ) ) ) + ".dat");
      std::string outfilePath =  ( basePath2 + (UNBLIND?"UNBLINDED_":"")+ corr +"_" + (withPOL?"withPOL_":"") +"FittedAsymmetry_Octet" +
				   itos(octet) + "_AnaCh" + anaChoice + "_" + itos((int)Elow) + std::string("-") +itos((int)Ehigh) +
				   ( key/10==1 ? "_Quadrant"+itos(key-10) : 
				     ( key/10==2 ? "_RadialRing"+itos(key-20) :
				       ( key/10==3 ? "_Strip"+itos(key-30) :"" ) ) ) + ".dat");
      
      //Remove old output files in case they were created on accident and aren't filled with good values
      std::string command = "rm " + outfilePath; 
      system(command.c_str());

      //First check that Octet was good
      std::string checkStatus;
      infile.open(infilePath.c_str());

      if ( infile.is_open() ) {
	infile >> checkStatus; 
	infile.close();
	if (checkStatus=="BAD") continue;
      }
      else {
	std::cout << "Could not open binned Asymmetries for Octet " << octet << " so skipping...\n";
	continue; 
      }

      std::vector < std::vector <Double_t > > AsymAndError;
      std::vector < std::vector <Double_t > > RawAsymAndError;
      std::vector < Double_t > enBinMedian;
      
      AsymAndError.resize(2,std::vector <Double_t> (0));
      RawAsymAndError.resize(2,std::vector <Double_t> (0));
	
      Double_t eBinLow, Asym, AsymError;

      infile.open(infilePath.c_str());

      Int_t bin = 0;

      while (infile >> eBinLow >> Asym >> AsymError) {
	Double_t En = eBinLow+enBinWidth/2.;
	enBinMedian.push_back(En);
	RawAsymAndError[0].push_back(Asym);
	RawAsymAndError[1].push_back(AsymError);
	AsymAndError[0].push_back(0.); // These are just being initialized 
	AsymAndError[1].push_back(0.);
	bin++;
      }
      
      infile.close();
      
      //Read in and apply the systematic corrections for this octet
      //std::vector < std::vector <Double_t> > deltaSys = LoadOctetSystematics(octet,anaChoice,enBinMedian);
      
      //Loading theory systematics... 
      std::vector <Double_t> theoryCorr = LoadTheoryCorrections(enBinMedian);
      std::vector <Double_t> angleCorr = LoadAngleCorrections(enBinMedian,octet,anaChoice);
      std::vector <Double_t> bsCorr = LoadBackscCorrections(enBinMedian,octet,anaChoice);
      
      
      for (UInt_t i=0 ; i<AsymAndError[0].size() ; i++) {
	
	if (withPOL) {
	  if (dualPol) {
	    Double_t P_p =  POL_plus2011; 
	    Double_t P_m =  POL_minus2011;
	    if (octet>59) { P_p = POL_plus2012; P_m =  POL_minus2012;}
	    Double_t R = backoutRfromA(RawAsymAndError[0][i]);
	    Double_t delta_R = backoutDeltaRfromA(RawAsymAndError[0][i],RawAsymAndError[1][i]);
	    Double_t gam = (R+1.)/(R-1.);
	    Double_t delta_gam = TMath::Abs(2.*delta_R/TMath::Power((R-1.),2.));
	    
	    //std::cout << enBinMedian[i] << "\t" << R << "\t" << delta_R << "\t"<< gam << "\t"<< delta_gam << "\n";  
	    
	    RawAsymAndError[0][i] = (-gam*(P_p+P_m)+TMath::Sqrt(gam*gam*(P_p+P_m)*(P_p+P_m)-4.*P_p*P_m))/(2.*P_p*P_m);
	    RawAsymAndError[1][i] = TMath::Abs( (-(P_p+P_m)/(2.*P_p*P_m) + gam*(P_p+P_m)*(P_p+P_m)/(2.*P_p*P_m*TMath::Sqrt(gam*gam*(P_p+P_m)*(P_p+P_m)-4.*P_p*P_m)) )*delta_gam);	  
	  } else {
	    double POL_ave=POL_ave2011;
	    if (octet>59)  POL_ave = POL_ave2012;
	    RawAsymAndError[0][i] /= ( POL_ave ); 
	    RawAsymAndError[1][i] /= ( POL_ave );
	  }
	}
	
	Double_t Beta = returnBeta(enBinMedian[i]);
	AsymAndError[0][i] = RawAsymAndError[0][i]*2./Beta; 
	AsymAndError[1][i] = RawAsymAndError[1][i]*2./Beta; 
	
	AsymAndError[0][i] *= (bsCorr[i]*angleCorr[i] / theoryCorr[i]); //Here is where the corrections to Ameas are made.. Need to add in delta theory
	AsymAndError[1][i] *= (bsCorr[i]*angleCorr[i] / theoryCorr[i]);
      }
  

        //Writing out corrected bin-by-bin asymmetries for Xuan.
      if (realOutput && key==0 && !withPOL) {
	std::string xuanOutFile = ( basePath2 + (UNBLIND?"UNBLINDED_":"") + corr + "_" +
				    groupType +"AsymmetriesByBin_Octet"+itos(octet)+"_AnaCh" + anaChoice +".dat" );
	std::ofstream xuanFile(xuanOutFile.c_str());
	xuanFile << std::setprecision(15);
	for (unsigned i=0;i<AsymAndError[0].size();++i) xuanFile << enBinMedian[i] << "\t" << AsymAndError[0][i] << "\t" << AsymAndError[1][i] << "\n";
	xuanFile.close();
      }

      
      std::string pdfPath = ( basePath2 + (UNBLIND?"UNBLINDED_":"") + corr + "_" +  "Asymmetry_Octet" + itos(octet) +
			      "_AnaCh" + anaChoice + "_" + itos((int)Elow) + std::string("-") +itos((int)Ehigh) +
			      ( key/10==1 ? "_Quadrant"+itos(key-10) : 
				( key/10==2 ? "_RadialRing"+itos(key-20) :
				  ( key/10==3 ? "_Strip"+itos(key-30) :"" ) ) ) + ".pdf");
      TCanvas *c1 = new TCanvas("c1", "c1",800, 400);
      gStyle->SetOptFit(1111);
      gStyle->SetTitleX(0.25);
      
      TGraphErrors *gOct = new TGraphErrors(enBinMedian.size(), &enBinMedian[0], &RawAsymAndError[0][0], 0, &RawAsymAndError[1][0]);
      std::string title = std::string("Raw Asymmetry A_{SR} ") + itos((int)Elow) + std::string("-") +itos((int)Ehigh) + std::string(" keV Window");
      gOct->SetTitle(title.c_str());
      gOct->SetMarkerStyle(20);
      gOct->SetLineWidth(2);
      gOct->GetXaxis()->SetLimits(0., 800.);
      gOct->GetXaxis()->SetTitle("Energy (keV)");
      gOct->GetYaxis()->SetTitle("Uncorrected Asymmetry");
      gOct->GetXaxis()->CenterTitle();
      gOct->GetYaxis()->CenterTitle();
	
      TF1 *fitOct = new TF1("fitOct","[0]",Elow, Ehigh);
      fitOct->SetLineColor(kRed);
      fitOct->SetLineWidth(3);
      fitOct->SetParameter(0,-0.05);
	
      gOct->Fit("fitOct","R");
	
      std::ofstream ofile(outfilePath.c_str());
      ofile << std::setprecision(10);

      ofile << "RawA_SR " << (fitOct->GetParameter(0)) << " " << fitOct->GetParError(0) << std::endl;
	
      gOct->Draw("AP");
      gOct->SetMinimum(-0.1);
      gOct->SetMaximum(0.);
      c1->Update();
      c1->Print(pdfPath.c_str());
	
	
      delete c1; delete gOct; delete fitOct;
	
      pdfPath = ( basePath2 + (UNBLIND?"UNBLINDED_":"") + corr + "_"  + (withPOL?"withPOL_":"") + "BetaCorrectedAsymmetry_Octet" +
		  itos(octet) + "_AnaCh" + anaChoice + "_" + itos((int)Elow) + std::string("-") +itos((int)Ehigh) +
		  ( key/10==1 ? "_Quadrant"+itos(key-10) : 
		    ( key/10==2 ? "_RadialRing"+itos(key-20) :
		      ( key/10==3 ? "_Strip"+itos(key-30) :"" ) ) ) + ".pdf");
      c1 = new TCanvas("c1", "c1",800, 400);
	
      gOct = new TGraphErrors(enBinMedian.size(), &enBinMedian[0], &AsymAndError[0][0], 0, &AsymAndError[1][0]);
      title = std::string("#frac{2A_{SR}}{#beta}#upoint(1+#Delta_{exp}) ") + itos((int)Elow) + std::string("-") +itos((int)Ehigh) + std::string(" keV Window");
      gOct->SetTitle(title.c_str());
      gOct->SetMarkerStyle(20);
      gOct->SetLineWidth(2);
      gOct->GetXaxis()->SetLimits(0., 800.);
      gOct->GetXaxis()->SetTitle("Energy (keV)");
      gOct->GetYaxis()->SetTitle("Corrected Asymmetry");
      gOct->GetXaxis()->CenterTitle();
      gOct->GetYaxis()->CenterTitle();
	
      fitOct = new TF1("fitOct","[0]",Elow, Ehigh);
      fitOct->SetLineColor(kRed);
      fitOct->SetLineWidth(3);
      fitOct->SetParameter(0,-0.12);
	
      gOct->Fit("fitOct","R");
	
      gOct->Draw("AP");
      gOct->SetMinimum(-0.6);
      gOct->SetMaximum(0.4);
      c1->Update();
      c1->Print(pdfPath.c_str());
	
	
      ofile << "2*A/Beta " << fitOct->GetParameter(0) << " " << fitOct->GetParError(0);
      ofile.close();
	
      delete c1; delete gOct; delete fitOct;
      
    }

    //Now the Quartets
    else if (groupType=="Quartet")
    {	  

      std::string basePath2 =basePath+ "Octet_"+itos(octet)+"/" + groupType + "Asymmetry/"; 

      std::string infilePath[2];
      infilePath[0] = basePath2 + (UNBLIND?"UNBLINDED_":"") + "rawAsymmetry_Octet" + itos(octet) + "_AnaCh" + anaChoice + "_Quartet_A" + ".dat";
      infilePath[1] = basePath2 + (UNBLIND?"UNBLINDED_":"") +  "rawAsymmetry_Octet" + itos(octet) + "_AnaCh" + anaChoice + "_Quartet_B" + ".dat";

      std::string outfilePath[2];
      outfilePath[0] = basePath2 + (UNBLIND?"UNBLINDED_":"") + std::string("UnCorr_") + "FittedAsymmetry_Octet" + itos(octet) + "_AnaCh" + anaChoice + "_Quartet_A" + "_" + itos((int)Elow) + std::string("-") +itos((int)Ehigh) +".dat";
      outfilePath[1] = basePath2 + (UNBLIND?"UNBLINDED_":"") + std::string("UnCorr_") + "FittedAsymmetry_Octet" + itos(octet) + "_AnaCh" + anaChoice + "_Quartet_B" + "_" + itos((int)Elow) + std::string("-") +itos((int)Ehigh) +".dat";

      std::string pdfPathRaw[2];
      pdfPathRaw[0] = basePath2 + (UNBLIND?"UNBLINDED_":"") + std::string("UnCorr_") + "Asymmetry_Octet" + itos(octet) + "_AnaCh" + anaChoice + "_Quartet_A" + "_" + itos((int)Elow) + std::string("-") +itos((int)Ehigh) + ".pdf";
      pdfPathRaw[1] = basePath2 + (UNBLIND?"UNBLINDED_":"") + std::string("UnCorr_") + "Asymmetry_Octet" + itos(octet) + "_AnaCh" + anaChoice + "_Quartet_B" + "_" + itos((int)Elow) + std::string("-") +itos((int)Ehigh) + ".pdf";

      std::string pdfPathCorr[2];
      pdfPathCorr[0] = basePath2 + (UNBLIND?"UNBLINDED_":"") + std::string("UnCorr_") + "BetaCorrectedAsymmetry_Octet" + itos(octet) + "_AnaCh" + anaChoice + "_Quartet_A" + "_" + itos((int)Elow) + std::string("-") +itos((int)Ehigh) + ".pdf";
      pdfPathCorr[1] = basePath2 + (UNBLIND?"UNBLINDED_":"") + std::string("UnCorr_") + "BetaCorrectedAsymmetry_Octet" + itos(octet) + "_AnaCh" + anaChoice + "_Quartet_B" + "_" + itos((int)Elow) + std::string("-") +itos((int)Ehigh) + ".pdf";
      
      for (int quart=0; quart<2; quart++) {

	std::string currentQuart = quart==0?"A":"B";

	//Remove old output files in case they were created on accident and aren't filled with good values
	std::string command = "rm " + outfilePath[quart]; 
	system(command.c_str());

	//First check that Quartet was good
	std::string checkStatus;
	infile.open(infilePath[quart].c_str());
	
	if ( infile.is_open() ) {
	  infile >> checkStatus; 
	  infile.close();
	  if (checkStatus=="BAD") continue;
	}
	else {
	  std::cout << "Could not open binned Asymmetries for Quartet " << currentQuart << " in Octet " << octet << " so skipping...\n";
	  continue; 
	}
	    
	infile.open(infilePath[quart].c_str());
      
	std::vector < std::vector <Double_t > > AsymAndError;
	std::vector < std::vector <Double_t > > RawAsymAndError;
	std::vector < Double_t > enBinMedian;
	
	AsymAndError.resize(2,std::vector <Double_t> (0));
	RawAsymAndError.resize(2,std::vector <Double_t> (0));
	
	Double_t eBinLow, Asym, AsymError;
	
	while (infile >> eBinLow >> Asym >> AsymError) {
	  Double_t En = eBinLow+enBinWidth/2.;
	  enBinMedian.push_back(En);
	  Double_t Beta = returnBeta(En);
	  AsymAndError[0].push_back(Asym*2./Beta);
	  AsymAndError[1].push_back(AsymError*2./Beta);
	  RawAsymAndError[0].push_back(Asym);
	  RawAsymAndError[1].push_back(AsymError);
	}

	infile.close();
	
	TCanvas *c1 = new TCanvas("c1", "c1",800, 400);
	gStyle->SetOptFit(1111);
	gStyle->SetTitleX(0.25);
      
	TGraphErrors *gOct = new TGraphErrors(enBinMedian.size()-5, &enBinMedian[5], &RawAsymAndError[0][5], 0, &RawAsymAndError[1][5]);
	std::string title = std::string("Raw Asymmetry A_{SR} ") + itos((int)Elow) + std::string("-") +itos((int)Ehigh) + std::string(" keV Window");
	gOct->SetTitle(title.c_str());
	gOct->SetMarkerStyle(20);
	gOct->SetLineWidth(2);
	gOct->GetXaxis()->SetLimits(0., 800.);
	gOct->GetXaxis()->SetTitle("Energy (keV)");
	gOct->GetYaxis()->SetTitle("Uncorrected Asymmetry");
	gOct->GetXaxis()->CenterTitle();
	gOct->GetYaxis()->CenterTitle();

	TF1 *fitOct = new TF1("fitOct","[0]",Elow, Ehigh);
	fitOct->SetLineColor(kRed);
	fitOct->SetLineWidth(3);
	fitOct->SetParameter(0,-0.05);
	
	gOct->Fit("fitOct","R");
	
	std::ofstream ofile(outfilePath[quart].c_str());
	ofile << "RawA_SR " << (fitOct->GetParameter(0)) << " " << fitOct->GetParError(0) << std::endl;
      
	gOct->Draw("AP");
	gOct->SetMinimum(-0.1);
	gOct->SetMaximum(0.);
	c1->Update();
	c1->Print(pdfPathRaw[quart].c_str());
	
	delete c1; delete gOct; delete fitOct;
	
	c1 = new TCanvas("c1", "c1",800, 400);
	
	gOct = new TGraphErrors(enBinMedian.size()-5, &enBinMedian[5], &AsymAndError[0][5], 0, &AsymAndError[1][5]);
	title = std::string("#frac{2A_{SR}}{#beta} ") + itos((int)Elow) + std::string("-") +itos((int)Ehigh) + std::string(" keV Window");
	gOct->SetTitle(title.c_str());
	gOct->SetMarkerStyle(20);
	gOct->SetLineWidth(2);
	gOct->GetXaxis()->SetLimits(0., 800.);
	gOct->GetXaxis()->SetTitle("Energy (keV)");
	gOct->GetYaxis()->SetTitle("Uncorrected Asymmetry");
	gOct->GetXaxis()->CenterTitle();
	gOct->GetYaxis()->CenterTitle();
	
        fitOct = new TF1("fitOct","[0]",Elow, Ehigh);
	fitOct->SetLineColor(kRed);
	fitOct->SetLineWidth(3);
	fitOct->SetParameter(0,-0.12);
	
	gOct->Fit("fitOct","R");
	
	gOct->Draw("AP");
	gOct->SetMinimum(-0.6);
	gOct->SetMaximum(0.4);
	c1->Update();
	c1->Print(pdfPathCorr[quart].c_str());

	ofile << "2*A/Beta " << fitOct->GetParameter(0) << " " << fitOct->GetParError(0);
	ofile.close();

	delete c1; delete gOct; delete fitOct;
      }
    }

    //And last of all the Pairs
    else if (groupType=="Pair")
    {
      std::string basePath2 =basePath+ "Octet_"+itos(octet)+"/" + groupType + "Asymmetry/"; 

      for (int pair=0; pair<2; pair++) {
	std::string infilePath[2];
	infilePath[0] = basePath2 + (UNBLIND?"UNBLINDED_":"") + "rawAsymmetry_Octet" + itos(octet) + "_AnaCh" + anaChoice + "_Pair_A" + itos(pair) + ".dat";
	infilePath[1] = basePath2 + (UNBLIND?"UNBLINDED_":"") + "rawAsymmetry_Octet" + itos(octet) + "_AnaCh" + anaChoice + "_Pair_B" + itos(pair) + ".dat";

	std::string outfilePath[2];
	outfilePath[0] = basePath2 + (UNBLIND?"UNBLINDED_":"") + std::string("UnCorr_") + "FittedAsymmetry_Octet" + itos(octet) + "_AnaCh" + anaChoice + "_Pair_A" + itos(pair) + "_" + itos((int)Elow) + std::string("-") +itos((int)Ehigh) + ".dat";
	outfilePath[1] = basePath2 + (UNBLIND?"UNBLINDED_":"") + std::string("UnCorr_") + "FittedAsymmetry_Octet" + itos(octet) + "_AnaCh" + anaChoice + "_Pair_B" + itos(pair) + "_" + itos((int)Elow) + std::string("-") +itos((int)Ehigh) + ".dat";

	std::string pdfPathRaw[2];
	pdfPathRaw[0] = basePath2 + (UNBLIND?"UNBLINDED_":"") + std::string("UnCorr_") + "Asymmetry_Octet" + itos(octet) + "_AnaCh" + anaChoice + "_Pair_A" + itos(pair) + "_" + itos((int)Elow) + std::string("-") +itos((int)Ehigh) + ".pdf";
	pdfPathRaw[1] = basePath2 + (UNBLIND?"UNBLINDED_":"") + std::string("UnCorr_") + "Asymmetry_Octet" + itos(octet) + "_AnaCh" + anaChoice + "_Pair_B" + itos(pair) + "_" + itos((int)Elow) + std::string("-") +itos((int)Ehigh) + ".pdf";
	
	std::string pdfPathCorr[2];
	pdfPathCorr[0] = basePath2 + (UNBLIND?"UNBLINDED_":"") + std::string("UnCorr_") + "BetaCorrectedAsymmetry_Octet" + itos(octet) + "_AnaCh" + anaChoice + "_Pair_A" + itos(pair) + "_" + itos((int)Elow) + std::string("-") +itos((int)Ehigh) + ".pdf";
	pdfPathCorr[1] = basePath2 + (UNBLIND?"UNBLINDED_":"") + std::string("UnCorr_") + "BetaCorrectedAsymmetry_Octet" + itos(octet) + "_AnaCh" + anaChoice + "_Pair_B" + itos(pair) + "_" + itos((int)Elow) + std::string("-") +itos((int)Ehigh) + ".pdf";
	

	//Remove old output files in case they were created on accident and aren't filled with good values
	std::string command = "rm " + outfilePath[0]; 
	system(command.c_str());
	command = "rm " + outfilePath[1]; 
	system(command.c_str());
	
	for (int quart=0; quart<2; quart++) {
	  
	  std::string currentPair = (quart==0?"A":"B") + itos(pair);	  

	  //First check that pair was good
	  std::string checkStatus;
	  infile.open(infilePath[quart].c_str());
	  
	  if ( infile.is_open() ) {
	    infile >> checkStatus; 
	    infile.close();
	    if (checkStatus=="BAD") continue;
	  }
	  else {
	    std::cout << "Could not open binned Asymmetries for Pair " << currentPair << " in Octet " << octet << " so skipping...\n";
	    continue; 
	  }


	 
	  
	  infile.open(infilePath[quart].c_str());
	  
	  std::vector < std::vector <Double_t > > AsymAndError;
	  std::vector < std::vector <Double_t > > RawAsymAndError;
	  std::vector < Double_t > enBinMedian;
	  
	  AsymAndError.resize(2,std::vector <Double_t> (0));
	  RawAsymAndError.resize(2,std::vector <Double_t> (0));
	  
	  Double_t eBinLow, Asym, AsymError;
	
	  while (infile >> eBinLow >> Asym >> AsymError) {
	    Double_t En = eBinLow+enBinWidth/2.;
	    enBinMedian.push_back(En);
	    Double_t Beta = returnBeta(En);
	    AsymAndError[0].push_back(Asym*2./Beta);
	    AsymAndError[1].push_back(AsymError*2./Beta);
	    RawAsymAndError[0].push_back(Asym);
	    RawAsymAndError[1].push_back(AsymError);
	  }

	  infile.close();
	  
	  
	  TCanvas *c1 = new TCanvas("c1", "c1",800, 400);
	  gStyle->SetOptFit(1111);
	  gStyle->SetTitleX(0.25);
      
	  TGraphErrors *gOct = new TGraphErrors(enBinMedian.size()-5, &enBinMedian[5], &RawAsymAndError[0][5], 0, &RawAsymAndError[1][5]);
	  std::string title = std::string("Raw Asymmetry A_{SR} ") + itos((int)Elow) + std::string("-") +itos((int)Ehigh) + std::string(" keV Window");
	  gOct->SetTitle(title.c_str());
	  gOct->SetMarkerStyle(20);
	  gOct->SetLineWidth(2);
	  gOct->GetXaxis()->SetLimits(0., 800.);
	  gOct->GetXaxis()->SetTitle("Energy (keV)");
	  gOct->GetYaxis()->SetTitle("Uncorrected Asymmetry");
	  gOct->GetXaxis()->CenterTitle();
	  gOct->GetYaxis()->CenterTitle();

	  TF1 *fitOct = new TF1("fitOct","[0]",Elow, Ehigh);
	  fitOct->SetLineColor(kRed);
	  fitOct->SetLineWidth(3);
	  fitOct->SetParameter(0,-0.05);
	  
	  gOct->Fit("fitOct","R");
	  
	  std::ofstream ofile(outfilePath[quart].c_str());
	  ofile << "RawA_SR " << (fitOct->GetParameter(0)) << " " << fitOct->GetParError(0) << std::endl;
	  
	  gOct->Draw("AP");
	  gOct->SetMinimum(-0.1);
	  gOct->SetMaximum(0.);
	  c1->Update();
	  c1->Print(pdfPathRaw[quart].c_str());
	  
	  delete c1; delete gOct; delete fitOct;
	  
	  c1 = new TCanvas("c1", "c1",800, 400);
	  
	  gOct = new TGraphErrors(enBinMedian.size()-5, &enBinMedian[5], &AsymAndError[0][5], 0, &AsymAndError[1][5]);
	  title = std::string("#frac{2A_{SR}}{#beta} ") + itos((int)Elow) + std::string("-") +itos((int)Ehigh) + std::string(" keV Window");
	  gOct->SetTitle(title.c_str());
	  gOct->SetMarkerStyle(20);
	  gOct->SetLineWidth(2);
	  gOct->GetXaxis()->SetLimits(0., 800.);
	  gOct->GetXaxis()->SetTitle("Energy (keV)");
	  gOct->GetYaxis()->SetTitle("Uncorrected Asymmetry");
	  gOct->GetXaxis()->CenterTitle();
	  gOct->GetYaxis()->CenterTitle();
	  
	  fitOct = new TF1("fitOct","[0]",Elow, Ehigh);
	  fitOct->SetLineColor(kRed);
	  fitOct->SetLineWidth(3);
	  fitOct->SetParameter(0,-0.12);
	  
	  gOct->Fit("fitOct","R");
	  
	  gOct->Draw("AP");
	  gOct->SetMinimum(-0.6);
	  gOct->SetMaximum(0.4);
	  c1->Update();
	  c1->Print(pdfPathCorr[quart].c_str());

	  ofile << "2*A/Beta " << fitOct->GetParameter(0) << " " << fitOct->GetParError(0);
	  ofile.close();
	  
	  delete c1; delete gOct; delete fitOct;
	}
      }
    }
  }
};
      


//These are summed over the energy range given as Elow-Ehigh

void ProduceRawAsymmetries(Int_t octBegin, Int_t octEnd, std::string anaChoice, Double_t Elow, Double_t Ehigh, Double_t enBinWidth, bool UKdata, bool simulation) {
 
    //Looking at octet evts and asymmetries
  std::string basePath = simulation ? getenv("SIM_ANALYSIS_RESULTS") : UKdata ? getenv("ANALYSIS_RESULTS") : getenv("MPM_ANALYSIS_RESULTS");
  
  std::string pairFile = basePath + std::string("Asymmetries/Pair_RawAsymmetries_AnaCh") + anaChoice 
    + std::string("_") + itos((int)Elow) + std::string("-") + itos((int)Ehigh) +std::string(".dat");
  std::string quartetFile = basePath + std::string("Asymmetries/Quartet_RawAsymmetries_AnaCh") + anaChoice
    + std::string("_") + itos((int)Elow) + std::string("-") + itos((int)Ehigh) +std::string(".dat");
  std::string octetFile = basePath + std::string("Asymmetries/Octet_RawAsymmetries_AnaCh") + anaChoice 
    + std::string("_") + itos((int)Elow) + std::string("-") + itos((int)Ehigh) +std::string(".dat");

    
    //unsigned int octetNum = 1;
  std::ofstream octAsym(octetFile.c_str());
  std::ofstream quartAsym(quartetFile.c_str());
  std::ofstream pairAsym(pairFile.c_str());
 
  Int_t numPair=0, numQuart=0;
  std::vector < std::vector <Double_t > > octetRawAsymAndError;
  std::vector < std::vector <Double_t > > quartetRawAsymAndError;
  std::vector < std::vector <Double_t > > pairRawAsymAndError;
  octetRawAsymAndError.resize(3,std::vector <Double_t> (0));
  quartetRawAsymAndError.resize(3,std::vector <Double_t> (0));
  pairRawAsymAndError.resize(3,std::vector <Double_t> (0)); 
  
  Double_t Asym=0., AsymError=0.;

  for (Int_t octet=octBegin; octet<=octEnd; octet++) {

    if (std::find(badOct.begin(), badOct.end(),octet) != badOct.end()) { numPair+=4; numQuart+=2; continue; } //Checking if octet should be ignored for data quality reasons

    try {
      //OCTETS
      OctetAsymmetry oct(octet,anaChoice,enBinWidth, 50., UKdata, simulation);     
      oct.calcTotalAsymmetry(Elow,Ehigh);
      //oct.calcAsymmetryBinByBin(1);
      //oct.calcTotalAsymmetry(170.,630.,1);
      
      Asym = oct.returnTotalAsymmetry();
      AsymError = oct.returnTotalAsymmetryError();
      octetRawAsymAndError[0].push_back(octet);
      octetRawAsymAndError[1].push_back(Asym);
      octetRawAsymAndError[2].push_back(AsymError);

      std::cout << "Asymmetry for Octet " << octet << ":\n";
      std::cout << Asym << " +- " << AsymError << std::endl;
      octAsym << octet << " " << Asym << " " << AsymError << "\n";
     
    }
    catch(const char* ex){
      std::cerr << "Error: " << ex << std::endl;
    }
    
    try {
      // QUARTETS
      QuartetAsymmetry quart(octet,anaChoice,enBinWidth, 50., UKdata, simulation); 
      quart.calcTotalAsymmetry(Elow,Ehigh);
	//quart.calcAsymmetryBinByBin(1);
	//quart.calcTotalAsymmetry(170.,630.,1);
      
      if (quart.boolGoodQuartet(0)) {
	
	Asym = quart.returnTotalAsymmetry_QuartetA();
	AsymError = quart.returnTotalAsymmetryError_QuartetA();
	quartetRawAsymAndError[0].push_back(numQuart);
	quartetRawAsymAndError[1].push_back(Asym);
	quartetRawAsymAndError[2].push_back(AsymError);
  
	std::cout << "Asymmetry for Quartet " << numQuart << ":\n";
	std::cout << Asym << " +- " << AsymError << std::endl;
	quartAsym << numQuart << " " << Asym << " " << AsymError << "\n";
      }
      numQuart++;

      if (quart.boolGoodQuartet(1)) {
	Asym = quart.returnTotalAsymmetry_QuartetB();
	AsymError = quart.returnTotalAsymmetryError_QuartetB();
	quartetRawAsymAndError[0].push_back(numQuart);
	quartetRawAsymAndError[1].push_back(Asym);
	quartetRawAsymAndError[2].push_back(AsymError);

	std::cout << "Asymmetry for Quartet " << numQuart << ":\n";
	std::cout << Asym << " +- " << AsymError << std::endl;
	quartAsym << numQuart << " " << Asym << " " << AsymError << "\n";
      }
      numQuart++;
    }
    catch(const char* ex){
      numQuart+=2;
      std::cerr << "Error: " << ex << std::endl;
    }

    try {
      // PAIRS
      PairAsymmetry pair(octet,anaChoice,enBinWidth, 50., UKdata, simulation); 
      pair.calcTotalAsymmetry(Elow,Ehigh);
	//pair.calcAsymmetryBinByBin(1);
	//pair.calcTotalAsymmetry(170.,630.,1);
      
      if (pair.boolGoodPair(0,0)) {
	
	Asym = pair.returnTotalAsymmetry_PairA0();
	AsymError = pair.returnTotalAsymmetryError_PairA0();
	pairRawAsymAndError[0].push_back(numPair);
	pairRawAsymAndError[1].push_back(Asym);
	pairRawAsymAndError[2].push_back(AsymError);
	
	std::cout << "Asymmetry for Pair " << numPair << ":\n";
	std::cout << Asym << " +- " << AsymError << std::endl;
	pairAsym << numPair << " " << Asym << " " << AsymError << "\n";
      }
      numPair++;

      if (pair.boolGoodPair(0,1)) {
	Asym = pair.returnTotalAsymmetry_PairA1();
	AsymError = pair.returnTotalAsymmetryError_PairA1();
	pairRawAsymAndError[0].push_back(numPair);
	pairRawAsymAndError[1].push_back(Asym);
	pairRawAsymAndError[2].push_back(AsymError);

	std::cout << "Asymmetry for Pair " << numPair << ":\n";
	std::cout << Asym << " +- " << AsymError << std::endl;
	pairAsym << numPair << " " << Asym << " " << AsymError << "\n";
      }
      numPair++;

      if (pair.boolGoodPair(1,0)) {
	
	Asym = pair.returnTotalAsymmetry_PairB0();
	AsymError = pair.returnTotalAsymmetryError_PairB0();
	pairRawAsymAndError[0].push_back(numPair);
	pairRawAsymAndError[1].push_back(Asym);
	pairRawAsymAndError[2].push_back(AsymError);

	std::cout << "Asymmetry for Pair " << numPair << ":\n";
	std::cout << Asym << " +- " << AsymError << std::endl;
	pairAsym << numPair << " " << Asym << " " << AsymError << "\n";
      }
      numPair++;

      if (pair.boolGoodPair(1,1)) {
	Asym = pair.returnTotalAsymmetry_PairB1();
	AsymError = pair.returnTotalAsymmetryError_PairB1();
	pairRawAsymAndError[0].push_back(numPair);
	pairRawAsymAndError[1].push_back(Asym);
	pairRawAsymAndError[2].push_back(AsymError);

	std::cout << "Asymmetry for Pair " << numPair << ":\n";
	std::cout << Asym << " +- " << AsymError << std::endl;
	pairAsym << numPair << " " << Asym << " " << AsymError << "\n";
      }
      numPair++;
    }
    catch(const char* ex){
      numPair+=4;
      std::cerr << "Error: " << ex << std::endl;
    }
  }

  octAsym.close();
  quartAsym.close();
  pairAsym.close();

  std::string pdfFile = basePath + std::string("Asymmetries/RawAsymmetries_AnaCh") + anaChoice 
    + std::string("_") + itos((int)Elow) + std::string("-") + itos((int)Ehigh) +std::string(".pdf");

  TCanvas *c1 = new TCanvas("c1", "c1", 1000., 1000.);
  c1->Divide(1,3);

  gStyle->SetOptFit(1111);
  gStyle->SetTitleX(0.25);
 
  std::vector <double> errorX(pairRawAsymAndError[0].size(),0.);
  std::string title;

  c1->cd(1);

  TGraphErrors *gOct = new TGraphErrors(octetRawAsymAndError[0].size(), &octetRawAsymAndError[0][0],&octetRawAsymAndError[1][0],&errorX[0], &octetRawAsymAndError[2][0]);
  title = std::string("Octet-By-Octet Raw Measured Asymmetry ") + itos((int)Elow) + std::string("-") +itos((int)Ehigh) + std::string(" keV Window");
  gOct->SetTitle(title.c_str());
  gOct->SetMarkerStyle(20);
  gOct->SetLineWidth(2);
  gOct->GetXaxis()->SetLimits(-2., octetRawAsymAndError[0][octetRawAsymAndError[0].size()-1]+2.);
  gOct->GetXaxis()->SetTitle("Octet Number");
  gOct->GetYaxis()->SetTitle("Raw Asymmetry");
  gOct->GetXaxis()->CenterTitle();
  gOct->GetYaxis()->CenterTitle();
  
  TF1 *fitOct = new TF1("fitOct","[0]",octetRawAsymAndError[0][0]-1., octetRawAsymAndError[0][octetRawAsymAndError[0].size()-1]+1.);
  fitOct->SetLineColor(kRed);
  fitOct->SetLineWidth(3);
  fitOct->SetParameter(0,0.05);

  gOct->Fit("fitOct","R");
  
  gOct->Draw("AP");
  gOct->SetMinimum(0.03);
  gOct->SetMaximum(0.07);
  c1->Update();

  c1->cd(2);

  TGraphErrors *gQuart = new TGraphErrors(quartetRawAsymAndError[0].size(), &quartetRawAsymAndError[0][0],&quartetRawAsymAndError[1][0],&errorX[0], &quartetRawAsymAndError[2][0]);
  title = std::string("Quartet-By-Quartet Raw Measured Asymmetry ") + itos((int)Elow) + std::string("-") +itos((int)Ehigh) + std::string(" keV Window");
  gQuart->SetTitle(title.c_str());
  gQuart->SetMarkerStyle(20);
  gQuart->SetLineWidth(2);
  gQuart->GetXaxis()->SetLimits(-2., quartetRawAsymAndError[0][quartetRawAsymAndError[0].size()-1]+2.);
  gQuart->GetXaxis()->SetTitle("Quartet Number");
  gQuart->GetYaxis()->SetTitle("Raw Asymmetry");
  gQuart->GetXaxis()->CenterTitle();
  gQuart->GetYaxis()->CenterTitle();
  
  TF1 *fitQuart = new TF1("fitQuart","[0]",quartetRawAsymAndError[0][0]-1., quartetRawAsymAndError[0][quartetRawAsymAndError[0].size()-1]+1.);
  fitQuart->SetLineColor(kRed);
  fitQuart->SetLineWidth(3);
  fitQuart->SetParameter(0,0.05);


  gQuart->Fit("fitQuart","R");
  
  gQuart->Draw("AP");
  gQuart->SetMinimum(0.03);
  gQuart->SetMaximum(0.07);
  c1->Update();

  c1->cd(3);

  TGraphErrors *gPair = new TGraphErrors(pairRawAsymAndError[0].size(), &pairRawAsymAndError[0][0],&pairRawAsymAndError[1][0],&errorX[0], &pairRawAsymAndError[2][0]);
  title = std::string("Pair-By-Pair Raw Measured Asymmetry ") + itos((int)Elow) + std::string("-") +itos((int)Ehigh) + std::string(" keV Window");
  gPair->SetTitle(title.c_str());
  gPair->SetMarkerStyle(20);
  gPair->SetLineWidth(2);
  gPair->GetXaxis()->SetLimits(-2., pairRawAsymAndError[0][pairRawAsymAndError[0].size()-1]+2.);
  gPair->GetXaxis()->SetTitle("Pair Number");
  gPair->GetYaxis()->SetTitle("Raw Asymmetry");
  gPair->GetXaxis()->CenterTitle();
  gPair->GetYaxis()->CenterTitle();
  
  TF1 *fitPair = new TF1("fitPair","[0]",pairRawAsymAndError[0][0]-1., pairRawAsymAndError[0][pairRawAsymAndError[0].size()-1]+1.);
  fitPair->SetLineColor(kRed);
  fitPair->SetLineWidth(3);
  fitPair->SetParameter(0,0.05);


  gPair->Fit("fitPair","R");
  
  gPair->Draw("AP");
  gPair->SetMinimum(0.03);
  gPair->SetMaximum(0.07);
  c1->Update();
  
  c1->Print(pdfFile.c_str());

  delete gPair; delete gQuart; delete gOct; delete fitPair; delete fitQuart; delete fitOct; delete c1;
																     

};

/////////////////////// Corrections ////////////////////////////

std::vector <Double_t> LoadTheoryCorrections(std::vector <Double_t> enBinMidpoint) {
  std::vector <Double_t> syst(enBinMidpoint.size(), 1.);

  if (corr.size()<7) return syst;
  if ( corr.substr(0,7)!=std::string("AllCorr") && corr.substr(0,7)!=std::string("DeltaTh") ) return syst;

  std::string c = corr.substr(0,7)==std::string("DeltaTh")?corr.substr(11,3):"ALL"; // options are DeltaTheoryRecoil,DeltaTheoryRadiative,DeltaTheoryALL

  double WE=0.;
  for (UInt_t i=0; i<syst.size(); i++) {
    WE = (enBinMidpoint[i]+m_e)/m_e;
    if (c==std::string("Rec") || c==std::string("ALL") ) syst[i] *= (1.+WilkinsonACorrection(WE));
    if (c==std::string("Rad") || c==std::string("ALL") ) syst[i] *= (1.+shann_h_minus_g_a2pi(WE));
  }

  return syst;
};

std::vector <Double_t> LoadAngleCorrections(std::vector <Double_t> enBinMidpoint,Int_t oct, std::string anaCh) {

  bool useDelta3i = true;

  std::vector <Double_t> syst(enBinMidpoint.size(), 1.);
  if (corr.size()<7) return syst;
  if ( corr.substr(0,7)!=std::string("AllCorr") && corr.substr(0,7)!=std::string("DeltaAn") ) return syst;

  std::string c;
  TString filename;
  
  if ( useDelta3i ) { //anaCh==std::string("C") ) {
    c = corr.substr(0,7)==std::string("DeltaAn")?corr.substr(10,1):"ALL";
    filename = TString::Format("../systematics/AngleCorrections/%s_delta_3%s_anaCh%s.txt",
			       oct<60?"2011-2012":"2012-2013",
			       (c==std::string("0")?"0":
				(c==std::string("1")?"1":
				 (c==std::string("2")?"2":
				  (c==std::string("3")?"3":"")))),anaCh.c_str());
  } else {
    filename = TString::Format("../systematics/AngleCorrections/%s_DeltaAngle_anaCh%s.txt",oct<60?"2011-2012":"2012-2013",anaCh.c_str());
  }

  std::cout << filename << std::endl;                                                                                                 

  std::ifstream infile(filename.Data());

  if (!infile.is_open()) throw "Couldn't open file in LoadAngleCorrections!";

  //  std::cout << "\n";                                                                                                                       

  //Read in the systematics                                                                                                                    
  Double_t mid; //midpoint of bin                                                                                                              
  std::string hold1;
  std::string hold2;  //systematics correction (Apure/Aproc)                                                                                     

  Int_t it = 0;

  while (infile >> mid >> hold1 >> hold2) {
    if (mid==enBinMidpoint[it]) {
      if ( useDelta3i ) {//anaCh==std::string("C") ) {
	syst[it] = (hold1!=std::string("nan") && hold1!=std::string("-nan")) ? atof(hold1.c_str())+1. : 1.;
	std::cout << mid << " " << atof(hold1.c_str())+1. << "\n";
      } else {
	syst[it] = (hold1!=std::string("inf") && hold2!=std::string("inf")) ? atof(hold2.c_str())+1. : 1.;
	std::cout << mid << " " << atof(hold2.c_str())+1. << "\n";
      }
      it++;
    }
  }
    
  return syst;
  

};



std::vector <Double_t> LoadBackscCorrections(std::vector <Double_t> enBinMidpoint, Int_t oct, std::string anaCh) {
  std::vector <Double_t> syst(enBinMidpoint.size(), 1.);
  if (corr.size()<7) return syst;
  if ( corr.substr(0,7)!=std::string("AllCorr") && corr.substr(0,7)!=std::string("DeltaBa") ) return syst;

  std::string c = corr.substr(0,7)==std::string("DeltaBa")?corr.substr(11,1):"ALL";
  TString filename = TString::Format("../systematics/OldMCCorrection/%s_delta_2%s_anaCh%s.txt",
				     oct<60?"2011-2012":"2012-2013",
				     (c==std::string("0")?"0":
				      c==std::string("1")?"1":
				      c==std::string("2")?"2":
				      c==std::string("3")?"3":""),
				     anaCh.c_str());
  //std::cout << filename.Data() << std::endl;
  
  std::ifstream infile(filename.Data());

  if (!infile.is_open()) throw "Couldn't open file in LoadFieldDipCorrections!";

  //  std::cout << "\n";

  //Read in the systematics
  Double_t mid; //midpoint of bin
  std::string sys;  //systematics correction (Apure/Aproc)
  std::string hold;

  Int_t it = 0;

  while (infile >> mid >> sys >> hold) {
    if (mid==enBinMidpoint[it]) {
      syst[it] = (sys!=std::string("inf") && sys!=std::string("-nan") && sys!=std::string("nan")) ? atof(sys.c_str())+1. : 1.;
      std::cout << mid << " " << syst[it]<< "\n";
      it++;
    }
  }

  return syst;


};

/*std::vector < std::vector <Double_t> > LoadOctetSystematics(Int_t octet, std::string anaChoice, std::vector <Double_t> enBinMidpoint) {

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
*/
