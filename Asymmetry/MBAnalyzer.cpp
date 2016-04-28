/*
Compilable code which houses the heart of the analysis. 
Here are classes which handle calculating asymmetries of either data 
or simulation. There will be flags to be passed dictating what type of 
asymmetry, what type of data, and what octet/runs to use. 

The code should be able to do BG subtraction, construct asymmetries, 
or do both depending on the flags passed.

Maybe even add in writing the final answer to the database if the user wants to

*/

#include "MBAnalyzer.hh"
#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <cmath>
#include <fstream>

void OctetProcessor(Int_t octBegin, Int_t octEnd, Double_t enBinWidth, bool simulation, bool UKdata);
void QuartetProcessor(Int_t octBegin, Int_t octEnd, Double_t enBinWidth, bool simulation, bool UKdata);
void PairProcessor(Int_t octBegin, Int_t octEnd, Double_t enBinWidth, bool simulation, bool UKdata);

void OctetRawAsymmetry(Int_t octBegin, Int_t octEnd, Double_t enBinWidth, bool simulation, bool UKdata);
void QuartetRawAsymmetry(Int_t octBegin, Int_t octEnd, Double_t enBinWidth, bool simulation, bool UKdata);
void PairRawAsymmetry(Int_t octBegin, Int_t octEnd, Double_t enBinWidth, bool simulation, bool UKdata);

int main()
{
  double enWinLow = 220.; //170
  double enWinHigh = 680.; //630
  
  PairAsymmetry *pair;
  for (unsigned int octet=0; octet<=0;octet++) {
   
    try {
      pair = new PairAsymmetry(octet,10.,50.,true, false);
      
      pair->calcTotalAsymmetry(enWinLow,enWinHigh,1);
      std::cout<< "Finished Total Asymmetry"<<std::endl;
      pair->calcAsymmetryBinByBin(1);
      std::cout<<std::endl;
      pair->calcSuperSum(1);
      std::cout<<std::endl;
      //oct->calcTotalAsymmetry(170.,630.,1);
      
      double AsymA0 = pair->returnTotalAsymmetry_PairA0();
      double AsymErrorA0 = pair->returnTotalAsymmetryError_PairA0();
      std::cout << "Asymmetry for Pair A0 in " << octet << ":\n";
      std::cout << AsymA0 << " +- " << AsymErrorA0 << std::endl;
      
      double AsymA1 = pair->returnTotalAsymmetry_PairA1();
      double AsymErrorA1 = pair->returnTotalAsymmetryError_PairA1();
      std::cout << "Asymmetry for Pair A1 in " << octet << ":\n";
      std::cout << AsymA1 << " +- " << AsymErrorA1 << std::endl;

      double AsymB0 = pair->returnTotalAsymmetry_PairB0();
      double AsymErrorB0 = pair->returnTotalAsymmetryError_PairB0();
      std::cout << "Asymmetry for Pair B0 in " << octet << ":\n";
      std::cout << AsymB0 << " +- " << AsymErrorB0 << std::endl;
      
      double AsymB1 = pair->returnTotalAsymmetry_PairB1();
      double AsymErrorB1 = pair->returnTotalAsymmetryError_PairB1();
      std::cout << "Asymmetry for Pair B1 in " << octet << ":\n";
      std::cout << AsymB1 << " +- " << AsymErrorB1 << std::endl;
     
      delete pair;
    }
    catch(const char* ex){
      std::cerr << "Error: " << ex << std::endl;
    }
  }
 
  return 0;
}


void OctetProcessor(Int_t octBegin, Int_t octEnd, Double_t enBinWidth=10., bool simulation=false, bool UKdata=true) {
  

};

void QuartetProcessor(Int_t octBegin, Int_t octEnd, Double_t enBinWidth=10., bool simulation=false, bool UKdata=true) {
  double enWinLow = 220.; //170
  double enWinHigh = 680.; //630
  
  QuartetAsymmetry *quart;
  for (unsigned int octet=0; octet<=0;octet++) {
      quart = new QuartetAsymmetry(octet,10.,50.,true, true);
      
      quart->calcTotalAsymmetry(enWinLow,enWinHigh,1);
      std::cout<< "Finished Total Asymmetry"<<std::endl;
      quart->calcAsymmetryBinByBin(1);
      std::cout<<std::endl;
      quart->calcSuperSum(1);
      std::cout<<std::endl;
      //oct->calcTotalAsymmetry(170.,630.,1);
      
      double AsymA = quart->returnTotalAsymmetry_QuartetA();
      double AsymErrorA = quart->returnTotalAsymmetryError_QuartetA();
      std::cout << "Asymmetry for Quartet A in " << octet << ":\n";
      std::cout << AsymA << " +- " << AsymErrorA << std::endl;

      double AsymB = quart->returnTotalAsymmetry_QuartetB();
      double AsymErrorB = quart->returnTotalAsymmetryError_QuartetB();
      std::cout << "Asymmetry for Quaret B in " << octet << ":\n";
      std::cout << AsymB << " +- " << AsymErrorB << std::endl;
      
      delete quart;
  }
  
};


void OctetRawAsymmetry() {
  std::string inDir = std::string(getenv("REPLAY_PASS4"));
 
    

    //Looking at octet evts and asymmetries

    //unsigned int octetNum = 1;
  ofstream rawAsym("rawAsymmetryByOctet_0-59_AnaCh1_50mm.dat");
    
  double enWinLow = 220.; //170
  double enWinHigh = 680.; //630
  std::vector <double> evts(4,0.);
  std::vector <double> totalEvts(4,0.);
  std::vector <double> evtVecHold; //Holds the events returned from Asymmetries.cc
  OctetAsymmetry *oct;
  for (unsigned int octet=0; octet<=59;octet++) {
    rawAsym << octet << " ";
    try {
      oct = new OctetAsymmetry(octet,10.,50.,true);
      
      oct->calcTotalAsymmetry(enWinLow,enWinHigh,1);
      std::cout<<std::endl;
      //oct->calcAsymmetryBinByBin(1);
      std::cout<<std::endl;
      //oct->calcTotalAsymmetry(170.,630.,1);
      
      double Asym = oct->returnTotalAsymmetry();
      double AsymError = oct->returnTotalAsymmetryError();
      std::cout << "Asymmetry for Octet " << octet << ":\n";
      std::cout << Asym << " +- " << AsymError << std::endl;
      rawAsym << Asym << " " << AsymError << "\n";
      
      delete oct;
    }
    catch(const char* ex){
      rawAsym << "-nan -nan\n";
      std::cerr << "Error: " << ex << std::endl;
    }
  }
  
  std::cout << "Total Event counts for energy window " << enWinLow << " - " << enWinHigh << ":\n"
	    <<"Type0: " << totalEvts[0] 
	    <<"\nType1: " << totalEvts[1]
	    <<"\nType2: " << totalEvts[2]
	    <<"\nType3: " << totalEvts[3]
	    <<"\nTotal Betas: " << totalEvts[0]+totalEvts[1]+totalEvts[2]+totalEvts[3]
	    << std::endl; 

   rawAsym.close();

};
  
  
