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

int main()
{

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
      
      //delete oct;*/
      
      
      
      /* oct->calcBGsubtractedEvts();
      for (unsigned int type=0;type<4;type++) {
	evtVecHold = oct->getNumBGsubtrEvts(enWinLow,enWinHigh,type);
	evts[type] = evtVecHold[0]+evtVecHold[1];
	totalEvts[type]+=evts[type];
      }
      
      std::cout << "Type0: " << evts[0] 
		<<"\nType1: " << evts[1]
		<<"\nType2: " << evts[2]
		<<"\nType3: " << evts[3]
		<< std::endl;*/
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
  return 0;
  
  /*std::string dbAddress = "localhost";
    std::string dbname = "UCNADB_full_quad1";
    std::string port = "3306";
    std::string dbAddressFull = "mysql://"+dbAddress+":"+port+"/"+dbname;
    std::string dbUser = "ucn";
    std::string dbPass = "T1td4ctc";
    try {
    SQLdatabase *db = new SQLdatabase(dbname,dbAddress,dbUser,dbPass);
    db->query = "Select * from posmap_set where posmap_set_id=283";
    db->fetchQuery();
    std::string result = db->returnQueryEntry();
    std::cout << result << std::endl;
    }
    catch(const char* ex){
    std::cerr << "Error: " << ex << std::endl;
    }*/

  
}

///// Snippet for output of rate histograms   
/*
  char temp[200];
  if (!outputFile) {
    outputFile = std::string(getenv("groupAnaDir"))+"/"+std::string(sprintf(temp,"%i_rawRate.root",runNumber));
  }
  file = new TFile(outputFile.c_str(),"RECREATE");
*/


