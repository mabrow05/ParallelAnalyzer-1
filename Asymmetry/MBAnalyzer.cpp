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
  int run = 17126;

  std::string inDir = std::string(getenv("REPLAY_PASS4"));
  try {
    /*std::cout << "Run " << run << std::endl;
    EvtRateHandler *rt = new EvtRateHandler(run,inDir);
    rt->polarization(run);
    std::cout << rt->pol  << std::endl;
    rt->CalcRates(10,60.);
    TFile *f = new TFile("test.root","RECREATE");
    TH1D hisE = rt->getRateHist(0,0);
    TH1D hisW = rt->getRateHist(1,0);
    f->Write();
    f->Close();
    
    std::vector< std::vector <double> > evec = rt->getRateVectors(0);
    std::vector< std::vector <double> > wvec = rt->getRateVectors(1);

    double timeE = rt->returnRunLength(0);
    double timeW = rt->returnRunLength(1);
    delete rt;

    std::cout << timeE << "\t" << timeW << std::endl << "/////////////////////////////////////////\n";;
      //for (unsigned int i=0; i<evec[0].size(); i++) {
      //std::cout << evec[0][i] << " " << wvec[0][i] << std::endl;
      //}
      
    //testing BG subtracted rate

    BGSubtractedRate *bg = new BGSubtractedRate(run,10.,60.,true,false);

    std::cout << "initialized BGStubtractedRate\n";
    std::cout << bg->getBackgroundRun(17126) << std::endl;
    bg->calcBGSubtRates();
    std::vector <double> evecbg = bg->returnBGSubtRate(0,0);
    std::vector <double> wvecbg = bg->returnBGSubtRate(1,0);

    std::vector <double> evecbgErr = bg->returnBGSubtRateError(0,0);
    std::vector <double> wvecbgErr = bg->returnBGSubtRateError(1,0);

    std::vector<double> runLengthBeta = bg->returnRunLengths(true);
    std::vector<double> runLengthBG = bg->returnRunLengths(false);

    for (unsigned int i=0; i<evecbg.size(); i++) {
      std::cout << evecbg[i] << " " << evecbgErr[i] << " " << wvecbg[i] << " " << wvecbgErr[i] << std::endl;
    }
    
    std::cout << "RunLength: E      W\n" 
    	      << runLengthBeta[0] << " " << runLengthBeta[1] << std::endl
    	      << runLengthBG[0] << " " << runLengthBG[1] << std::endl;

    delete bg;*/

    //Looking at octet evts and asymmetries

    //unsigned int octetNum = 1;
    ofstream rawAsym("rawAsymmetryByOctet_0-59_AnaCh7_50mm.dat");
    
    double enWinLow = 170.; //170
    double enWinHigh = 630.; //630
    double evts[4]={0.};
    double totalEvts[4]={0.};
    OctetAsymmetry *oct;
    for (unsigned int octet=0; octet<60;octet++) {
      oct = new OctetAsymmetry(octet,10.,50.,true);
      oct->loadRates();
      oct->calcTotalAsymmetry(enWinLow,enWinHigh,7);
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
    rawAsym.close();
	
      /*oct->calcBGsubtractedEvts();
      for (unsigned int type=0;type<4;type++) {
	evts[type] = oct->getNumBGsubtrEvts(enWinLow,enWinHigh,type);
	totalEvts[type]+=evts[type];
      }
	  
      std::cout << "Type0: " << evts[0] 
		<<"\nType1: " << evts[1]
		<<"\nType2: " << evts[2]
		<<"\nType3: " << evts[3]
		<< std::endl;
      delete oct;
    }
    
    std::cout << "Total Event counts for energy window " << enWinLow << " - " << enWinHigh << ":\n"
	      <<"Type0: " << totalEvts[0] 
	      <<"\nType1: " << totalEvts[1]
	      <<"\nType2: " << totalEvts[2]
	      <<"\nType3: " << totalEvts[3]
	      <<"\nTotal Betas: " << totalEvts[0]+totalEvts[1]+totalEvts[2]+totalEvts[3]
	      << std::endl;
      */
      
  }
  catch(const char* ex){
    std::cerr << "Error: " << ex << std::endl;
  }
    
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


