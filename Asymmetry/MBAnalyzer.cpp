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

int main()
{
  int run = 17800;
  std::string inDir = std::string(getenv("UCNAOUTPUTDIR"))+"/hists";
  try {
    EvtRateHandler *rt = new EvtRateHandler(run,inDir);
  
    std::cout << rt->polarization(run) << std::endl;
    rt->CalcRates(0,50,50.);
    TFile *f = new TFile("test.root","RECREATE");
    TH1D hisE = rt->getRateHist(0);
    TH1D hisW = rt->getRateHist(1);
    f->Write();
    f->Close();
    
    std::vector <double> evec = rt->getRateVector(0);
    std::vector <double> wvec = rt->getRateVector(1);

    delete rt;

    for (unsigned int i=0; i<evec.size(); i++) {
      std::cout << evec[i] << " " << wvec[i] << std::endl;
    }

    //testing BG subtracted rate

    BGSubtractedRate *bg = new BGSubtractedRate(run,50., 50.,0);

    std::cout << "initialized BGStubtractedRate\n";
    std::cout << bg->getBackgroundRun(18137) << std::endl;
    std::vector <double> evecbg = bg->ReturnBGSubtRate(0);
    std::vector <double> wvecbg = bg->ReturnBGSubtRate(1);

    for (unsigned int i=0; i<evecbg.size(); i++) {
      std::cout << evecbg[i] << " " << wvecbg[i] << std::endl;
      }
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
