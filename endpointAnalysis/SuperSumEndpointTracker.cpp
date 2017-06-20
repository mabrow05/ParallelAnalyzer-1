#include "BetaDecayTools.hh"
#include "DataTree.hh"
#include <TString.h>
#include <TH1D.h>
#include <TGraphErrors.h>
#include <TMath.h>

#include <iostream>
#include <fstream>
#include <algorithm>

std::vector <Int_t> badOct {7,9,59,60,61,62,63,64,65,66,67,91,93,101,107,121};

std::vector <int>  readOctetFile(int octet); 

int main(int argc, char *argv[]) {

  if (argc!=4) {
    std::cout << "Usage: ./endpointTracker.exe [octet start] [octet end] [bool simulation]\n";
    //std::cout << "The code will produce comparisons for every octet in the range given,\non an octet-by-octet basis, and as a whole, using the Super-Sum\n";
    exit(0);
    
  }
  
  //Show that we can do the Kurie fit...

  int octetMin = atoi(argv[1]);
  int octetMax = atoi(argv[2]);

  bool sim = ( std::string(argv[3])=="true" ) ? true : false;
  std::string prefix = sim ? "SIM" : "UK";

  int nBins = 100;
  double minRange = 0.;
  double maxRange = 1000.;

  ///// Loop over all octets in range
  
  

  for (int octetNum=octetMin; octetNum<octetMax+1; octetNum++) {
    
    if (std::find(badOct.begin(), badOct.end(),octetNum) != badOct.end()) { continue; } //Checking if octet should be ignored for data quality reasons


    Double_t totalTimeON_E, totalTimeON_W, totalTimeOFF_E, totalTimeOFF_W;
    totalTimeON_E = totalTimeON_W = totalTimeOFF_E = totalTimeOFF_W = 0.;

    Double_t totalBGTimeON_E, totalBGTimeON_W, totalBGTimeOFF_E, totalBGTimeOFF_W;
    totalBGTimeON_E = totalBGTimeON_W = totalBGTimeOFF_E = totalBGTimeOFF_W = 0.;
 
    //Count vectors for holding all of the total (BG subtracted) counts for each spin state

    //type0=0, type1=1, type2=2, type3=3 ALL=4
    std::vector <double> W_TotalCountsON(nBins,0.); 
    std::vector <double> W_TotalErrorON(nBins,0.);
    std::vector <double> W_TotalCountsOFF(nBins,0.); 
    std::vector <double> W_TotalErrorOFF(nBins,0.);

    std::vector <double> E_TotalCountsON(nBins,0.); 
    std::vector <double> E_TotalErrorON(nBins,0.);
    std::vector <double> E_TotalCountsOFF(nBins,0.); 
    std::vector <double> E_TotalErrorOFF(nBins,0.);

    std::vector <double> W_BGTotalCountsON(nBins,0.); 
    std::vector <double> W_BGTotalErrorON(nBins,0.);
    std::vector <double> W_BGTotalCountsOFF(nBins,0.); 
    std::vector <double> W_BGTotalErrorOFF(nBins,0.);

    std::vector <double> E_BGTotalCountsON(nBins,0.); 
    std::vector <double> E_BGTotalErrorON(nBins,0.);
    std::vector <double> E_BGTotalCountsOFF(nBins,0.); 
    std::vector <double> E_BGTotalErrorOFF(nBins,0.);
    
    // read in runs in octet and then open the rate files and fill the vectors according
    // to the type of run
    
    std::ifstream octFile(TString::Format("%s/All_Octets/octet_list_%i.dat",
					  getenv("OCTET_LIST"),octetNum));
    Double_t en, fg, bg, time, bg_time;
    std::string txt_hold;
    Int_t bin = 0;
    
    Int_t runNumber = 0;
    std::string runType = "";
    
    while ( octFile >> runType >> runNumber ) {
      
      std::ifstream rateFile;
      
      /////////////// Spin flipper off///////////////////
      if ( runType == "A2" || runType == "A10" || runType == "B5" || runType == "B7" ) {	
	
	//////////////// EAST RATES //////////////////
	
	
	// Type 0
	rateFile.open(TString::Format("../Asymmetry/BinByBinComparison/%s_Erun%i_anaChD.dat", prefix.c_str(),runNumber));
	rateFile >> txt_hold >> txt_hold >> txt_hold >> time;
	rateFile >> txt_hold >> txt_hold >> txt_hold >> bg_time;
	
	totalTimeOFF_E = time;
	totalBGTimeOFF_E = bg_time; 

	bin = 0; 
	
	while ( rateFile >> en >> fg >> bg && bin<nBins)  { 
	  E_TotalCountsOFF[bin] += (fg - bg)*time;
	  E_BGTotalCountsOFF[bin] += (bg)*bg_time;
	  //	  std::cout << E_TotalCountsOFF[bin] << std:: endl;
	  ++bin;
	}
	
	rateFile.close();
	
	
	
	///////////////// WEST RATES ///////////////////
	
	
	
	// Type 0 only
	rateFile.open(TString::Format("../Asymmetry/BinByBinComparison/%s_Wrun%i_anaChD.dat", prefix.c_str(), runNumber));
	
	rateFile >> txt_hold >> txt_hold >> txt_hold >> time;
	rateFile >> txt_hold >> txt_hold >> txt_hold >> bg_time;
	
	totalTimeOFF_W = time;
	totalBGTimeOFF_W = bg_time; 
	
        bin = 0;
	
	while ( rateFile >> en >> fg >> bg && bin<nBins)  { 
	  W_TotalCountsOFF[bin] += (fg - bg)*time;
	  W_BGTotalCountsOFF[bin] += (bg)*bg_time;
	  ++bin;
	}
	
	rateFile.close();
      }
	
	
	
	
	
      /////////////// Spin flipper ON ///////////////////
      if ( runType == "B2" || runType == "B10" || runType == "A5" || runType == "A7" ) {	
	
	//////////////// EAST RATES //////////////////
	
	
	
	
	// Type 0 only
	rateFile.open(TString::Format("../Asymmetry/BinByBinComparison/%s_Erun%i_anaChD.dat", prefix.c_str(), runNumber));
	
	rateFile >> txt_hold >> txt_hold >> txt_hold >> time;
	rateFile >> txt_hold >> txt_hold >> txt_hold >> bg_time;
	
	totalTimeON_E = time;
	totalBGTimeON_E = bg_time; 
	
	bin = 0;
	
	while ( rateFile >> en >> fg >> bg && bin<nBins)  { 
	  E_TotalCountsON[bin] += (fg - bg)*time;
	  E_BGTotalCountsON[bin] += (bg)*bg_time;
	  ++bin;
	}
	
	rateFile.close();
	
	
	
	
	///////////////// WEST RATES ///////////////////
	
	
	// Type 0 only
	rateFile.open(TString::Format("../Asymmetry/BinByBinComparison/%s_Wrun%i_anaChD.dat", prefix.c_str(), runNumber));
	
	rateFile >> txt_hold >> txt_hold >> txt_hold >> time;
	rateFile >> txt_hold >> txt_hold >> txt_hold >> bg_time;

	totalTimeON_W = time;
	totalBGTimeON_W = bg_time; 
	
	bin = 0;
	
	while ( rateFile >> en >> fg >> bg && bin<nBins)  { 
	  W_TotalCountsON[bin] += (fg - bg)*time;
	  W_BGTotalCountsON[bin] += (bg)*bg_time;
	  ++bin;
	}
	
	rateFile.close();
	
	
      }

    }
    

    // Calculate the simple errors
    for ( unsigned int b = 0; b<100; ++b ) {
      
      W_TotalErrorON[b] = sqrt( fabs(W_TotalCountsON[b]) );
      W_TotalErrorOFF[b] = sqrt( fabs(W_TotalCountsOFF[b]) );
      E_TotalErrorON[b] = sqrt( fabs(E_TotalCountsON[b]) );
      E_TotalErrorOFF[b] = sqrt( fabs(E_TotalCountsOFF[b]) );
      
      W_BGTotalErrorON[b] = sqrt( fabs(W_BGTotalCountsON[b]) );
      W_BGTotalErrorOFF[b] = sqrt( fabs(W_BGTotalCountsOFF[b]) );
      E_BGTotalErrorON[b] = sqrt( fabs(E_BGTotalCountsON[b]) );
      E_BGTotalErrorOFF[b] = sqrt( fabs(E_BGTotalCountsOFF[b]) );
      
    }
  
  
    // Calculate the super-sum in every bin
    std::vector <double> data_SuperSum(nBins,0.); 
    std::vector <double> data_SuperSumError(nBins,0.); 
    
    std::vector <double> BGdata_SuperSum(nBins,0.); 
    std::vector <double> BGdata_SuperSumError(nBins,0.); 
    
    
    // First for data
    for ( unsigned int b = 0; b<nBins; ++b ) {
      
      double sfON[2]={0.};
      double sfOFF[2]={0.};
      double sfON_err[2]={0.};
      double sfOFF_err[2]={0.};
      
      
      // AFP Off
      sfOFF[0] = E_TotalCountsOFF[b] / totalTimeOFF_E;
      sfOFF_err[0] = E_TotalErrorOFF[b] / totalTimeOFF_E;
      
      sfOFF[1] = W_TotalCountsOFF[b] / totalTimeOFF_W;
      sfOFF_err[1] = W_TotalErrorOFF[b] / totalTimeOFF_W;
      
      // AFP ON
      sfON[0] = E_TotalCountsON[b] / totalTimeON_E;
      sfON_err[0] = E_TotalErrorON[b] / totalTimeON_E;
      
      sfON[1] = W_TotalCountsON[b] / totalTimeON_W;
      sfON_err[1] = W_TotalErrorON[b] / totalTimeON_W;
      
      
      
      //The geometric mean as defined by http://www.arpapress.com/Volumes/Vol11Issue3/IJRRAS_11_3_08.pdf
      double R1 = (sfON[1]>0. && sfOFF[0]>0.) ? sqrt(sfOFF[0]*sfON[1]) : (sfON[1]<0. && sfOFF[0]<0.) ? -sqrt(sfOFF[0]*sfON[1]) : 0.5*(sfOFF[0]+sfON[1]); 
      double R2 = (sfON[0]>0. && sfOFF[1]>0.) ? sqrt(sfOFF[1]*sfON[0]) : (sfON[0]<0. && sfOFF[1]<0.) ? -sqrt(sfOFF[1]*sfON[0]) : 0.5*(sfOFF[1]+sfON[0]);
      double deltaR1 = ( ( (sfON[1]>0. && sfOFF[0]>0.) || (sfON[1]<0. && sfOFF[0]<0.) ) ? 0.5*sqrt((TMath::Power(sfOFF[0]*sfON_err[1],2)+TMath::Power(sfON[1]*sfOFF_err[0],2))/(sfON[1]*sfOFF[0])) : 
			 0.5*sqrt( TMath::Power(sfON_err[1],2) + TMath::Power(sfOFF_err[0],2) ) );
      double deltaR2 = ((sfON[0]>0. && sfOFF[1]>0.) || (sfON[0]<0. && sfOFF[1]<0.)) ? 0.5*sqrt((TMath::Power(sfOFF[1]*sfON_err[0],2)+TMath::Power(sfON[0]*sfOFF_err[1],2))/(sfON[0]*sfOFF[1])) : 
	0.5*sqrt(TMath::Power(sfON_err[0],2) + TMath::Power(sfOFF_err[1],2));
      
      data_SuperSum[b] = R1 + R2;
      //std::cout << sfOFF[0] << " " << sfOFF_err[0] << std::endl;
      data_SuperSumError[b] = sqrt(TMath::Power(deltaR1,2) + TMath::Power(deltaR2,2));
      
      std::cout << data_SuperSum[b] << "\t" << data_SuperSumError[b] << "\n";
      
    }
  

    // Now for BG
    
    if ( !sim ) {
      for ( unsigned int b = 0; b<nBins; ++b ) {
	
	double sfON[2]={0.};
	double sfOFF[2]={0.};
	double sfON_err[2]={0.};
	double sfOFF_err[2]={0.};
	
	
	// AFP Off
	sfOFF[0] = E_BGTotalCountsOFF[b] / totalBGTimeOFF_E;
	sfOFF_err[0] = E_BGTotalErrorOFF[b] / totalBGTimeOFF_E;
	
	sfOFF[1] = W_BGTotalCountsOFF[b] / totalBGTimeOFF_W;
	sfOFF_err[1] = W_BGTotalErrorOFF[b] / totalBGTimeOFF_W;
	
	// AFP ON
	sfON[0] = E_BGTotalCountsON[b] / totalBGTimeON_E;
	sfON_err[0] = E_BGTotalErrorON[b] / totalBGTimeON_E;
	
	sfON[1] = W_BGTotalCountsON[b] / totalBGTimeON_W;
	sfON_err[1] = W_BGTotalErrorON[b] / totalBGTimeON_W;
	
	
	
	//The geometric mean as defined by http://www.arpapress.com/Volumes/Vol11Issue3/IJRRAS_11_3_08.pdf
	double R1 = (sfON[1]>0. && sfOFF[0]>0.) ? sqrt(sfOFF[0]*sfON[1]) : (sfON[1]<0. && sfOFF[0]<0.) ? -sqrt(sfOFF[0]*sfON[1]) : 0.5*(sfOFF[0]+sfON[1]); 
	double R2 = (sfON[0]>0. && sfOFF[1]>0.) ? sqrt(sfOFF[1]*sfON[0]) : (sfON[0]<0. && sfOFF[1]<0.) ? -sqrt(sfOFF[1]*sfON[0]) : 0.5*(sfOFF[1]+sfON[0]);
	double deltaR1 = ( ( (sfON[1]>0. && sfOFF[0]>0.) || (sfON[1]<0. && sfOFF[0]<0.) ) ? 0.5*sqrt((TMath::Power(sfOFF[0]*sfON_err[1],2)+TMath::Power(sfON[1]*sfOFF_err[0],2))/(sfON[1]*sfOFF[0])) : 
			   0.5*sqrt( TMath::Power(sfON_err[1],2) + TMath::Power(sfOFF_err[0],2) ) );
	double deltaR2 = ((sfON[0]>0. && sfOFF[1]>0.) || (sfON[0]<0. && sfOFF[1]<0.)) ? 0.5*sqrt((TMath::Power(sfOFF[1]*sfON_err[0],2)+TMath::Power(sfON[0]*sfOFF_err[1],2))/(sfON[0]*sfOFF[1])) : 
	  0.5*sqrt(TMath::Power(sfON_err[0],2) + TMath::Power(sfOFF_err[1],2));
	
	BGdata_SuperSum[b] = R1 + R2;
	//std::cout << sfOFF[0] << " " << sfOFF_err[0] << std::endl;
	BGdata_SuperSumError[b] = sqrt(TMath::Power(deltaR1,2) + TMath::Power(deltaR2,2));
	
	//std::cout << data_SuperSum[b] << "\t" << data_SuperSumError[b] << "\n";
	
      }
    }
  
    //nRun++;
    

    std::cout << "\nMoving to Kurie Fits...\n\n";
    TFile *fout = new TFile(TString::Format("%s/FinalEndpoints/%s_octet_%i_SuperSum.root",getenv("ENDPOINT_ANALYSIS"),prefix.c_str(),octetNum),"RECREATE");  
    std::cout << "Made output rootfile...\n\n";
    
    TH1D *spec = new TH1D("spec","Super-Sum Spectrum", nBins, minRange, maxRange);
    
    // loop over each bin for each run and make weighted average
    
    for ( int bin=1; bin<=nBins; ++bin ) {
      
      spec->SetBinContent(bin, data_SuperSum[bin-1]);
      spec->SetBinError(bin, data_SuperSumError[bin-1]);
      
    }
    
    
    
    TGraphErrors kurie; //These will hold the kurie plots
    
    //Now do some Kurie Fitting and write endpoints to file
    
    // open file for writing out all endpoints for each detector in a txt file
    std::ofstream epfile(TString::Format("%s/FinalEndpoints/%s_finalEndpoints_Octet_%i_SuperSum.txt",
					 getenv("ENDPOINT_ANALYSIS"),prefix.c_str(),octetNum).Data());
    
    KurieFitter kf;
    
    kf.FitSpectrum(spec, 300., 550., 1.);
    epfile << "East_Endpoint: " << kf.returnK0() << " +/- " << kf.returnK0err() << "\n";
    
    kurie = kf.returnKuriePlot();
    kurie.SetName("kurie");
    kurie.Write();
    
    fout->Write();
    delete fout;
  }



  return 0;
  
}

  
std::vector <int>  readOctetFile(int octet) {

  std::vector <int> runs;
  
  char fileName[200];
  sprintf(fileName,"%s/All_Octets/octet_list_%i.dat",getenv("OCTET_LIST"),octet);
  std::ifstream infile(fileName);
  std::string runTypeHold;
  int runNumberHold;
  int numRuns = 0;
  // Populate the map with the runType and runNumber from this octet
  while (infile >> runTypeHold >> runNumberHold) {

    if (runTypeHold=="A2" || runTypeHold=="A5" || runTypeHold=="A7" || runTypeHold=="A10" || 
	runTypeHold=="B2" || runTypeHold=="B5" || runTypeHold=="B7" || runTypeHold=="B10" )  {
     
      runs.push_back(runNumberHold);

    }
    numRuns++;
  }

  infile.close();
 
  std::cout << "Read in octet file for octet " << octet << "\n";
  return runs;

};
