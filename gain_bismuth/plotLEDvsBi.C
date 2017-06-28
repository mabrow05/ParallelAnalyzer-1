//Plots the gain as a function of run number (time)

//Should make an array of all the calibration periods and the
// runs they apply to, and then put vertical lines at each of these
#include <vector>
#include <utility>

std::vector <Int_t> getPMTEreconQuality(Int_t runNumber) {
  //Read in PMT quality file
  //cout << "Reading in PMT Quality file ...\n";
  vector <Int_t>  pmtQuality (8,0);
  Char_t temp[200];
  sprintf(temp,"%s/residuals/PMT_EreconQuality_master.dat",getenv("ANALYSIS_CODE")); 
  ifstream pmt;
  //std::cout << temp << std::endl;
  pmt.open(temp);
  Int_t run_hold;
  while (pmt >> run_hold >> pmtQuality[0] >> pmtQuality[1] >> pmtQuality[2]
	 >> pmtQuality[3] >> pmtQuality[4] >> pmtQuality[5]
	 >> pmtQuality[6] >> pmtQuality[7]) {
    if (run_hold==runNumber) break;
    if (pmt.fail()) break;
  }
  pmt.close();
  if (run_hold!=runNumber) {
    cout << "Run not found in PMT quality file!" << endl;
    exit(0);
  }
  return pmtQuality;
}

std::vector <Double_t> getLEDGainAdjuster(Int_t rn) {

  std::ifstream infile("LEDgainAdjusters.txt");

  std::vector <Double_t> adj(2,1.);
  Int_t run_hold;
  while ( infile >> run_hold >> adj[0] >> adj[1] ) {
    if ( run_hold==rn ) break;
  }
  return adj;

};

void plotLEDvsBi(TString year) {

  gStyle->SetTitleSize(0.08,"t");
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetTitleYSize(0.06);
  gStyle->SetTitleYOffset(0.6);
  gStyle->SetTitleXSize(0.06);
  gStyle->SetTitleXOffset(0.9);
  gStyle->SetLabelSize(0.06,"xyz");
  
  std::vector <Int_t> calPeriodRunEnd; // TO BE IMPLEMENTED LATER
  std::vector <Int_t> calPeriodRefRun;



  Int_t runMin, runMax;

  if (year==TString("2011-2012") ) {
    runMin = 17000;
    runMax = 20000;
    
    calPeriodRefRun.push_back(17238); calPeriodRunEnd.push_back(17298);
    calPeriodRefRun.push_back(17370); calPeriodRunEnd.push_back(17440);
    calPeriodRefRun.push_back(17521); calPeriodRunEnd.push_back(17735);
    calPeriodRefRun.push_back(17892); calPeriodRunEnd.push_back(17956);
    calPeriodRefRun.push_back(18361); calPeriodRunEnd.push_back(18387);
    calPeriodRefRun.push_back(18621); calPeriodRunEnd.push_back(18684);
    calPeriodRefRun.push_back(18749); calPeriodRunEnd.push_back(18995);
    calPeriodRefRun.push_back(19232); calPeriodRunEnd.push_back(19240);
    calPeriodRefRun.push_back(19359); calPeriodRunEnd.push_back(19545);
    calPeriodRefRun.push_back(19857); calPeriodRunEnd.push_back(20000);
  }
  else {
    runMin = 21000;
    runMax = 23300;

    calPeriodRefRun.push_back(20519); calPeriodRunEnd.push_back(20742);
    calPeriodRefRun.push_back(20820); calPeriodRunEnd.push_back(20838);
    calPeriodRefRun.push_back(21091); calPeriodRunEnd.push_back(21238);
    calPeriodRefRun.push_back(21315); calPeriodRunEnd.push_back(21606);
    calPeriodRefRun.push_back(21683); calPeriodRunEnd.push_back(21864);
    calPeriodRefRun.push_back(21918); calPeriodRunEnd.push_back(22119);
    calPeriodRefRun.push_back(22219); calPeriodRunEnd.push_back(22239);
    calPeriodRefRun.push_back(22441); calPeriodRunEnd.push_back(22631);
    calPeriodRefRun.push_back(22771); calPeriodRunEnd.push_back(23174);
  }

  std::vector < std::vector <Double_t> > runNumber( 8, std::vector<Double_t>(4000) );
  std::vector < std::vector <Double_t> > eastGain( 4, std::vector<Double_t>(4000) );
  std::vector < std::vector <Double_t> > westGain( 4, std::vector<Double_t>(4000) );
  std::vector < std::vector <Double_t> > eastGainFactor( 4, std::vector<Double_t>(4000) );
  std::vector < std::vector <Double_t> > westGainFactor( 4, std::vector<Double_t>(4000) );
  std::vector < Int_t > numRuns(8,0);

  //These are hold the reference runs specifically to be drawn in a different color
  std::vector < std::vector <Double_t> > runNumber_ref( 8, std::vector<Double_t>(20) );
  std::vector < std::vector <Double_t> > eastGain_ref( 4, std::vector<Double_t>(20) );
  std::vector < std::vector <Double_t> > westGain_ref( 4, std::vector<Double_t>(20) );
  std::vector < std::vector <Double_t> > eastGainFactor_ref( 4, std::vector<Double_t>(20) );
  std::vector < std::vector <Double_t> > westGainFactor_ref( 4, std::vector<Double_t>(20) );
  std::vector < Int_t > numRuns_ref( 8,0 );
  

  double e1=0., e2=0., e3=0., e4=0., w1=0., w2=0., w3=0., w4=0.;
  double e1_norm=0., e2_norm=0., e3_norm=0., e4_norm=0., w1_norm=0., w2_norm=0., w3_norm=0., w4_norm=0.;
  
  ifstream gainFile;

  for (Int_t rn=runMin; rn<runMax; rn++) {
    
    if (rn%500==0) std::cout << "NOW AT RUN " << rn << std::endl;

    TString filename = TString::Format("%s/gain_LED_%i.dat",getenv("GAIN_LED"),rn);
    gainFile.open(filename.Data());
    
    if (gainFile.is_open()) {
      if ( rn==17588 || (rn>19583 && rn<19616) || (rn>=16000 && rn<=17078) ) { gainFile.close(); continue; } //Groups of runs to ignore due to bad crap
      // 17588 is empty, the ranges are unused Xe runs which weren't processed, 
      // (rn>20515 && rn<21087) is all the garbage which isn't really usable from the beginning of 2012-2013

      std::vector <Int_t> pmtQuality = getPMTEreconQuality(rn);
      //std::vector <Double_t> adj = getLEDGainAdjuster(rn);
      std::vector <Double_t> adj(2,1.);

      //std::cout << "Opened " << filename << std::endl;
      gainFile >> e1 >> e1_norm
	       >> e2 >> e2_norm
	       >> e3 >> e3_norm
	       >> e4 >> e4_norm
	       >> w1 >> w1_norm
	       >> w2 >> w2_norm
	       >> w3 >> w3_norm
	       >> w4 >> w4_norm;
      
      if (pmtQuality[0]) { 
	runNumber[0][numRuns[0]]=rn; eastGain[0][numRuns[0]]=e1; eastGainFactor[0][numRuns[0]]=e1_norm*adj[0]; numRuns[0]++; 
	for ( std::vector<Int_t>::iterator it = calPeriodRefRun.begin(); it!=calPeriodRefRun.end(); ++it ) {
	  if ( (int)rn==*it ) {
	    runNumber_ref[0][numRuns_ref[0]]=rn; eastGain_ref[0][numRuns_ref[0]]=e1; eastGainFactor_ref[0][numRuns_ref[0]]=e1_norm; numRuns_ref[0]++; 
	    continue;
	  }
	}
      }
      else e1=1000.; //so it doesn't trigger the check below
      if (pmtQuality[1]) { 
	runNumber[1][numRuns[1]]=rn; eastGain[1][numRuns[1]]=e2; eastGainFactor[1][numRuns[1]]=e2_norm*adj[0]; numRuns[1]++; 
	for ( std::vector<Int_t>::iterator it = calPeriodRefRun.begin(); it!=calPeriodRefRun.end(); ++it ) {
	  if ( (int)rn==*it ) {
	    runNumber_ref[1][numRuns_ref[1]]=rn; eastGain_ref[1][numRuns_ref[1]]=e2; eastGainFactor_ref[1][numRuns_ref[1]]=e2_norm; numRuns_ref[1]++; 
	    continue;
	  }
	}
      }
      else e2=1000.; //so it doesn't trigger the check below
      if (pmtQuality[2]) {
	runNumber[2][numRuns[2]]=rn; eastGain[2][numRuns[2]]=e3; eastGainFactor[2][numRuns[2]]=e3_norm*adj[0]; numRuns[2]++; 
	for ( std::vector<Int_t>::iterator it = calPeriodRefRun.begin(); it!=calPeriodRefRun.end(); ++it ) {
	  if ( (int)rn==*it ) {
	    runNumber_ref[2][numRuns_ref[2]]=rn; eastGain_ref[2][numRuns_ref[2]]=e3; eastGainFactor_ref[2][numRuns_ref[2]]=e3_norm; numRuns_ref[2]++; 
	    continue;
	  }
	}
      }
      else e3=1000.; //so it doesn't trigger the check below
      if (pmtQuality[3]) { 
	runNumber[3][numRuns[3]]=rn; eastGain[3][numRuns[3]]=e4; eastGainFactor[3][numRuns[3]]=e4_norm*adj[0]; numRuns[3]++; 
	for ( std::vector<Int_t>::iterator it = calPeriodRefRun.begin(); it!=calPeriodRefRun.end(); ++it ) {
	  if ( (int)rn==*it ) {
	    runNumber_ref[3][numRuns_ref[3]]=rn; eastGain_ref[3][numRuns_ref[3]]=e4; eastGainFactor_ref[3][numRuns_ref[3]]=e4_norm; numRuns_ref[3]++;
	    continue;; 
	  }
	}
      }
      else e4=1000.; //so it doesn't trigger the check below
      if (pmtQuality[4]) { 
	runNumber[4][numRuns[4]]=rn; westGain[0][numRuns[4]]=w1; westGainFactor[0][numRuns[4]]=w1_norm*adj[1]; numRuns[4]++; 
	for ( std::vector<Int_t>::iterator it = calPeriodRefRun.begin(); it!=calPeriodRefRun.end(); ++it ) {
	  if ( (int)rn==*it ) {
	    runNumber_ref[4][numRuns_ref[4]]=rn; westGain_ref[0][numRuns_ref[4]]=w1; westGainFactor_ref[0][numRuns_ref[4]]=w1_norm; numRuns_ref[4]++; 
	    continue;
	  }
	}
      }
      else w1=1000.; //so it doesn't trigger the check below
      if (pmtQuality[5]) { 
	runNumber[5][numRuns[5]]=rn; westGain[1][numRuns[5]]=w2; westGainFactor[1][numRuns[5]]=w2_norm*adj[1]; numRuns[5]++; 
	for ( std::vector<Int_t>::iterator it = calPeriodRefRun.begin(); it!=calPeriodRefRun.end(); ++it ) {
	  if ( (int)rn==*it ) {
	    runNumber_ref[5][numRuns_ref[5]]=rn; westGain_ref[1][numRuns_ref[5]]=w2; westGainFactor_ref[1][numRuns_ref[5]]=w2_norm; numRuns_ref[5]++; 
	    continue;
	  }
	}
      }
      else w2=1000.; //so it doesn't trigger the check below
      if (pmtQuality[6]) { 
	runNumber[6][numRuns[6]]=rn; westGain[2][numRuns[6]]=w3; westGainFactor[2][numRuns[6]]=w3_norm*adj[1]; numRuns[6]++; 
	for ( std::vector<Int_t>::iterator it = calPeriodRefRun.begin(); it!=calPeriodRefRun.end(); ++it ) {
	  if ( (int)rn==*it ) {
	    runNumber_ref[6][numRuns_ref[6]]=rn; westGain_ref[2][numRuns_ref[6]]=w3; westGainFactor_ref[2][numRuns_ref[6]]=w3_norm; numRuns_ref[6]++;
	    continue;
	  }
	}
      }
      else w3=1000.; //so it doesn't trigger the check below
      if (pmtQuality[7]) { 
	runNumber[7][numRuns[7]]=rn; westGain[3][numRuns[7]]=w4; westGainFactor[3][numRuns[7]]=w4_norm*adj[1]; numRuns[7]++; 
	for ( std::vector<Int_t>::iterator it = calPeriodRefRun.begin(); it!=calPeriodRefRun.end(); ++it ) {
	  if ( (int)rn==*it ) {
	    runNumber_ref[7][numRuns_ref[7]]=rn; westGain_ref[3][numRuns_ref[7]]=w4; westGainFactor_ref[3][numRuns_ref[7]]=w4_norm; numRuns_ref[7]++; 
	    continue;
	  }
	}
      }
      else w4=1000.; //so it doesn't trigger the check below
      
      
      gainFile.close();
    }
  }
  

  std::vector < Int_t > numRunsBi(8,0);
  std::vector < Int_t > numRuns_refBi(8,0);

  std::vector < std::pair<Int_t,Double_t> > gainAdjusterE;
  std::vector < std::pair<Int_t,Double_t> > gainAdjusterW;

  for (Int_t rn=runMin; rn<runMax; rn++) {
    
    if (rn%500==0) std::cout << "NOW AT RUN " << rn << std::endl;

    TString filename = TString::Format("%s/gain_bismuth_%i.dat",getenv("GAIN_BISMUTH"),rn);
    gainFile.open(filename.Data());
    
    if (gainFile.is_open()) {

      std::vector <Double_t> adjE;
      std::vector <Double_t> adjW;
      
      if ( rn==17588 || (rn>19583 && rn<19616) || (rn>=16000 && rn<=17078) ) { gainFile.close(); continue; } //Groups of runs to ignore due to bad crap
      // 17588 is empty, the ranges are unused Xe runs which weren't processed, 
      // (rn>20515 && rn<21087) is all the garbage which isn't really usable from the beginning of 2012-2013

      std::vector <Int_t> pmtQuality = getPMTEreconQuality(rn);
      
      //std::cout << "Opened " << filename << std::endl;
      gainFile >> e1 >> e1_norm
	       >> e2 >> e2_norm
	       >> e3 >> e3_norm
	       >> e4 >> e4_norm
	       >> w1 >> w1_norm
	       >> w2 >> w2_norm
	       >> w3 >> w3_norm
	       >> w4 >> w4_norm;
      
      if (pmtQuality[0]) { 
        eastGain[0][numRunsBi[0]]/=e1; eastGainFactor[0][numRunsBi[0]]/=e1_norm; 
	adjE.push_back(1./eastGainFactor[0][numRunsBi[0]]);
	numRunsBi[0]++;	
	for ( std::vector<Int_t>::iterator it = calPeriodRefRun.begin(); it!=calPeriodRefRun.end(); ++it ) {
	  if ( (int)rn==*it ) {
	    eastGain_ref[0][numRuns_refBi[0]]/=e1; eastGainFactor_ref[0][numRuns_refBi[0]]/=e1_norm; numRuns_refBi[0]++;
	    continue;
	  }
	}
      }
      else e1=1000.; //so it doesn't trigger the check below
      if (pmtQuality[1]) { 
	eastGain[1][numRunsBi[1]]/=e2; eastGainFactor[1][numRunsBi[1]]/=e2_norm; 
	adjE.push_back(1./eastGainFactor[1][numRunsBi[1]]);
	numRunsBi[1]++;
	for ( std::vector<Int_t>::iterator it = calPeriodRefRun.begin(); it!=calPeriodRefRun.end(); ++it ) {
	  if ( (int)rn==*it ) {
	    eastGain_ref[1][numRuns_refBi[1]]/=e2; eastGainFactor_ref[1][numRuns_refBi[1]]/=e2_norm; numRuns_refBi[1]++;
	    continue;
	  }
	}
      }
      else e2=1000.; //so it doesn't trigger the check below
      if (pmtQuality[2]) {
	eastGain[2][numRunsBi[2]]/=e3; eastGainFactor[2][numRunsBi[2]]/=e3_norm; 
	adjE.push_back(1./eastGainFactor[2][numRunsBi[2]]);
	numRunsBi[2]++;
	for ( std::vector<Int_t>::iterator it = calPeriodRefRun.begin(); it!=calPeriodRefRun.end(); ++it ) {
	  if ( (int)rn==*it ) {
	    eastGain_ref[2][numRuns_refBi[2]]/=e3; eastGainFactor_ref[2][numRuns_refBi[2]]/=e3_norm; numRuns_refBi[2]++;
	    continue;
	  }
	}
      }
      else e3=1000.; //so it doesn't trigger the check below
      if (pmtQuality[3]) { 
	eastGain[3][numRunsBi[3]]/=e4; eastGainFactor[3][numRunsBi[3]]/=e4_norm; 
	numRunsBi[3]++;
	for ( std::vector<Int_t>::iterator it = calPeriodRefRun.begin(); it!=calPeriodRefRun.end(); ++it ) {
	  if ( (int)rn==*it ) {
	    eastGain_ref[3][numRuns_refBi[3]]/=e4; eastGainFactor_ref[3][numRuns_refBi[3]]/=e4_norm; numRuns_refBi[3]++;
	    continue;
	  }
	}
      }
      else e4=1000.; //so it doesn't trigger the check below
      if (pmtQuality[4]) { 
	westGain[0][numRunsBi[4]]/=w1; westGainFactor[0][numRunsBi[4]]/=w1_norm; 
	adjW.push_back(1./westGainFactor[0][numRunsBi[4]]);
	numRunsBi[4]++;
	for ( std::vector<Int_t>::iterator it = calPeriodRefRun.begin(); it!=calPeriodRefRun.end(); ++it ) {
	  if ( (int)rn==*it ) {
	    westGain_ref[0][numRuns_refBi[4]]/=w1; westGainFactor_ref[0][numRuns_refBi[4]]/=w1_norm; numRuns_refBi[4]++;
	    continue;
	  }
	}
      }
      else w1=1000.; //so it doesn't trigger the check below
      if (pmtQuality[5]) { 
	westGain[1][numRunsBi[5]]/=w2; westGainFactor[1][numRunsBi[5]]/=w2_norm; 
	adjW.push_back(1./westGainFactor[1][numRunsBi[5]]);
	numRunsBi[5]++;
	for ( std::vector<Int_t>::iterator it = calPeriodRefRun.begin(); it!=calPeriodRefRun.end(); ++it ) {
	  if ( (int)rn==*it ) {
	    westGain_ref[1][numRuns_refBi[5]]/=w2; westGainFactor_ref[1][numRuns_refBi[5]]/=w2_norm; numRuns_refBi[5]++;
	    continue;
	  }
	}
      }
      else w2=1000.; //so it doesn't trigger the check below
      if (pmtQuality[6]) { 
	westGain[2][numRunsBi[6]]/=w3; westGainFactor[2][numRunsBi[6]]/=w3_norm; 
	adjW.push_back(1./westGainFactor[2][numRunsBi[6]]);
	numRunsBi[6]++;
	for ( std::vector<Int_t>::iterator it = calPeriodRefRun.begin(); it!=calPeriodRefRun.end(); ++it ) {
	  if ( (int)rn==*it ) {
	    westGain_ref[2][numRuns_refBi[6]]/=w3; westGainFactor_ref[2][numRuns_refBi[6]]/=w3_norm; numRuns_refBi[6]++;
	    continue;
	  }
	}
      }
      else w3=1000.; //so it doesn't trigger the check below
      if (pmtQuality[7]) { 
	westGain[3][numRunsBi[7]]/=w4; westGainFactor[3][numRunsBi[7]]/=w4_norm; 
	numRunsBi[7]++;
	for ( std::vector<Int_t>::iterator it = calPeriodRefRun.begin(); it!=calPeriodRefRun.end(); ++it ) {
	  if ( (int)rn==*it ) {
	    westGain_ref[3][numRuns_refBi[7]]/=w4; westGainFactor_ref[3][numRuns_refBi[7]]/=w4_norm; numRuns_refBi[7]++;
	    continue;
	  }
	}
      }
      else w4=1000.; //so it doesn't trigger the check below
      
      Double_t ave = 0.;
      for (UInt_t i=0;i<adjE.size();++i) ave+=adjE[i];
      ave = (adjE.size()>0?ave/(double)adjE.size():1.);
      std::pair<Int_t,Double_t> pE(rn,ave);
      gainAdjusterE.push_back(pE);

      //std::cout << rn << "\t" <<  ave << "\t";

      ave = 0.;
      for (UInt_t ii=0;ii<adjW.size();++ii) ave+=adjW[ii];
      ave = (adjW.size()>0?ave/(double)adjW.size():1.);
      std::pair<Int_t,Double_t> pW(rn,ave);
      gainAdjusterW.push_back(pW);

      //std::cout <<  ave << "\n";
      
      gainFile.close();
    }
  }

  // Write all the gain adjusters
  std::ofstream gainAdjusterFile("LEDgainAdjusters.txt");
  
  gainAdjusterFile << std::setprecision(10);

  for (UInt_t j=0; j<gainAdjusterE.size();++j) {
    gainAdjusterFile << gainAdjusterE[j].first << "\t" << gainAdjusterE[j].second << "\t" << gainAdjusterW[j].second << "\n";
  }

  gainAdjusterFile.close();



  // making all of the TLines to reperesent the calibration periods

  std::vector <TLine*> linesE1(calPeriodRunEnd.size(),0);
  std::vector <TLine*> linesE2(calPeriodRunEnd.size(),0);
  std::vector <TLine*> linesE3(calPeriodRunEnd.size(),0);
  std::vector <TLine*> linesE4(calPeriodRunEnd.size(),0);
  std::vector <TLine*> linesW1(calPeriodRunEnd.size(),0);
  std::vector <TLine*> linesW2(calPeriodRunEnd.size(),0);
  std::vector <TLine*> linesW3(calPeriodRunEnd.size(),0);
  std::vector <TLine*> linesW4(calPeriodRunEnd.size(),0);

  std::vector <TLine*> linesE1_fac(calPeriodRunEnd.size(),0);
  std::vector <TLine*> linesE2_fac(calPeriodRunEnd.size(),0);
  std::vector <TLine*> linesE3_fac(calPeriodRunEnd.size(),0);
  std::vector <TLine*> linesE4_fac(calPeriodRunEnd.size(),0);
  std::vector <TLine*> linesW1_fac(calPeriodRunEnd.size(),0);
  std::vector <TLine*> linesW2_fac(calPeriodRunEnd.size(),0);
  std::vector <TLine*> linesW3_fac(calPeriodRunEnd.size(),0);
  std::vector <TLine*> linesW4_fac(calPeriodRunEnd.size(),0);

  Double_t min = 0., max = 0.;
  

  
    ///////////////////////////// Plotting the gain factors //////////////////////////////////


    //EAST SIDE
  TCanvas *cEast_Fac = new TCanvas("cEast_Fac","East PMTs",1000, 1200);
  cEast_Fac->Divide(1,4);

  cEast_Fac->cd(1);

  TGraph *east_fac1 = new TGraph(numRuns[0],&runNumber[0][0],&eastGainFactor[0][0]);
  east_fac1->SetTitle("PMT East 1");
  east_fac1->GetXaxis()->SetTitle("Run Number");
  east_fac1->GetYaxis()->SetTitle("LEDGain / BiGain");
  east_fac1->GetXaxis()->SetLimits(runMin, runMax);
  east_fac1->SetMarkerStyle(8);
  // east_fac1->SetMinimum(0.);
  //east_fac1->SetMaximum(4096.);
  east_fac1->Draw("AP");

  cEast_Fac->Update();
  
  min = gPad->GetUymin();
  max = gPad->GetUymax();

  for (UInt_t i=0; i<calPeriodRunEnd.size(); i++) {
    
    linesE1_fac[i] = new TLine((double)calPeriodRunEnd[i],min,(double)calPeriodRunEnd[i],max);
    linesE1_fac[i]->SetLineStyle(2);
    linesE1_fac[i]->SetLineColor(kRed);
    linesE1_fac[i]->SetLineWidth(2);
    if ( (double)calPeriodRunEnd[i] > runNumber[0][0] && (double)calPeriodRunEnd[i] < runNumber[0][numRuns[0]-1] ) {
      linesE1_fac[i]->Draw();
    }
  } 

  TGraph *east_fac1_ref = new TGraph(numRuns_ref[0],&runNumber_ref[0][0],&eastGainFactor_ref[0][0]);
  east_fac1_ref->SetTitle("PMT East 1");
  east_fac1_ref->GetXaxis()->SetTitle("Run Number");
  east_fac1_ref->GetYaxis()->SetTitle("LEDGain / BiGain");
  east_fac1_ref->SetMarkerStyle(8);
  east_fac1_ref->SetMarkerColor(kRed);
  // east_fac1_ref->SetMinimum(0.);
  //east_fac1_ref->SetMaximum(4096.);
  east_fac1_ref->Draw("P SAME");



  cEast_Fac->cd(2);

  TGraph *east_fac2 = new TGraph(numRuns[1],&runNumber[1][0],&eastGainFactor[1][0]);
  east_fac2->SetTitle("PMT East 2");
  east_fac2->GetXaxis()->SetTitle("Run Number");
  east_fac2->GetYaxis()->SetTitle("LEDGain / BiGain");
  east_fac2->GetXaxis()->SetLimits(runMin, runMax);
  east_fac2->SetMarkerStyle(8);
  //east_fac2->SetMinimum(0.);
  //east_fac2->SetMaximum(4096.);
  east_fac2->Draw("AP");

  cEast_Fac->Update();
  
  min = gPad->GetUymin();
  max = gPad->GetUymax();

  for (UInt_t i=0; i<calPeriodRunEnd.size(); i++) {
    
    linesE2_fac[i] = new TLine((double)calPeriodRunEnd[i],min,(double)calPeriodRunEnd[i],max);
    linesE2_fac[i]->SetLineStyle(2);
    linesE2_fac[i]->SetLineColor(kRed);
    linesE2_fac[i]->SetLineWidth(2);
    if ( (double)calPeriodRunEnd[i] > runNumber[1][0] && (double)calPeriodRunEnd[i] < runNumber[1][numRuns[1]-1] ) {
      linesE2_fac[i]->Draw();
    }
  } 

  TGraph *east_fac2_ref = new TGraph(numRuns_ref[1],&runNumber_ref[1][0],&eastGainFactor_ref[1][0]);
  east_fac2_ref->SetTitle("PMT East 2");
  east_fac2_ref->GetXaxis()->SetTitle("Run Number");
  east_fac2_ref->GetYaxis()->SetTitle("LEDGain / BiGain");
  east_fac2_ref->SetMarkerStyle(8);
  east_fac2_ref->SetMarkerColor(kRed);
  // east_fac2_ref->SetMinimum(0.);
  //east_fac2_ref->SetMaximum(4096.);
  east_fac2_ref->Draw("P SAME");

  

  cEast_Fac->cd(3);

  TGraph *east_fac3 = new TGraph(numRuns[2],&runNumber[2][0],&eastGainFactor[2][0]);
  east_fac3->SetTitle("PMT East 3");
  east_fac3->GetXaxis()->SetTitle("Run Number");
  east_fac3->GetYaxis()->SetTitle("LEDGain / BiGain");
  east_fac3->GetXaxis()->SetLimits(runMin, runMax);
  east_fac3->SetMarkerStyle(8);
  //east_fac3->SetMinimum(0.);
  //east_fac3->SetMaximum(4096.);
  east_fac3->Draw("AP");

  cEast_Fac->Update();
  
  min = gPad->GetUymin();
  max = gPad->GetUymax();

  for (UInt_t i=0; i<calPeriodRunEnd.size(); i++) {
    
    linesE3_fac[i] = new TLine((double)calPeriodRunEnd[i],min,(double)calPeriodRunEnd[i],max);
    linesE3_fac[i]->SetLineStyle(2);
    linesE3_fac[i]->SetLineColor(kRed);
    linesE3_fac[i]->SetLineWidth(2);
    if ( (double)calPeriodRunEnd[i] > runNumber[2][0] && (double)calPeriodRunEnd[i] < runNumber[2][numRuns[2]-1] ) {
      linesE3_fac[i]->Draw();
    }
  } 

  TGraph *east_fac3_ref = new TGraph(numRuns_ref[2],&runNumber_ref[2][0],&eastGainFactor_ref[2][0]);
  east_fac3_ref->SetTitle("PMT East 3");
  east_fac3_ref->GetXaxis()->SetTitle("Run Number");
  east_fac3_ref->GetYaxis()->SetTitle("LEDGain / BiGain");
  east_fac3_ref->SetMarkerStyle(8);
  east_fac3_ref->SetMarkerColor(kRed);
  // east_fac3_ref->SetMinimum(0.);
  //east_fac3_ref->SetMaximum(4096.);
  east_fac3_ref->Draw("P SAME");


  cEast_Fac->cd(4);

  TGraph *east_fac4 = new TGraph(numRuns[3],&runNumber[3][0],&eastGainFactor[3][0]);
  east_fac4->SetTitle("PMT East 4");
  east_fac4->GetXaxis()->SetTitle("Run Number");
  east_fac4->GetYaxis()->SetTitle("LEDGain / BiGain");
  east_fac4->GetXaxis()->SetLimits(runMin, runMax);
  east_fac4->SetMarkerStyle(8);
  //east_fac4->SetMinimum(0.);
  //east_fac4->SetMaximum(4096.);
  east_fac4->Draw("AP");

  cEast_Fac->Update();
  
  min = gPad->GetUymin();
  max = gPad->GetUymax();

  for (UInt_t i=0; i<calPeriodRunEnd.size(); i++) {
    
    linesE4_fac[i] = new TLine((double)calPeriodRunEnd[i],min,(double)calPeriodRunEnd[i],max);
    linesE4_fac[i]->SetLineStyle(2);
    linesE4_fac[i]->SetLineColor(kRed);
    linesE4_fac[i]->SetLineWidth(2);
    if ( (double)calPeriodRunEnd[i] > runNumber[3][0] && (double)calPeriodRunEnd[i] < runNumber[3][numRuns[3]-1] ) {
      linesE4_fac[i]->Draw();
    }
  } 

  TGraph *east_fac4_ref = new TGraph(numRuns_ref[3],&runNumber_ref[3][0],&eastGainFactor_ref[3][0]);
  east_fac4_ref->SetTitle("PMT East 4");
  east_fac4_ref->GetXaxis()->SetTitle("Run Number");
  east_fac4_ref->GetYaxis()->SetTitle("LEDGain / BiGain");
  east_fac4_ref->SetMarkerStyle(8);
  east_fac4_ref->SetMarkerColor(kRed);
  // east_fac4_ref->SetMinimum(0.);
  //east_fac4_ref->SetMaximum(4096.);
  east_fac4_ref->Draw("P SAME");


  //WEST SIDE
  TCanvas *cWest_Fac = new TCanvas("cWest_Fac","West PMTs",1000, 1200);
  cWest_Fac->Divide(1,4);

  cWest_Fac->cd(1);

  TGraph *west_fac1 = new TGraph(numRuns[4],&runNumber[4][0],&westGainFactor[0][0]);
  west_fac1->SetTitle("PMT West 1");
  west_fac1->GetXaxis()->SetTitle("Run Number");
  west_fac1->GetYaxis()->SetTitle("LEDGain / BiGain");
  west_fac1->GetXaxis()->SetLimits(runMin, runMax);
  west_fac1->SetMarkerStyle(8);
  //west_fac1->SetMinimum(0.);
  //west_fac1->SetMaximum(4096.);
  west_fac1->Draw("AP");

  cWest_Fac->Update();
  
  min = gPad->GetUymin();
  max = gPad->GetUymax();

  for (UInt_t i=0; i<calPeriodRunEnd.size(); i++) {
    
    linesW1_fac[i] = new TLine((double)calPeriodRunEnd[i],min,(double)calPeriodRunEnd[i],max);
    linesW1_fac[i]->SetLineStyle(2);
    linesW1_fac[i]->SetLineColor(kRed);
    linesW1_fac[i]->SetLineWidth(2);
    if ( (double)calPeriodRunEnd[i] > runNumber[4][0] && (double)calPeriodRunEnd[i] < runNumber[4][numRuns[4]-1] ) {
      linesW1_fac[i]->Draw();
    }
  } 

  TGraph *west_fac1_ref = new TGraph(numRuns_ref[4],&runNumber_ref[4][0],&westGainFactor_ref[0][0]);
  west_fac1_ref->SetTitle("PMT West 1");
  west_fac1_ref->GetXaxis()->SetTitle("Run Number");
  west_fac1_ref->GetYaxis()->SetTitle("LEDGain / BiGain");
  west_fac1_ref->SetMarkerStyle(8);
  west_fac1_ref->SetMarkerColor(kRed);
  // west_fac1_ref->SetMinimum(0.);
  //west_fac1_ref->SetMaximum(4096.);
  west_fac1_ref->Draw("P SAME");



  cWest_Fac->cd(2);

  TGraph *west_fac2 = new TGraph(numRuns[5],&runNumber[5][0],&westGainFactor[1][0]);
  west_fac2->SetTitle("PMT West 2");
  west_fac2->GetXaxis()->SetTitle("Run Number");
  west_fac2->GetYaxis()->SetTitle("LEDGain / BiGain");
  west_fac2->GetXaxis()->SetLimits(runMin, runMax);
  west_fac2->SetMarkerStyle(8);
  //west_fac2->SetMinimum(0.);
  //west_fac2->SetMaximum(4096.);
  west_fac2->Draw("AP");

  cWest_Fac->Update();
  
  min = gPad->GetUymin();
  max = gPad->GetUymax();

  for (UInt_t i=0; i<calPeriodRunEnd.size(); i++) {
    
    linesW2_fac[i] = new TLine((double)calPeriodRunEnd[i],min,(double)calPeriodRunEnd[i],max);
    linesW2_fac[i]->SetLineStyle(2);
    linesW2_fac[i]->SetLineColor(kRed);
    linesW2_fac[i]->SetLineWidth(2);
    if ( (double)calPeriodRunEnd[i] > runNumber[5][0] && (double)calPeriodRunEnd[i] < runNumber[5][numRuns[5]-1] ) {
      linesW2_fac[i]->Draw();
    }
  } 

  TGraph *west_fac2_ref = new TGraph(numRuns_ref[5],&runNumber_ref[5][0],&westGainFactor_ref[1][0]);
  west_fac2_ref->SetTitle("PMT West 2");
  west_fac2_ref->GetXaxis()->SetTitle("Run Number");
  west_fac2_ref->GetYaxis()->SetTitle("LEDGain / BiGain");
  west_fac2_ref->SetMarkerStyle(8);
  west_fac2_ref->SetMarkerColor(kRed);
  // west_fac2_ref->SetMinimum(0.);
  //west_fac2_ref->SetMaximum(4096.);
  west_fac2_ref->Draw("P SAME");


  cWest_Fac->cd(3);

  TGraph *west_fac3 = new TGraph(numRuns[6],&runNumber[6][0],&westGainFactor[2][0]);
  west_fac3->SetTitle("PMT West 3");
  west_fac3->GetXaxis()->SetTitle("Run Number");
  west_fac3->GetYaxis()->SetTitle("LEDGain / BiGain");
  west_fac3->GetXaxis()->SetLimits(runMin, runMax);
  west_fac3->SetMarkerStyle(8);
  //west_fac3->SetMinimum(0.);
  //west_fac3->SetMaximum(4096.);
  west_fac3->Draw("AP");

  cWest_Fac->Update();
  
  min = gPad->GetUymin();
  max = gPad->GetUymax();

  for (UInt_t i=0; i<calPeriodRunEnd.size(); i++) {
    
    linesW3_fac[i] = new TLine((double)calPeriodRunEnd[i],min,(double)calPeriodRunEnd[i],max);
    linesW3_fac[i]->SetLineStyle(2);
    linesW3_fac[i]->SetLineColor(kRed);
    linesW3_fac[i]->SetLineWidth(2);
    if ( (double)calPeriodRunEnd[i] > runNumber[6][0] && (double)calPeriodRunEnd[i] < runNumber[6][numRuns[6]-1] ) {
      linesW3_fac[i]->Draw();
    }
  } 

  TGraph *west_fac3_ref = new TGraph(numRuns_ref[6],&runNumber_ref[6][0],&westGainFactor_ref[2][0]);
  west_fac3_ref->SetTitle("PMT West 3");
  west_fac3_ref->GetXaxis()->SetTitle("Run Number");
  west_fac3_ref->GetYaxis()->SetTitle("LEDGain / BiGain");
  west_fac3_ref->SetMarkerStyle(8);
  west_fac3_ref->SetMarkerColor(kRed);
  // west_fac3_ref->SetMinimum(0.);
  //west_fac3_ref->SetMaximum(4096.);
  west_fac3_ref->Draw("P SAME");


  if (numRuns[7]>0) {
    cWest_Fac->cd(4);
    
    
    TGraph *west_fac4 = new TGraph(numRuns[7],&runNumber[7][0],&westGainFactor[3][0]);
    west_fac4->SetTitle("PMT West 4");
    west_fac4->GetXaxis()->SetTitle("Run Number");
    west_fac4->GetYaxis()->SetTitle("LEDGain / BiGain");
    west_fac4->GetXaxis()->SetLimits(runMin, runMax);
    west_fac4->SetMarkerStyle(8);
    //west_fac4->SetMinimum(0.);
    //west_fac4->SetMaximum(4096.);
    west_fac4->Draw("AP");

    cWest_Fac->Update();
  
    min = gPad->GetUymin();
    max = gPad->GetUymax();
    
    for (UInt_t i=0; i<calPeriodRunEnd.size(); i++) {
      
      linesW4_fac[i] = new TLine((double)calPeriodRunEnd[i],min,(double)calPeriodRunEnd[i],max);
      linesW4_fac[i]->SetLineStyle(2);
      linesW4_fac[i]->SetLineColor(kRed);
      linesW4_fac[i]->SetLineWidth(2);
      if ( (double)calPeriodRunEnd[i] > runNumber[7][0] && (double)calPeriodRunEnd[i] < runNumber[7][numRuns[7]-1] ) {
	linesW4_fac[i]->Draw();
      }
    } 

    TGraph *west_fac4_ref = new TGraph(numRuns_ref[7],&runNumber_ref[7][0],&westGainFactor_ref[3][0]);
    west_fac4_ref->SetTitle("PMT West 4");
    west_fac4_ref->GetXaxis()->SetTitle("Run Number");
    west_fac4_ref->GetYaxis()->SetTitle("LEDGain / BiGain");
    west_fac4_ref->SetMarkerStyle(8);
    west_fac4_ref->SetMarkerColor(kRed);
    // west_fac4_ref->SetMinimum(0.);
    //west_fac4_ref->SetMaximum(4096.);
    west_fac4_ref->Draw("P SAME");
    
  }
    

  TString fnamebase = TString::Format("%s_LEDvsBi.pdf",year.Data());
  //cEast->Print(TString::Format("%s(",fnamebase.Data()));
  //cWest->Print(TString::Format("%s",fnamebase.Data()));
  cEast_Fac->Print(TString::Format("%s(",fnamebase.Data()));
  cWest_Fac->Print(TString::Format("%s)",fnamebase.Data()));
  

}
