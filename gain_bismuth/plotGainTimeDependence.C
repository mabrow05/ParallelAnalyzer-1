//Plots the gain as a function of run number (time)

//Should make an array of all the calibration periods and the
// runs they apply to, and then put vertical lines at each of these
#include <vector>

std::vector <Int_t> getPMTQuality(Int_t runNumber) {
  //Read in PMT quality file
  //cout << "Reading in PMT Quality file ...\n";
  vector <Int_t>  pmtQuality (8,0);
  Char_t temp[200];
  sprintf(temp,"%s/residuals/PMT_runQuality_master.dat",getenv("ANALYSIS_CODE")); 
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

void plotGainTimeDependence(TString year) {

  gStyle->SetTitleSize(0.08,"t");
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetTitleYSize(0.06);
  gStyle->SetTitleYOffset(0.4);
  gStyle->SetTitleXSize(0.06);
  gStyle->SetTitleXOffset(0.9);
  gStyle->SetLabelSize(0.06,"xyz");
  
  std::vector <Int_t> calPeriods(10,0); // TO BE IMPLEMENTED LATER

  Int_t runMin, runMax;

  if (year==TString("2011-2012") ) {
    runMin = 16000;//Skipping first Xe run period due to not using it for anything
    runMax = 20000;
  }
  else {
    runMin = 20000;
    runMax = 23500;
  }

  std::vector < std::vector <Int_t> > runNumber( 8, std::vector<Int_t>(0) );
  std::vector < std::vector <Int_t> > eastGain( 4, std::vector<Int_t>(0) );
  std::vector < std::vector <Int_t> > westGain( 4, std::vector<Int_t>(0) );

  double e1=0., e2=0., e3=0., e4=0., w1=0., w2=0., w3=0., w4=0.;
  double e1_norm=0., e2_norm=0., e3_norm=0., e4_norm=0., w1_norm=0., w2_norm=0., w3_norm=0., w4_norm=0.;
  
  ifstream gainFile;

  for (Int_t rn=runMin; rn<runMax; rn++) {
    
    if (rn%500==0) std::cout << "NOW AT RUN " << rn << std::endl;

    TString filename = TString::Format("%s/gain_bismuth_%i.dat",getenv("GAIN_BISMUTH"),rn);
    gainFile.open(filename.Data());
    
    if (gainFile.is_open()) {
      if (rn==17588 || (rn>19583 && rn<19616) || (rn>16000 && rn<17079) ) { gainFile.close(); continue; } //Groups of runs to ignore due to bad crap
      // 17588 is empty, the ranges are unused Xe runs which weren't processed

      std::vector <Int_t> pmtQuality = getPMTQuality(rn);
      
      //std::cout << "Opened " << filename << std::endl;
      gainFile >> e1 >> e1_norm
	       >> e2 >> e2_norm
	       >> e3 >> e3_norm
	       >> e4 >> e4_norm
	       >> w1 >> w1_norm
	       >> w2 >> w2_norm
	       >> w3 >> w3_norm
	       >> w4 >> w4_norm;
      
      if (pmtQuality[0]) { runNumber[0].push_back(rn); eastGain[0].push_back(e1); }
      else e1=1000.; //so it doesn't trigger the check below
      if (pmtQuality[1]) { runNumber[1].push_back(rn); eastGain[1].push_back(e2); }
      else e2=1000.; //so it doesn't trigger the check below
      if (pmtQuality[2]) { runNumber[2].push_back(rn); eastGain[2].push_back(e3); }
      else e3=1000.; //so it doesn't trigger the check below
      if (pmtQuality[3]) { runNumber[3].push_back(rn); eastGain[3].push_back(e4); }
      else e4=1000.; //so it doesn't trigger the check below
      if (pmtQuality[4]) { runNumber[4].push_back(rn); westGain[0].push_back(w1); }
      else w1=1000.; //so it doesn't trigger the check below
      if (pmtQuality[5]) { runNumber[5].push_back(rn); westGain[1].push_back(w2); }
      else w2=1000.; //so it doesn't trigger the check below
      if (pmtQuality[6]) { runNumber[6].push_back(rn); westGain[2].push_back(w3); }
      else w3=1000.; //so it doesn't trigger the check below
      if (pmtQuality[7]) { runNumber[7].push_back(rn); westGain[3].push_back(w4); }
      else w4=1000.; //so it doesn't trigger the check below
      

      if (e1<1000. || e2<1000. || e3<1000. || e4<1000. || w1<1000. || w2<1000. || w3<1000. || w4<1000.) {

	std::cout << rn << "\t"
		  << e1 << "\t" << e2 << "\t" << e3 << "\t" << e4 << "\t"
		  << w1 << "\t" << w2 << "\t" << w3 << "\t" << w4 << "\n";
      }
      
      else if (e1>4000. || e2>4000. || e3>4000. || e4>4000. || w1>4000. || w2>4000. || w3>4000. || w4>4000.) {

	std::cout << rn << "\t"
		  << e1 << "\t" << e2 << "\t" << e3 << "\t" << e4 << "\t"
		  << w1 << "\t" << w2 << "\t" << w3 << "\t" << w4 << "\n";
      }
      
      gainFile.close();
    }
  }


  //EAST SIDE
  TCanvas *cEast = new TCanvas("cEast","East PMTs",1000, 1200);
  cEast->Divide(1,4);

  cEast->cd(1);

  TGraph *east1 = new TGraph(runNumber[0].size(),&runNumber[0][0],&eastGain[0][0]);
  east1->SetTitle("PMT East 1");
  east1->GetXaxis()->SetTitle("Run Number");
  east1->GetYaxis()->SetTitle("ADC Channels");
  east1->SetMarkerStyle(7);
  // east1->SetMinimum(0.);
  //east1->SetMaximum(4096.);
  east1->Draw("AP");


  cEast->cd(2);

  TGraph *east2 = new TGraph(runNumber[1].size(),&runNumber[1][0],&eastGain[1][0]);
  east2->SetTitle("PMT East 2");
  east2->GetXaxis()->SetTitle("Run Number");
  east2->GetYaxis()->SetTitle("ADC Channels");
  east2->SetMarkerStyle(7);
  //east2->SetMinimum(0.);
  //east2->SetMaximum(4096.);
  east2->Draw("AP");


  cEast->cd(3);

  TGraph *east3 = new TGraph(runNumber[2].size(),&runNumber[2][0],&eastGain[2][0]);
  east3->SetTitle("PMT East 3");
  east3->GetXaxis()->SetTitle("Run Number");
  east3->GetYaxis()->SetTitle("ADC Channels");
  east3->SetMarkerStyle(7);
  //east3->SetMinimum(0.);
  //east3->SetMaximum(4096.);
  east3->Draw("AP");


  cEast->cd(4);

  TGraph *east4 = new TGraph(runNumber[3].size(),&runNumber[3][0],&eastGain[3][0]);
  east4->SetTitle("PMT East 4");
  east4->GetXaxis()->SetTitle("Run Number");
  east4->GetYaxis()->SetTitle("ADC Channels");
  east4->SetMarkerStyle(7);
  //east4->SetMinimum(0.);
  //east4->SetMaximum(4096.);
  east4->Draw("AP");

  //WEST SIDE
  TCanvas *cWest = new TCanvas("cWest","West PMTs",1000, 1200);
  cWest->Divide(1,4);

  cWest->cd(1);

  TGraph *west1 = new TGraph(runNumber[4].size(),&runNumber[4][0],&westGain[0][0]);
  west1->SetTitle("PMT West 1");
  west1->GetXaxis()->SetTitle("Run Number");
  west1->GetYaxis()->SetTitle("ADC Channels");
  west1->SetMarkerStyle(7);
  //west1->SetMinimum(0.);
  //west1->SetMaximum(4096.);
  west1->Draw("AP");


  cWest->cd(2);

  TGraph *west2 = new TGraph(runNumber[5].size(),&runNumber[5][0],&westGain[1][0]);
  west2->SetTitle("PMT West 2");
  west2->GetXaxis()->SetTitle("Run Number");
  west2->GetYaxis()->SetTitle("ADC Channels");
  west2->SetMarkerStyle(7);
  //west2->SetMinimum(0.);
  //west2->SetMaximum(4096.);
  west2->Draw("AP");


  cWest->cd(3);

  TGraph *west3 = new TGraph(runNumber[6].size(),&runNumber[6][0],&westGain[2][0]);
  west3->SetTitle("PMT West 3");
  west3->GetXaxis()->SetTitle("Run Number");
  west3->GetYaxis()->SetTitle("ADC Channels");
  west3->SetMarkerStyle(7);
  //west3->SetMinimum(0.);
  //west3->SetMaximum(4096.);
  west3->Draw("AP");


  cWest->cd(4);

  TGraph *west4 = new TGraph(runNumber[7].size(),&runNumber[7][0],&westGain[3][0]);
  west4->SetTitle("PMT West 4");
  west4->GetXaxis()->SetTitle("Run Number");
  west4->GetYaxis()->SetTitle("ADC Channels");
  west4->SetMarkerStyle(7);
  //west4->SetMinimum(0.);
  //west4->SetMaximum(4096.);
  west4->Draw("AP");

}
