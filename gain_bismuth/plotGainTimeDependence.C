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
  gStyle->SetTitleYOffset(0.6);
  gStyle->SetTitleXSize(0.06);
  gStyle->SetTitleXOffset(0.9);
  gStyle->SetLabelSize(0.06,"xyz");
  
  std::vector <Int_t> calPeriods(10,0); // TO BE IMPLEMENTED LATER

  Int_t runMin, runMax;

  if (year==TString("2011-2012") ) {
    runMin = 16000;
    runMax = 20000;
  }
  else {
    runMin = 20000;
    runMax = 23500;
  }

  std::vector < std::vector <Double_t> > runNumber( 8, std::vector<Double_t>(4000) );
  std::vector < std::vector <Double_t> > eastGain( 4, std::vector<Double_t>(4000) );
  std::vector < std::vector <Double_t> > westGain( 4, std::vector<Double_t>(4000) );
  std::vector < std::vector <Double_t> > eastGainFactor( 4, std::vector<Double_t>(4000) );
  std::vector < std::vector <Double_t> > westGainFactor( 4, std::vector<Double_t>(4000) );
  std::vector < Int_t > numRuns(8,0);

  double e1=0., e2=0., e3=0., e4=0., w1=0., w2=0., w3=0., w4=0.;
  double e1_norm=0., e2_norm=0., e3_norm=0., e4_norm=0., w1_norm=0., w2_norm=0., w3_norm=0., w4_norm=0.;
  
  ifstream gainFile;

  for (Int_t rn=runMin; rn<runMax; rn++) {
    
    if (rn%500==0) std::cout << "NOW AT RUN " << rn << std::endl;

    TString filename = TString::Format("%s/gain_bismuth_%i.dat",getenv("GAIN_BISMUTH"),rn);
    gainFile.open(filename.Data());
    
    if (gainFile.is_open()) {
      if ( rn==17588 || (rn>19583 && rn<19616) || (rn>=16000 && rn<=17078) ) { gainFile.close(); continue; } //Groups of runs to ignore due to bad crap
      // 17588 is empty, the ranges are unused Xe runs which weren't processed, 
      // (rn>20515 && rn<21087) is all the garbage which isn't really usable from the beginning of 2012-2013

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
      
      if (pmtQuality[0]) { runNumber[0][numRuns[0]]=rn; eastGain[0][numRuns[0]]=e1; eastGainFactor[0][numRuns[0]]=e1_norm; numRuns[0]++; }
      else e1=1000.; //so it doesn't trigger the check below
      if (pmtQuality[1]) { runNumber[1][numRuns[1]]=rn; eastGain[1][numRuns[1]]=e2; eastGainFactor[1][numRuns[1]]=e2_norm; numRuns[1]++; }
      else e2=1000.; //so it doesn't trigger the check below
      if (pmtQuality[2]) { runNumber[2][numRuns[2]]=rn; eastGain[2][numRuns[2]]=e3; eastGainFactor[2][numRuns[2]]=e3_norm; numRuns[2]++; }
      else e3=1000.; //so it doesn't trigger the check below
      if (pmtQuality[3]) { runNumber[3][numRuns[3]]=rn; eastGain[3][numRuns[3]]=e4; eastGainFactor[3][numRuns[3]]=e4_norm; numRuns[3]++; }
      else e4=1000.; //so it doesn't trigger the check below
      if (pmtQuality[4]) { runNumber[4][numRuns[4]]=rn; westGain[0][numRuns[4]]=w1; westGainFactor[0][numRuns[4]]=w1_norm; numRuns[4]++; }
      else w1=1000.; //so it doesn't trigger the check below
      if (pmtQuality[5]) { runNumber[5][numRuns[5]]=rn; westGain[1][numRuns[5]]=w2; westGainFactor[1][numRuns[5]]=w2_norm; numRuns[5]++; }
      else w2=1000.; //so it doesn't trigger the check below
      if (pmtQuality[6]) { runNumber[6][numRuns[6]]=rn; westGain[2][numRuns[6]]=w3; westGainFactor[2][numRuns[6]]=w3_norm; numRuns[6]++; }
      else w3=1000.; //so it doesn't trigger the check below
      if (pmtQuality[7]) { runNumber[7][numRuns[7]]=rn; westGain[3][numRuns[7]]=w4; westGainFactor[3][numRuns[7]]=w4_norm; numRuns[7]++; }
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

  TGraph *east1 = new TGraph(numRuns[0],&runNumber[0][0],&eastGain[0][0]);
  east1->SetTitle("PMT East 1");
  east1->GetXaxis()->SetTitle("Run Number");
  east1->GetYaxis()->SetTitle("ADC Channels");
  east1->SetMarkerStyle(8);
  // east1->SetMinimum(0.);
  //east1->SetMaximum(4096.);
  east1->Draw("AP");


  cEast->cd(2);

  TGraph *east2 = new TGraph(numRuns[1],&runNumber[1][0],&eastGain[1][0]);
  east2->SetTitle("PMT East 2");
  east2->GetXaxis()->SetTitle("Run Number");
  east2->GetYaxis()->SetTitle("ADC Channels");
  east2->SetMarkerStyle(8);
  //east2->SetMinimum(0.);
  //east2->SetMaximum(4096.);
  east2->Draw("AP");


  cEast->cd(3);

  TGraph *east3 = new TGraph(numRuns[2],&runNumber[2][0],&eastGain[2][0]);
  east3->SetTitle("PMT East 3");
  east3->GetXaxis()->SetTitle("Run Number");
  east3->GetYaxis()->SetTitle("ADC Channels");
  east3->SetMarkerStyle(8);
  //east3->SetMinimum(0.);
  //east3->SetMaximum(4096.);
  east3->Draw("AP");


  cEast->cd(4);

  TGraph *east4 = new TGraph(numRuns[3],&runNumber[3][0],&eastGain[3][0]);
  east4->SetTitle("PMT East 4");
  east4->GetXaxis()->SetTitle("Run Number");
  east4->GetYaxis()->SetTitle("ADC Channels");
  east4->SetMarkerStyle(8);
  //east4->SetMinimum(0.);
  //east4->SetMaximum(4096.);
  east4->Draw("AP");

  //WEST SIDE
  TCanvas *cWest = new TCanvas("cWest","West PMTs",1000, 1200);
  cWest->Divide(1,4);

  cWest->cd(1);

  TGraph *west1 = new TGraph(numRuns[4],&runNumber[4][0],&westGain[0][0]);
  west1->SetTitle("PMT West 1");
  west1->GetXaxis()->SetTitle("Run Number");
  west1->GetYaxis()->SetTitle("ADC Channels");
  west1->SetMarkerStyle(8);
  //west1->SetMinimum(0.);
  //west1->SetMaximum(4096.);
  west1->Draw("AP");


  cWest->cd(2);

  TGraph *west2 = new TGraph(numRuns[5],&runNumber[5][0],&westGain[1][0]);
  west2->SetTitle("PMT West 2");
  west2->GetXaxis()->SetTitle("Run Number");
  west2->GetYaxis()->SetTitle("ADC Channels");
  west2->SetMarkerStyle(8);
  //west2->SetMinimum(0.);
  //west2->SetMaximum(4096.);
  west2->Draw("AP");


  cWest->cd(3);

  TGraph *west3 = new TGraph(numRuns[6],&runNumber[6][0],&westGain[2][0]);
  west3->SetTitle("PMT West 3");
  west3->GetXaxis()->SetTitle("Run Number");
  west3->GetYaxis()->SetTitle("ADC Channels");
  west3->SetMarkerStyle(8);
  //west3->SetMinimum(0.);
  //west3->SetMaximum(4096.);
  west3->Draw("AP");


  if (numRuns[7]>0) {
    cWest->cd(4);
    
    
    TGraph *west4 = new TGraph(numRuns[7],&runNumber[7][0],&westGain[3][0]);
    west4->SetTitle("PMT West 4");
    west4->GetXaxis()->SetTitle("Run Number");
    west4->GetYaxis()->SetTitle("ADC Channels");
    west4->SetMarkerStyle(8);
    //west4->SetMinimum(0.);
    //west4->SetMaximum(4096.);
    west4->Draw("AP");
    
  }
    
  
    ///////////////////////////// Plotting the gain factors //////////////////////////////////


    //EAST SIDE
  TCanvas *cEast_Fac = new TCanvas("cEast_Fac","East PMTs",1000, 1200);
  cEast_Fac->Divide(1,4);

  cEast_Fac->cd(1);

  TGraph *east_fac1 = new TGraph(numRuns[0],&runNumber[0][0],&eastGainFactor[0][0]);
  east_fac1->SetTitle("PMT East 1");
  east_fac1->GetXaxis()->SetTitle("Run Number");
  east_fac1->GetYaxis()->SetTitle("Gain Factor");
  east_fac1->SetMarkerStyle(8);
  // east_fac1->SetMinimum(0.);
  //east_fac1->SetMaximum(4096.);
  east_fac1->Draw("AP");


  cEast_Fac->cd(2);

  TGraph *east_fac2 = new TGraph(numRuns[1],&runNumber[1][0],&eastGainFactor[1][0]);
  east_fac2->SetTitle("PMT East 2");
  east_fac2->GetXaxis()->SetTitle("Run Number");
  east_fac2->GetYaxis()->SetTitle("ADC Channels");
  east_fac2->SetMarkerStyle(8);
  //east_fac2->SetMinimum(0.);
  //east_fac2->SetMaximum(4096.);
  east_fac2->Draw("AP");


  cEast_Fac->cd(3);

  TGraph *east_fac3 = new TGraph(numRuns[2],&runNumber[2][0],&eastGainFactor[2][0]);
  east_fac3->SetTitle("PMT East 3");
  east_fac3->GetXaxis()->SetTitle("Run Number");
  east_fac3->GetYaxis()->SetTitle("ADC Channels");
  east_fac3->SetMarkerStyle(8);
  //east_fac3->SetMinimum(0.);
  //east_fac3->SetMaximum(4096.);
  east_fac3->Draw("AP");


  cEast_Fac->cd(4);

  TGraph *east_fac4 = new TGraph(numRuns[3],&runNumber[3][0],&eastGainFactor[3][0]);
  east_fac4->SetTitle("PMT East 4");
  east_fac4->GetXaxis()->SetTitle("Run Number");
  east_fac4->GetYaxis()->SetTitle("ADC Channels");
  east_fac4->SetMarkerStyle(8);
  //east_fac4->SetMinimum(0.);
  //east_fac4->SetMaximum(4096.);
  east_fac4->Draw("AP");

  //WEST SIDE
  TCanvas *cWest_Fac = new TCanvas("cWest_Fac","West PMTs",1000, 1200);
  cWest_Fac->Divide(1,4);

  cWest_Fac->cd(1);

  TGraph *west_fac1 = new TGraph(numRuns[4],&runNumber[4][0],&westGainFactor[0][0]);
  west_fac1->SetTitle("PMT West 1");
  west_fac1->GetXaxis()->SetTitle("Run Number");
  west_fac1->GetYaxis()->SetTitle("ADC Channels");
  west_fac1->SetMarkerStyle(8);
  //west_fac1->SetMinimum(0.);
  //west_fac1->SetMaximum(4096.);
  west_fac1->Draw("AP");


  cWest_Fac->cd(2);

  TGraph *west_fac2 = new TGraph(numRuns[5],&runNumber[5][0],&westGainFactor[1][0]);
  west_fac2->SetTitle("PMT West 2");
  west_fac2->GetXaxis()->SetTitle("Run Number");
  west_fac2->GetYaxis()->SetTitle("ADC Channels");
  west_fac2->SetMarkerStyle(8);
  //west_fac2->SetMinimum(0.);
  //west_fac2->SetMaximum(4096.);
  west_fac2->Draw("AP");


  cWest_Fac->cd(3);

  TGraph *west_fac3 = new TGraph(numRuns[6],&runNumber[6][0],&westGainFactor[2][0]);
  west_fac3->SetTitle("PMT West 3");
  west_fac3->GetXaxis()->SetTitle("Run Number");
  west_fac3->GetYaxis()->SetTitle("ADC Channels");
  west_fac3->SetMarkerStyle(8);
  //west_fac3->SetMinimum(0.);
  //west_fac3->SetMaximum(4096.);
  west_fac3->Draw("AP");


  if (numRuns[7]>0) {
    cWest_Fac->cd(4);
    
    
    TGraph *west_fac4 = new TGraph(numRuns[7],&runNumber[7][0],&westGainFactor[3][0]);
    west_fac4->SetTitle("PMT West 4");
    west_fac4->GetXaxis()->SetTitle("Run Number");
    west_fac4->GetYaxis()->SetTitle("ADC Channels");
    west_fac4->SetMarkerStyle(8);
    //west_fac4->SetMinimum(0.);
    //west_fac4->SetMaximum(4096.);
    west_fac4->Draw("AP");
    
  }
    

  TString fnamebase = TString::Format("%s_gain.pdf",year.Data());
  cEast->Print(TString::Format("%s(",fnamebase.Data()));
  cWest->Print(TString::Format("%s",fnamebase.Data()));
  cEast_Fac->Print(TString::Format("%s",fnamebase.Data()));
  cWest_Fac->Print(TString::Format("%s)",fnamebase.Data()));
  

}
