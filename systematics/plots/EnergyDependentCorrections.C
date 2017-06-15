

#include <vector>

void EnergyDependentCorrections() {

  gStyle->SetOptStat(0);
  gStyle->SetTitleSize(0.08,"t");
  gStyle->SetPadLeftMargin(0.15);
  //gStyle->SetPadRightMargin(0.15);
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetTitleYSize(0.05);
  gStyle->SetTitleYOffset(1.4);
  gStyle->SetTitleXSize(0.05);
  gStyle->SetTitleXOffset(1.2);
  gStyle->SetTitleSize(0.04,"Z");
  gStyle->SetTitleOffset(1.6,"Z");
  gStyle->SetLabelSize(0.04,"xyz");
  gStyle->SetLegendBorderSize(0);
  gStyle->SetFillStyle(0);

  // read in the corrections and plot their deltaA/A for each bin
  // This can be done for both theory and MC

  //Int_t octmin = year==2011?0:60;
  //Int_t octmax = year==2011?59:121;

  std::vector <Double_t> energy;
  std::vector <Double_t> UnCorr2011;
  std::vector <Double_t> TheoryCorr2011;
  std::vector <Double_t> AllCorr2011;

  std::vector <Double_t> UnCorr2012;
  std::vector <Double_t> TheoryCorr2012;
  std::vector <Double_t> AllCorr2012;

  Double_t en, val, err;
  
  std::ifstream infile(TString::Format("%s/Asymmetries/UnCorr_OctetAsymmetries_AnaChC_Octets_%i-%i_BinByBin.txt",getenv("ANALYSIS_RESULTS"),0,59));
  
  while (infile >> en >> val >> err ) {
    energy.push_back(en);
    UnCorr2011.push_back(val);
  }

  infile.close();

  infile.open(TString::Format("%s/Asymmetries/DeltaTheoryOnly_OctetAsymmetries_AnaChC_Octets_%i-%i_BinByBin.txt",getenv("ANALYSIS_RESULTS"),0,59));
  
  while (infile >> en >> val >> err) {
    TheoryCorr2011.push_back(val);
  }

  infile.close();

  infile.open(TString::Format("%s/Asymmetries/AllCorr_OctetAsymmetries_AnaChC_Octets_%i-%i_BinByBin.txt",getenv("ANALYSIS_RESULTS"),0,59));
  
  while (infile >> en >> val >> err) {
    AllCorr2011.push_back(val);
  }

  infile.close();
    
  infile.open(TString::Format("%s/Asymmetries/UnCorr_OctetAsymmetries_AnaChC_Octets_%i-%i_BinByBin.txt",getenv("ANALYSIS_RESULTS"),60,121));
  
  while (infile >> en >> val >> err ) {
    energy.push_back(en);
    UnCorr2012.push_back(val);
  }

  infile.close();

  infile.open(TString::Format("%s/Asymmetries/DeltaTheoryOnly_OctetAsymmetries_AnaChC_Octets_%i-%i_BinByBin.txt",getenv("ANALYSIS_RESULTS"),60,121));
  
  while (infile >> en >> val >> err) {
    TheoryCorr2012.push_back(val);
  }

  infile.close();

  infile.open(TString::Format("%s/Asymmetries/AllCorr_OctetAsymmetries_AnaChC_Octets_%i-%i_BinByBin.txt",getenv("ANALYSIS_RESULTS"),60,121));
  
  while (infile >> en >> val >> err) {
    AllCorr2012.push_back(val);
  }

  infile.close();


  TH1D *theory = new TH1D("theory","Theory Corrections vs. Energy",
			  80,0.,800.);
  theory->SetMarkerStyle(kFullCircle);
  theory->SetMinimum(-5.);
  theory->SetMaximum(0.);
  theory->SetLineWidth(3);
  theory->SetLineColor(kBlack);
  theory->GetYaxis()->SetTitle("#DeltaA/A (%)");
  theory->GetYaxis()->CenterTitle();
  theory->GetXaxis()->SetTitle("Energy (keV)");
  theory->GetXaxis()->CenterTitle();
  
  
  TH1D *mc2011 = new TH1D("mc2011","MC Corrections vs. Energy",
			  80,0.,800.);
  mc2011->SetMarkerStyle(kFullCircle);
  mc2011->SetMinimum(-10.);
  mc2011->SetMaximum(10.);
  mc2011->SetLineWidth(3);
  mc2011->SetLineColor(kBlue);
  mc2011->GetYaxis()->SetTitle("#DeltaA/A (%)");
  mc2011->GetYaxis()->CenterTitle();
  mc2011->GetXaxis()->SetTitle("Energy (keV)");
  mc2011->GetXaxis()->CenterTitle();
  

  TH1D *mc2012 = new TH1D("mc2012","MC Corrections vs. Energy",
			  80,0.,800.);
  mc2012->SetMarkerStyle(kFullCircle);
  mc2012->SetMinimum(-10.);
  mc2012->SetMaximum(10.);
  mc2012->SetLineWidth(3);
  mc2012->SetLineColor(kGreen);
  mc2012->GetYaxis()->SetTitle("#DeltaA/A (%)");
  mc2012->GetYaxis()->CenterTitle();
  mc2012->GetXaxis()->SetTitle("Energy (keV)");
  mc2012->GetXaxis()->CenterTitle();
  
  TLegend *leg = new TLegend(0.55,0.75,0.85,0.85);
  leg->AddEntry(mc2011,"2011-2012","l");
  leg->AddEntry(mc2012,"2012-2013","l");
  
  
  for (int i=2;i<80;++i) {

    theory->SetBinContent(i+1,UnCorr2011[i]==0.?0.:
			  (100.*(TheoryCorr2011[i]/UnCorr2011[i]-1.)));
    mc2011->SetBinContent(i+1,TheoryCorr2011[i]==0.?0.:
		      (100.*(AllCorr2011[i]/TheoryCorr2011[i]-1.)));
    mc2012->SetBinContent(i+1,TheoryCorr2012[i]==0.?0.:
		      (100.*(AllCorr2012[i]/TheoryCorr2012[i]-1.)));
  }

  TCanvas *c1 = new TCanvas("c1","c1",1200,800);
  c1->Divide(2,1);

  c1->cd(1);

  theory->Draw("C0");
  
  c1->cd(2);

  mc2011->Draw("C0");
  mc2012->Draw("C0 SAME");
  leg->Draw("SAME");
}
