#include <fstream>

void betaAsymmetryFitter() {

  bool color = false;
  UInt_t nskip=4;

  gStyle->SetTitleSize(0.08,"t");
  gStyle->SetTitleSize(0.06,"x");
  gStyle->SetTitleSize(0.07,"y");
  gStyle->SetStatX(0.75);
  gStyle->SetStatY(0.9);
  gStyle->SetStatBorderSize(0);
  //gStyle->SetOptStat(0);
  //gStyle->SetTitleOffset(0.85,"y");
  gStyle->SetPadTopMargin(0.0);
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetPadLeftMargin(0.12);
  gStyle->SetLabelSize(0.05,"xyz");
  

  gStyle->SetOptFit(1);//(1111);
  gStyle->SetStatY(0.85);
  gStyle->SetStatX(0.7);
  //gStyle->SetStatW(.09);
  gStyle->SetFillStyle(0000); 
  gStyle->SetLegendBorderSize(0); 
  
  
  TString octets2011 = "0-59";
  TString octets2012 = "60-121";

  //Rest of the corrections...
  Double_t effCorr2011 = 0.0013;
  Double_t effCorr2012 = 0.0011;
  Double_t neutronBGcorr = 0.0001;
  
  std::vector <Double_t> energy2011;
  std::vector <Double_t> AsymWithEnDep2011;
  std::vector <Double_t> AsymWithEnDepErr2011;
  std::vector <Double_t> Asym2011;
  std::vector <Double_t> AsymErr2011;

  std::vector <Double_t> energy2012;
  std::vector <Double_t> AsymWithEnDep2012;
  std::vector <Double_t> AsymWithEnDepErr2012;
  std::vector <Double_t> Asym2012;
  std::vector <Double_t> AsymErr2012;

  Int_t en = 0;
  Double_t a = 0.;
  Double_t aErr = 0.;

  // Load the beta/2 corrected data first
  TString fname = TString::Format("%s/Asymmetries/UNBLINDED_AllCorr_withPOL_OctetAsymmetries_AnaChC_Octets_%s_BinByBin.txt",getenv("ANALYSIS_RESULTS"),octets2011.Data());
  //UNBLINDED_

  std::ifstream infile(fname.Data());
  
  while ( infile >> en >> a >> aErr ) {
    //cout << en << "\t" << a << "\t" << aErr << "\t\n";
    energy2011.push_back(en);
    Asym2011.push_back(a*(1.+neutronBGcorr+effCorr2011));
    AsymErr2011.push_back(aErr*(1.+neutronBGcorr+effCorr2011));
  }
  
  infile.close();

  
  
  fname = TString::Format("%s/Asymmetries/UNBLINDED_AllCorr_withPOL_OctetAsymmetries_AnaChC_Octets_%s_BinByBin_withEnergyDependence.txt",getenv("ANALYSIS_RESULTS"),octets2011.Data());

  infile.open(fname.Data());
  
  while ( infile >> en >> a >> aErr ) {
    AsymWithEnDep2011.push_back(a*(1.+neutronBGcorr+effCorr2011));
    AsymWithEnDepErr2011.push_back(aErr*(1.+neutronBGcorr+effCorr2011));
  }
  
  infile.close();

  // Load the beta/2 corrected data first
  fname = TString::Format("%s/Asymmetries/UNBLINDED_AllCorr_withPOL_OctetAsymmetries_AnaChC_Octets_%s_BinByBin.txt",getenv("ANALYSIS_RESULTS"),octets2012.Data());
  //UNBLINDED_

  infile.open(fname.Data());
  
  while ( infile >> en >> a >> aErr ) {
    //cout << en << "\t" << a << "\t" << aErr << "\t\n";
    energy2012.push_back(en);
    Asym2012.push_back(a*(1.+neutronBGcorr+effCorr2012));
    AsymErr2012.push_back(aErr*(1.+neutronBGcorr+effCorr2012));
  }
  
  infile.close();

  
  
  fname = TString::Format("%s/Asymmetries/UNBLINDED_AllCorr_withPOL_OctetAsymmetries_AnaChC_Octets_%s_BinByBin_withEnergyDependence.txt",getenv("ANALYSIS_RESULTS"),octets2012.Data());

  infile.open(fname.Data());
  
  while ( infile >> en >> a >> aErr ) {
    AsymWithEnDep2012.push_back(a*(1.+neutronBGcorr+effCorr2012));
    AsymWithEnDepErr2012.push_back(aErr*(1.+neutronBGcorr+effCorr2012));
  }
  
  infile.close();
  
  TCanvas *c1 = new TCanvas("c1","c1",700,700);
  c1->Divide(1,2);

  c1->cd(1);
  TGraphErrors *gr1 = new TGraphErrors(energy2011.size(),&energy2011[0],&AsymWithEnDep2011[0],NULL,&AsymWithEnDepErr2011[0]);
  gr1->SetTitle(TString::Format("A(E) vs. Energy for 2011-2012").Data());
  gr1->GetXaxis()->SetTitle("Energy (keV)");
  gr1->GetYaxis()->SetTitle("A(E)");
  gr1->GetXaxis()->CenterTitle();
  gr1->GetYaxis()->CenterTitle();
  gr1->SetMarkerStyle(20);
  gr1->SetLineWidth(2);
  gr1->GetXaxis()->SetLimits(0., 800.);
  //gr1->SetMarkerSize(0.6);
  gr1->SetMinimum(-0.1);
  gr1->SetMaximum(0.);
  gr1->GetXaxis()->SetTitleOffset(1.);
  gr1->GetYaxis()->SetTitleOffset(1.2);
  
  Double_t A0 = -0.1184;
  TF1 *f1 = new TF1("f1","([0]/2.*sqrt(x*x+2.*x*510.99892)/(x+510.99892))",190.,740.);
  f1->SetParameter(0,-0.12);
  f1->SetParName(0,"A_{0}");

  //f1->SetParLimits(4,1.,20.);

  gr1->Fit("f1","LR");
  

  gr1->Draw("AP");
  
  TF1 *f2 = new TF1("f2","(-0.1184/2.*sqrt(x*x+2.*x*510.99892)/(x+510.99892))",0., 1000.);
  f2->SetLineColor(kBlack);
  f2->SetLineStyle(2);

  f2->Draw("SAME");

  c1->cd(2);

  TGraphErrors *gr1_2012 = new TGraphErrors(energy2012.size(),&energy2012[0],&AsymWithEnDep2012[0],NULL,&AsymWithEnDepErr2012[0]);
  gr1_2012->SetTitle(TString::Format("A(E) vs. Energy for 2012-2013").Data());
  gr1_2012->GetXaxis()->SetTitle("Energy (keV)");
  gr1_2012->GetYaxis()->SetTitle("A(E)");
  gr1_2012->GetXaxis()->CenterTitle();
  gr1_2012->GetYaxis()->CenterTitle();
  gr1_2012->SetMarkerStyle(20);
  gr1_2012->SetLineWidth(2);
  gr1_2012->GetXaxis()->SetLimits(0., 800.);
  //gr1_2012->SetMarkerSize(0.6);
  gr1_2012->SetMinimum(-0.1);
  gr1_2012->SetMaximum(0.);
  gr1_2012->GetXaxis()->SetTitleOffset(1.);
  gr1_2012->GetYaxis()->SetTitleOffset(1.2);
  
  TF1 *f1_2012 = new TF1("f1_2012","([0]/2.*sqrt(x*x+2.*x*510.99892)/(x+510.99892))",190.,740.);
  f1_2012->SetParameter(0,-0.12);
  f1_2012->SetParName(0,"A_{0}");

  //f1->SetParLimits(4,1.,20.);

  gr1_2012->Fit("f1","LR");
  

  gr1_2012->Draw("AP");
  f2->Draw("SAME");


  /////////////////////////////////////////////////////////////////////////////////////////
  TCanvas *c2 = new TCanvas("c2","c2",700,700);
  c2->Divide(1,2);
  c2->cd(1);

  TGraphErrors *gr2 = new TGraphErrors(energy2011.size()-nskip,&energy2011[nskip],&Asym2011[nskip],NULL,&AsymErr2011[nskip]);
  gr2->SetTitle("2011-2012");//TString::Format("A_{0} vs. Energy for 2011-2012").Data());
  gr2->GetXaxis()->SetTitle("Energy (keV)");
  gr2->GetYaxis()->SetTitle("A_{0}");
  gr2->GetXaxis()->CenterTitle();
  gr2->GetYaxis()->CenterTitle();
  gr2->SetMarkerStyle(20);
  gr2->SetLineWidth(2);
  gr2->GetXaxis()->SetLimits(0., 800.);
  gr2->SetMarkerSize(0.8);
  gr2->GetXaxis()->SetTitleOffset(1.);
  gr2->GetYaxis()->SetTitleOffset(0.8);
  
  TF1 *f3 = new TF1("f3","[0]",190.,740.);
  f3->SetParameter(0,-0.1184);
  if (!color) f3->SetLineColor(1);
  f3->SetLineStyle(2);
  f3->SetLineWidth(3);
  f3->SetParName(0,"A_{0}");

  gr2->Fit("f3","LR");
  gr2->SetMinimum(f3->GetParameter(0)-0.03);
  gr2->SetMaximum(f3->GetParameter(0)+0.045);
  
  gr2->Draw("AP0");

  c2->cd(2);

  TGraphErrors *gr2_2012 = new TGraphErrors(energy2012.size()-nskip,&energy2012[nskip],&Asym2012[nskip],NULL,&AsymErr2012[nskip]);
  gr2_2012->SetTitle("2012-2013");//TString::Format("A_{0} vs. Energy for 2012-2013").Data());
  gr2_2012->GetXaxis()->SetTitle("Energy (keV)");
  gr2_2012->GetYaxis()->SetTitle("A_{0}");
  gr2_2012->GetXaxis()->CenterTitle();
  gr2_2012->GetYaxis()->CenterTitle();
  gr2_2012->SetMarkerStyle(20);
  gr2_2012->SetLineWidth(2);
  gr2_2012->GetXaxis()->SetLimits(0., 800.);
  gr2_2012->SetMarkerSize(0.8);
  gr2_2012->GetXaxis()->SetTitleOffset(1.);
  gr2_2012->GetYaxis()->SetTitleOffset(0.8);
  
  TF1 *f3_2012 = new TF1("f3_2012","[0]",190.,740.);
  f3_2012->SetParameter(0,-0.1184);
  if (!color) f3_2012->SetLineColor(1);
  f3_2012->SetLineStyle(2);
  f3_2012->SetLineWidth(3);
  f3_2012->SetParName(0,"A_{0}");

  gr2_2012->Fit("f3_2012","LR");
  gr2_2012->SetMinimum(f3_2012->GetParameter(0)-0.03);
  gr2_2012->SetMaximum(f3_2012->GetParameter(0)+0.045);
  
  gr2_2012->Draw("AP0");

  c2->Print(TString::Format("AsymmetryVsEnergy%s.pdf",(color?"_color":"")));
  
  

}

  

 
