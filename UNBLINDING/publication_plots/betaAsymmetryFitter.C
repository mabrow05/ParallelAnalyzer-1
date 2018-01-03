#include <fstream>

void betaAsymmetryFitter() {

  bool color = false;
  UInt_t nskip=4;

  //gStyle->SetTitleSize(0.1,"t");
  //gStyle->SetTitleSize(0.1,"x");
  //gStyle->SetTitleSize(0.1,"y");
  gStyle->SetOptStat(0);
  //gStyle->SetStatX(0.75);
  //gStyle->SetStatY(0.9);
  gStyle->SetStatBorderSize(0);
  //gStyle->SetOptStat(0);
  //gStyle->SetTitleOffset(0.85,"y");
  gStyle->SetPadTopMargin(0.06);
  gStyle->SetPadRightMargin(0.05);
  gStyle->SetPadBottomMargin(0.11);
  gStyle->SetPadLeftMargin(0.12);
  //gStyle->SetLabelSize(0.1,"xyz");
  

  gStyle->SetOptFit(0);//(1111);
  gStyle->SetStatY(0.85);
  gStyle->SetStatX(0.7);
  //gStyle->SetStatW(.09);
  gStyle->SetFillStyle(0000); 
  gStyle->SetLegendBorderSize(0);
  gStyle->SetErrorX(0);

  TH1::SetDefaultSumw2(true);
  
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
  gr1->GetXaxis()->SetTitle("E_{recon} (keV)");
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
  gr1_2012->GetXaxis()->SetTitle("E_{recon} (keV)");
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
  TCanvas *c2 = new TCanvas("c2","c2",700,900);
  TPad *p1 = new TPad("p1","p1",0.,0.6,1.,1.);
  TPad *p2 = new TPad("p2","p2",0.,0.42,1.,0.6);
  TPad *p3 = new TPad("p3","p3",0.,0.24,1.,0.42);
  TPad *p4 = new TPad("p4","p4",0.,0.06,1.,0.24);
  //p2->SetBottomMargin(0.3);
  p1->Draw();
  p2->Draw();
  p3->Draw();
  p4->Draw();

  p3->cd();

  TGraphErrors *gr2 = new TGraphErrors(energy2011.size()-nskip,&energy2011[nskip],&Asym2011[nskip],NULL,&AsymErr2011[nskip]);
  gr2->SetTitle("");//TString::Format("A_{0} vs. Energy for 2011-2012").Data());
  //gr2->GetXaxis()->SetTitle("Energy (keV)");
  gr2->GetYaxis()->SetTitle("A_{0}");
  gr2->GetYaxis()->SetLabelSize(0.12);
  gr2->GetYaxis()->SetTitleOffset(0.32);
  gr2->GetYaxis()->SetTitleSize(0.18);
  gr2->GetXaxis()->SetLabelSize(0.12);
  gr2->GetXaxis()->CenterTitle();
  gr2->GetYaxis()->CenterTitle();
  gr2->SetMarkerStyle(20);
  gr2->SetLineWidth(2);
  gr2->GetXaxis()->SetLimits(0., 800.);
  gr2->GetYaxis()->SetNdivisions(508);
  gr2->SetMarkerSize(1);
  gr2->GetXaxis()->SetTitleOffset(1.);
  
  
  TF1 *f3 = new TF1("f3","[0]",190.,740.);
  f3->SetParameter(0,-0.1184);
  if (!color) f3->SetLineColor(1);
  f3->SetLineStyle(1);
  f3->SetLineWidth(3);
  f3->SetParName(0,"A_{0}");

  
  gr2->Fit("f3","LR");
  gr2->SetMinimum(-0.138);//f3->GetParameter(0)-0.02);
  gr2->SetMaximum(-0.094);//f3->GetParameter(0)+0.03);
  
  gr2->Draw("AP0Z");

  TPaveText *pv2011 = new TPaveText(0.3,0.67,0.7,0.95,"nbNDC");
  pv2011->SetBorderSize(0);
  pv2011->AddText(TString::Format("2011-2012: A_{0} = %0.5f #pm %0.5f",
				  f3->GetParameter(0),f3->GetParError(0)));
  pv2011->GetLine(0)->SetTextSize(0.13);
  pv2011->GetLine(0)->SetTextFont(42);
  pv2011->Draw();

  p4->cd();
  //p3->SetBottomMargin(0.25);

  TGraphErrors *gr2_2012 = new TGraphErrors(energy2012.size()-nskip,&energy2012[nskip],&Asym2012[nskip],NULL,&AsymErr2012[nskip]);
  gr2_2012->SetTitle("");//TString::Format("A_{0} vs. Energy for 2012-2013").Data());
  //gr2_2012->GetXaxis()->SetTitle("Energy (keV)");
  gr2_2012->GetYaxis()->SetTitle("A_{0}");
  gr2_2012->GetYaxis()->SetLabelSize(0.12);
  gr2_2012->GetYaxis()->SetTitleOffset(0.32);
  gr2_2012->GetYaxis()->SetTitleSize(0.18);
  gr2_2012->GetXaxis()->SetLabelSize(0.12);
  gr2_2012->GetXaxis()->CenterTitle();
  gr2_2012->GetYaxis()->CenterTitle();
  gr2_2012->SetMarkerStyle(20);
  gr2_2012->SetLineWidth(2);
  gr2_2012->GetXaxis()->SetLimits(0., 800.);
  gr2_2012->GetYaxis()->SetNdivisions(508);
  gr2_2012->SetMarkerSize(1);
  gr2_2012->GetXaxis()->SetTitleOffset(1.);
  
  
  TF1 *f3_2012 = new TF1("f3_2012","[0]",190.,740.);
  f3_2012->SetParameter(0,-0.1184);
  if (!color) f3_2012->SetLineColor(1);
  f3_2012->SetLineStyle(1);
  f3_2012->SetLineWidth(3);
  f3_2012->SetParName(0,"A_{0}");

  gr2_2012->Fit("f3_2012","LR");
  gr2_2012->SetMinimum(-0.138);
  gr2_2012->SetMaximum(-0.094);
  //gr2_2012->SetMinimum(f3_2012->GetParameter(0)-0.02);
  //gr2_2012->SetMaximum(f3_2012->GetParameter(0)+0.03);
  
  gr2_2012->Draw("AP0Z");


  TPaveText *pv2012 = new TPaveText(0.3,0.67,0.7,0.95,"nbNDC");
  pv2012->SetBorderSize(0);
  pv2012->AddText(TString::Format("2012-2013: A_{0} = %0.5f #pm %0.5f",
				  f3_2012->GetParameter(0),f3_2012->GetParError(0)));
  pv2012->GetLine(0)->SetTextSize(0.13);
  pv2012->GetLine(0)->SetTextFont(42);
  pv2012->Draw();

  
  
  
  ////////////////////////////////////////////////////////////////////////////
  
  Double_t normLow = 190.;
  Double_t normHigh = 740.;
  
  Double_t xAxisMax = 800.;

  Double_t minEnBin = 19; //Energy bins to integrate event types over
  Double_t maxEnBin = 73;

  //Storing event fractions for data, E/W and sfON/OFF
  Double_t t0E_sfON, t0E_sfOFF, t1E_sfON, t1E_sfOFF, t2E_sfON, t2E_sfOFF, t3E_sfON, t3E_sfOFF;
  Double_t t0W_sfON, t0W_sfOFF, t1W_sfON, t1W_sfOFF, t2W_sfON, t2W_sfOFF, t3W_sfON, t3W_sfOFF;
  t0E_sfON = t0E_sfOFF = t1E_sfON = t1E_sfOFF = t2E_sfOFF = t3E_sfOFF = 0;
  t0W_sfON = t0W_sfOFF = t1W_sfON = t1W_sfOFF = t2W_sfOFF = t3E_sfOFF = 0.;

  //Storing event fractions for sim, E/W and sfON/OFF
  Double_t sim_t0E_sfON, sim_t0E_sfOFF, sim_t1E_sfON, sim_t1E_sfOFF, sim_t2E_sfON, sim_t2E_sfOFF, sim_t3E_sfON, sim_t3E_sfOFF;
  Double_t sim_t0W_sfON, sim_t0W_sfOFF, sim_t1W_sfON, sim_t1W_sfOFF, sim_t2W_sfON, sim_t2W_sfOFF, sim_t3W_sfON, sim_t3W_sfOFF;
  sim_t0E_sfON = sim_t0E_sfOFF = sim_t1E_sfON = sim_t1E_sfOFF  = sim_t2E_sfOFF = sim_t3E_sfOFF = 0;
  sim_t0W_sfON = sim_t0W_sfOFF = sim_t1W_sfON = sim_t1W_sfOFF =  sim_t2W_sfOFF = sim_t3E_sfOFF = 0.;

  
  TFile *data_file = new TFile(TString::Format("../../Asymmetry/SpectraComparisons/Octets_0-121_DATA.root"),"READ");
  TFile *sim_file = new TFile(TString::Format("../../Asymmetry/SpectraComparisons/Octets_0-121_SIM.root"),"READ");


  //Load the histograms

  TH1D *dataALL = (TH1D*)data_file->Get("EreconALL");
  TH1D *data0 = (TH1D*)data_file->Get("Erecon0");
  TH1D *data1 = (TH1D*)data_file->Get("Erecon1");
  TH1D *data2 = (TH1D*)data_file->Get("Erecon2");
  TH1D *data3 = (TH1D*)data_file->Get("Erecon3");

  TH1D *BG_dataALL = (TH1D*)data_file->Get("BG_EreconALL");
  TH1D *BG_data0 = (TH1D*)data_file->Get("BG_Erecon0");
  TH1D *BG_data1 = (TH1D*)data_file->Get("BG_Erecon1");
  TH1D *BG_data2 = (TH1D*)data_file->Get("BG_Erecon2");
  TH1D *BG_data3 = (TH1D*)data_file->Get("BG_Erecon3");
  
  TH1D *simALL = (TH1D*)sim_file->Get("EreconALL");
  TH1D *sim0 = (TH1D*)sim_file->Get("Erecon0");
  TH1D *sim1 = (TH1D*)sim_file->Get("Erecon1");
  TH1D *sim2 = (TH1D*)sim_file->Get("Erecon2");
  TH1D *sim3 = (TH1D*)sim_file->Get("Erecon3");
  
  TH1D *dataALL_sfOFF_E = (TH1D*)data_file->Get("EreconALL_sfOFF_E");
  TH1D *data0_sfOFF_E = (TH1D*)data_file->Get("Erecon0_sfOFF_E");
  TH1D *data1_sfOFF_E = (TH1D*)data_file->Get("Erecon1_sfOFF_E");
  TH1D *data2_sfOFF_E = (TH1D*)data_file->Get("Erecon2_sfOFF_E");
  TH1D *data3_sfOFF_E = (TH1D*)data_file->Get("Erecon3_sfOFF_E");

  TH1D *BG_dataALL_sfOFF_E = (TH1D*)data_file->Get("BG_EreconALL_sfOFF_E");
  TH1D *BG_data0_sfOFF_E = (TH1D*)data_file->Get("BG_Erecon0_sfOFF_E");
  TH1D *BG_data1_sfOFF_E = (TH1D*)data_file->Get("BG_Erecon1_sfOFF_E");
  TH1D *BG_data2_sfOFF_E = (TH1D*)data_file->Get("BG_Erecon2_sfOFF_E");
  TH1D *BG_data3_sfOFF_E = (TH1D*)data_file->Get("BG_Erecon3_sfOFF_E");
  
  TH1D *simALL_sfOFF_E = (TH1D*)sim_file->Get("EreconALL_sfOFF_E");
  TH1D *sim0_sfOFF_E = (TH1D*)sim_file->Get("Erecon0_sfOFF_E");
  TH1D *sim1_sfOFF_E = (TH1D*)sim_file->Get("Erecon1_sfOFF_E");
  TH1D *sim2_sfOFF_E = (TH1D*)sim_file->Get("Erecon2_sfOFF_E");
  TH1D *sim3_sfOFF_E = (TH1D*)sim_file->Get("Erecon3_sfOFF_E");

  TH1D *dataALL_sfOFF_W = (TH1D*)data_file->Get("EreconALL_sfOFF_W");
  TH1D *data0_sfOFF_W = (TH1D*)data_file->Get("Erecon0_sfOFF_W");
  TH1D *data1_sfOFF_W = (TH1D*)data_file->Get("Erecon1_sfOFF_W");
  TH1D *data2_sfOFF_W = (TH1D*)data_file->Get("Erecon2_sfOFF_W");
  TH1D *data3_sfOFF_W = (TH1D*)data_file->Get("Erecon3_sfOFF_W");

  TH1D *BG_dataALL_sfOFF_W = (TH1D*)data_file->Get("BG_EreconALL_sfOFF_W");
  TH1D *BG_data0_sfOFF_W = (TH1D*)data_file->Get("BG_Erecon0_sfOFF_W");
  TH1D *BG_data1_sfOFF_W = (TH1D*)data_file->Get("BG_Erecon1_sfOFF_W");
  TH1D *BG_data2_sfOFF_W = (TH1D*)data_file->Get("BG_Erecon2_sfOFF_W");
  TH1D *BG_data3_sfOFF_W = (TH1D*)data_file->Get("BG_Erecon3_sfOFF_W");
  
  TH1D *simALL_sfOFF_W = (TH1D*)sim_file->Get("EreconALL_sfOFF_W");
  TH1D *sim0_sfOFF_W = (TH1D*)sim_file->Get("Erecon0_sfOFF_W");
  TH1D *sim1_sfOFF_W = (TH1D*)sim_file->Get("Erecon1_sfOFF_W");
  TH1D *sim2_sfOFF_W = (TH1D*)sim_file->Get("Erecon2_sfOFF_W");
  TH1D *sim3_sfOFF_W = (TH1D*)sim_file->Get("Erecon3_sfOFF_W");
  
  TH1D *dataALL_sfON_E = (TH1D*)data_file->Get("EreconALL_sfON_E");
  TH1D *data0_sfON_E = (TH1D*)data_file->Get("Erecon0_sfON_E");
  TH1D *data1_sfON_E = (TH1D*)data_file->Get("Erecon1_sfON_E");
  TH1D *data2_sfON_E = (TH1D*)data_file->Get("Erecon2_sfON_E");
  TH1D *data3_sfON_E = (TH1D*)data_file->Get("Erecon3_sfON_E");

  TH1D *BG_dataALL_sfON_E = (TH1D*)data_file->Get("BG_EreconALL_sfON_E");
  TH1D *BG_data0_sfON_E = (TH1D*)data_file->Get("BG_Erecon0_sfON_E");
  TH1D *BG_data1_sfON_E = (TH1D*)data_file->Get("BG_Erecon1_sfON_E");
  TH1D *BG_data2_sfON_E = (TH1D*)data_file->Get("BG_Erecon2_sfON_E");
  TH1D *BG_data3_sfON_E = (TH1D*)data_file->Get("BG_Erecon3_sfON_E");
  
  TH1D *simALL_sfON_E = (TH1D*)sim_file->Get("EreconALL_sfON_E");
  TH1D *sim0_sfON_E = (TH1D*)sim_file->Get("Erecon0_sfON_E");
  TH1D *sim1_sfON_E = (TH1D*)sim_file->Get("Erecon1_sfON_E");
  TH1D *sim2_sfON_E = (TH1D*)sim_file->Get("Erecon2_sfON_E");
  TH1D *sim3_sfON_E = (TH1D*)sim_file->Get("Erecon3_sfON_E");

  TH1D *dataALL_sfON_W = (TH1D*)data_file->Get("EreconALL_sfON_W");
  TH1D *data0_sfON_W = (TH1D*)data_file->Get("Erecon0_sfON_W");
  TH1D *data1_sfON_W = (TH1D*)data_file->Get("Erecon1_sfON_W");
  TH1D *data2_sfON_W = (TH1D*)data_file->Get("Erecon2_sfON_W");
  TH1D *data3_sfON_W = (TH1D*)data_file->Get("Erecon3_sfON_W");

  TH1D *BG_dataALL_sfON_W = (TH1D*)data_file->Get("BG_EreconALL_sfON_W");
  TH1D *BG_data0_sfON_W = (TH1D*)data_file->Get("BG_Erecon0_sfON_W");
  TH1D *BG_data1_sfON_W = (TH1D*)data_file->Get("BG_Erecon1_sfON_W");
  TH1D *BG_data2_sfON_W = (TH1D*)data_file->Get("BG_Erecon2_sfON_W");
  TH1D *BG_data3_sfON_W = (TH1D*)data_file->Get("BG_Erecon3_sfON_W");
  
  TH1D *simALL_sfON_W = (TH1D*)sim_file->Get("EreconALL_sfON_W");
  TH1D *sim0_sfON_W = (TH1D*)sim_file->Get("Erecon0_sfON_W");
  TH1D *sim1_sfON_W = (TH1D*)sim_file->Get("Erecon1_sfON_W");
  TH1D *sim2_sfON_W = (TH1D*)sim_file->Get("Erecon2_sfON_W");
  TH1D *sim3_sfON_W = (TH1D*)sim_file->Get("Erecon3_sfON_W");
  

  t0E_sfOFF = data0_sfOFF_E->Integral(minEnBin,maxEnBin);
  t0E_sfON = data0_sfON_E->Integral(minEnBin,maxEnBin);
  t1E_sfOFF = data1_sfOFF_E->Integral(minEnBin,maxEnBin);
  t1E_sfON = data1_sfON_E->Integral(minEnBin,maxEnBin);
  t2E_sfOFF = data2_sfOFF_E->Integral(minEnBin,maxEnBin);
  t2E_sfON = data2_sfON_E->Integral(minEnBin,maxEnBin);
  t3E_sfOFF = data3_sfOFF_E->Integral(minEnBin,maxEnBin);
  t3E_sfON = data3_sfON_E->Integral(minEnBin,maxEnBin);
  t0W_sfOFF = data0_sfOFF_W->Integral(minEnBin,maxEnBin);
  t0W_sfON = data0_sfON_W->Integral(minEnBin,maxEnBin);
  t1W_sfOFF = data1_sfOFF_W->Integral(minEnBin,maxEnBin);
  t1W_sfON = data1_sfON_W->Integral(minEnBin,maxEnBin);
  t2W_sfOFF = data2_sfOFF_W->Integral(minEnBin,maxEnBin);
  t2W_sfON = data2_sfON_W->Integral(minEnBin,maxEnBin);
  t3W_sfOFF = data3_sfOFF_W->Integral(minEnBin,maxEnBin);
  t3W_sfON = data3_sfON_W->Integral(minEnBin,maxEnBin);
  

  double totalE_sfOFF = t0E_sfOFF;// + t1E_sfOFF + t2E_sfOFF + t3E_sfOFF;
  double totalE_sfON = t0E_sfON;// + t1E_sfON + t2E_sfON + t3E_sfON;
  double totalW_sfOFF = t0W_sfOFF;// + t1W_sfOFF + t2W_sfOFF + t3W_sfOFF;
  double totalW_sfON = t0W_sfON;// + t1W_sfON + t2W_sfON + t3W_sfON;

  double fracData1 = 100.*data1->Integral(minEnBin,maxEnBin)/data0->Integral(minEnBin,maxEnBin);
  double fracData2 = 100.*data2->Integral(minEnBin,maxEnBin)/data0->Integral(minEnBin,maxEnBin);
  double fracData3 = 100.*data3->Integral(minEnBin,maxEnBin)/data0->Integral(minEnBin,maxEnBin);  
  double fracSim1 = 100.*sim1->Integral(minEnBin,maxEnBin)/sim0->Integral(minEnBin,maxEnBin);
  double fracSim2 = 100.*sim2->Integral(minEnBin,maxEnBin)/sim0->Integral(minEnBin,maxEnBin);
  double fracSim3 = 100.*sim3->Integral(minEnBin,maxEnBin)/sim0->Integral(minEnBin,maxEnBin);  
  
  std::cout << std::fixed;
  std::cout << std::setprecision(2);
  std::cout << "********************************************\n"
	    << "          Event ratios (220-670 keV)\n\n"
	    << "\t\tsfON\t\t\tsfOFF\n"
	    << "Type 0E:\t" << t0E_sfON/totalE_sfON << "\t\t" << t0E_sfOFF/totalE_sfOFF << "\n"
	    << "Type 0W:\t" << t0W_sfON/totalW_sfON << "\t\t" << t0W_sfOFF/totalW_sfOFF << "\n\n"
	    << "Type 1E:\t" << t1E_sfON/totalE_sfON << "\t\t" << t1E_sfOFF/totalE_sfOFF << "\n"
	    << "Type 1W:\t" << t1W_sfON/totalW_sfON << "\t\t" << t1W_sfOFF/totalW_sfOFF << "\n\n"
	    << "Type 2E:\t" << t2E_sfON/totalE_sfON << "\t\t" << t2E_sfOFF/totalE_sfOFF << "\n"
	    << "Type 2W:\t" << t2W_sfON/totalW_sfON << "\t\t" << t2W_sfOFF/totalW_sfOFF << "\n\n"
	    << "Type 3E:\t" << t3E_sfON/totalE_sfON << "\t\t" << t3E_sfOFF/totalE_sfOFF << "\n"
	    << "Type 3W:\t" << t3W_sfON/totalW_sfON << "\t\t" << t3W_sfOFF/totalW_sfOFF << "\n\n"
	    << "********************************************\n\n"
	    << "          Total Fractions (as % of type0)\n"
	    << "Type\t\t1\t2\t3\n"
	    << "Data\t\t"<<fracData1<<"\t"<<fracData2<<"\t"<<fracData3<<"\n"
	    << "Sim\t\t"<<fracSim1<<"\t"<<fracSim2<<"\t"<<fracSim3<<"\n"
	    << "% Diff\t\t"
	    << (fracSim1-fracData1)/fracData1*100.<<"\t"
	    << (fracSim2-fracData2)/fracData2*100.<<"\t"
	    << (fracSim3-fracData3)/fracData3*100.<<"\n\n***************************************\n\n";
 
    
  //Normalize
  Double_t normFactor = 1.;

  normFactor = dataALL->Integral(dataALL->GetXaxis()->FindFixBin(normLow),dataALL->GetXaxis()->FindFixBin(normHigh))/simALL->Integral(simALL->GetXaxis()->FindFixBin(normLow),simALL->GetXaxis()->FindFixBin(normHigh));
  
  
  simALL->Scale(normFactor);
  

  ////////////////////// SUPER_SUMS ////////////////////////////////////////////

  //////// Super-Sum of all data ///////

  // All Types

  /*
  //Residuals...
  TCanvas *c1 = new TCanvas("c1","c1", 500, 600);
  //c1->Divide(1,2);
  //c1->cd(2);
  TPad *p1 = new TPad("p1","p1",0.,0.3,1.,1.);
  TPad *p2 = new TPad("p2","p2",0.,0.,1.,0.3);
  p2->SetBottomMargin(0.3);
  p1->Draw();
  p2->Draw();

  p2->cd();
  Int_t nBins = dataALL->GetNbinsX();
  Double_t Min = dataALL->GetXaxis()->GetBinLowEdge(dataALL->GetXaxis()->GetFirst());
  Double_t Max = dataALL->GetXaxis()->GetBinUpEdge(dataALL->GetXaxis()->GetLast());

  //residual
  TH1D *residALL = new TH1D("residALL","Residuals: MC-Data", nBins, Min, Max);
  residALL->Add(simALL,dataALL,1,-1);
  residALL->GetXaxis()->SetRangeUser(0., xAxisMax);
  residALL->GetXaxis()->SetTitleOffset(1.);
  residALL->GetYaxis()->SetTitleOffset(0.5);
  residALL->GetYaxis()->SetTitle("event rate (mHz/keV)");
  residALL->GetXaxis()->SetTitle("Energy (keV)");
  //residALL->SetMaximum(2.);
  //residALL->SetMinimum(-2.);
  residALL->SetLineWidth(2);
  residALL->SetMarkerColor(color?kBlue:kBlack);
  residALL->SetMarkerStyle(34);
  residALL->SetMarkerSize(1.);

  TH1D *perc_residALL = new TH1D("perc_residALL","Fractional Residuals: (MC-Data)/MC", nBins, Min, Max);
  perc_residALL->Add(simALL,dataALL,1,-1);
  perc_residALL->Divide(perc_residALL,dataALL,1.,1.);
  perc_residALL->GetXaxis()->SetRangeUser(0., xAxisMax);
  //perc_residALL->GetYaxis()->SetTitle("event rate (mHz/keV)");
  perc_residALL->SetMaximum(.1);
  perc_residALL->SetMinimum(-.1);
  perc_residALL->SetLineWidth(2);
  


  //perc_residALL->Draw();
  residALL->Draw("E0");
  //c1->Update();
  TLine *l = new TLine(Min, 0., xAxisMax, 0.);
  l->SetLineStyle(8);
  l->Draw();
  */
  //c1->cd(1);

  p2->cd();

  Int_t nBins = dataALL->GetNbinsX();
  Double_t Min = dataALL->GetXaxis()->GetBinLowEdge(dataALL->GetXaxis()->GetFirst());
  Double_t Max = dataALL->GetXaxis()->GetBinUpEdge(dataALL->GetXaxis()->GetLast());
  
  TH1D *residALL = new TH1D("residALL","", nBins, Min, Max);
  residALL->Add(simALL,dataALL,1,-1);
  residALL->SetFillColor(1);
  residALL->GetXaxis()->SetRangeUser(0., xAxisMax);
  residALL->GetXaxis()->SetTitleOffset(1.);
  residALL->GetYaxis()->SetTitleOffset(0.5);
  residALL->GetYaxis()->SetTitle("event rate (mHz/keV)");
  //residALL->GetXaxis()->SetTitle("Energy (keV)");
  //residALL->SetMaximum(2.);
  //residALL->SetMinimum(-2.);
  residALL->SetLineWidth(2);
  residALL->SetMarkerColor(color?kBlue:kBlack);
  residALL->SetMarkerStyle(0);
  residALL->SetMarkerSize(1.);
  residALL->GetYaxis()->SetTitle("MC-Data");
  residALL->GetYaxis()->SetLabelSize(0.12);
  residALL->GetYaxis()->SetTitleOffset(0.4);
  residALL->GetYaxis()->SetTitleSize(0.14);
  residALL->GetXaxis()->SetLabelSize(0.12);
  residALL->GetXaxis()->CenterTitle();
  residALL->GetYaxis()->CenterTitle();
  residALL->GetYaxis()->SetNdivisions(508);

  

  residALL->Draw("E3");
  TLine *l = new TLine(Min, 0., xAxisMax, 0.);
  l->SetLineStyle(2);
  l->Draw();

  p1->cd();
  p1->SetBottomMargin(0.06);
  
  dataALL->SetTitle("");
  dataALL->SetMarkerColor(color?kBlue:kBlack);
  dataALL->SetMarkerStyle(24);
  dataALL->SetMarkerSize(0.75);
  dataALL->SetLineColor(color?kBlue:kBlack); 
  //dataALL->SetFillColor(kBlue);
  dataALL->SetLineWidth(3);
  dataALL->GetXaxis()->SetRangeUser(0., xAxisMax);
  simALL->SetTitle("");
  simALL->SetLineWidth(2);
  simALL->SetMarkerColor(kBlack);
  simALL->SetLineColor(kBlack);
  simALL->SetMarkerSize(0);
  simALL->SetMarkerStyle(24);
  simALL->SetMaximum(dataALL->GetMaximum()*1.2);
  simALL->GetYaxis()->SetTitle("Event Rate (mHz/keV)");
  simALL->GetYaxis()->SetLabelSize(0.06);
  simALL->GetXaxis()->SetLabelSize(0.055);
  simALL->GetYaxis()->SetTitleOffset(0.7);
  simALL->GetYaxis()->SetTitleSize(0.07);
  simALL->GetXaxis()->SetTitle("");
  simALL->GetXaxis()->SetRangeUser(0., xAxisMax);
  simALL->Draw("HIST C");
  dataALL->Draw("SAMEE0");

  BG_dataALL->SetMarkerColor(1);
  BG_dataALL->SetMarkerStyle(20);
  BG_dataALL->SetMarkerSize(0.75);
  BG_dataALL->SetLineColor(1);
  BG_dataALL->SetLineWidth(3);
  BG_dataALL->Draw("SAMEE0");

  TLegend *leg = new TLegend(0.61,0.60,0.88,0.85);
  leg->AddEntry(dataALL,"Data","p");
  leg->AddEntry(simALL,"Monte Carlo","l");
  leg->AddEntry(BG_dataALL,"Background","p");
  leg->Draw("SAME");

  c2->Update();

  c2->cd();
  TPaveText *pv_xTitle = new TPaveText(0.3,0.01,0.7,0.04,"nbNDC");
  pv_xTitle->SetBorderSize(0);
  pv_xTitle->AddText("E_{recon} (keV)");
  pv_xTitle->GetLine(0)->SetTextSize(0.040);
  pv_xTitle->GetLine(0)->SetTextFont(42);
  pv_xTitle->Draw();
  
  c2->Print(TString::Format("AsymmetryVsEnergy%s.pdf",(color?"_color":"")));
  
}

  

 
