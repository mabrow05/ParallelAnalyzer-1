{
  cout.setf(ios::fixed, ios::floatfield);
  cout.precision(12);

  // Style options
  gROOT->SetStyle("Plain");
  //gStyle->SetOptStat(11);
  gStyle->SetOptStat(0);
  gStyle->SetStatFontSize(0.020);
  gStyle->SetOptFit(0); // 1111
  gStyle->SetOptTitle(1);
  gStyle->SetTitleFontSize(0.05);
  //gStyle->SetTitleX(0.17);
  //gStyle->SetTitleAlign(13);
  gStyle->SetTitleOffset(1.20, "x");
  gStyle->SetTitleOffset(1.10, "y");
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetNdivisions(510,"X");
  gStyle->SetNdivisions(510,"Y");
  gStyle->SetPadLeftMargin(0.13);
  gStyle->SetPadRightMargin(0.05);
  gStyle->SetPadBottomMargin(0.12);

  //Run Range for this envelope
  //Int_t runLow = 17359;
  //Int_t runHigh = 19959;
  Int_t calibrationPeriod = 8;
  // Source peaks from simulation
  Double_t peakCe = 98.2;
  Double_t peakSn = 331.2;
  Double_t peakBiLow = 443.0;
  Double_t peakBiHigh = 928.0;

  // Read data file
  char temp[500];
  sprintf(temp, "../residuals/source_runs_RunPeriod_%i.dat",calibrationPeriod);
  ifstream filein(temp);

  Int_t i = 0;
  Double_t run[500], EQ[500];
  Double_t ADCE1[500], ADCE2[500], ADCE3[500], ADCE4[500];
  Double_t ADCW1[500], ADCW2[500], ADCW3[500], ADCW4[500];
  Double_t resE1[500], resE2[500], resE3[500], resE4[500];
  Double_t resW1[500], resW2[500], resW3[500], resW4[500];
  Double_t err[500];
  while (!filein.eof()) {
    filein >> run[i] >> EQ[i]
           >> ADCE1[i] >> ADCE2[i] >> ADCE3[i] >> ADCE4[i]
           >> ADCW1[i] >> ADCW2[i] >> ADCW3[i] >> ADCW4[i];
    if (filein.fail()) break;
    i++;
  }
  Int_t num = i;
  cout << "Number of data points: " << num << endl;

  //TFile *outfile = new TFile("EnergyCal_18745-18756.root","RECREATE");

  // Fit function
  TF1 *fitADC = new TF1("fitADC", "([0] + [1]*x + [2]*x*x)", 0.0, 2500.0);
  fitADC->SetParameter(0, 0.0);
  fitADC->SetParameter(1, 1.0);
  fitADC->SetParameter(2, 0.0);

  fitADC->SetParLimits(2, -0.0001, 0.0001);

  fitADC->SetNpx(100000);
  fitADC->SetLineColor(2);

  // East 1
  c1 = new TCanvas("c1", "c1");
  c1->SetLogy(0);

  TGraphErrors *grE1 = new TGraphErrors(num,ADCE1,EQ,err,err);
  grE1->SetTitle("");
  grE1->SetMarkerColor(1);
  grE1->SetLineColor(1);
  grE1->SetMarkerStyle(20);
  grE1->SetMarkerSize(0.75);
  grE1->GetXaxis()->SetTitle("g(t)#upoint#eta(x,y)#upointADC [channels]");
  grE1->GetXaxis()->SetTitleOffset(1.2);
  grE1->GetXaxis()->CenterTitle();
  grE1->GetYaxis()->SetTitle("E_{Q} [keV]");
  grE1->GetYaxis()->SetTitleOffset(1.6);
  grE1->GetYaxis()->CenterTitle();
  grE1->GetXaxis()->SetLimits(0.0,2400.0);
  grE1->SetMinimum(0.0);
  grE1->SetMaximum(1000.0);
  grE1->Draw("AP");

  grE1->Fit("fitADC", "R");
  Double_t offsetE1, slopeE1, quadE1;
  offsetE1 = fitADC->GetParameter(0);
  slopeE1 = fitADC->GetParameter(1);
  quadE1 = fitADC->GetParameter(2);

  Double_t x1_text = 1200;
  Double_t y1_text = 100;

  TPaveText *pt1 = new TPaveText(x1_text,y1_text,x1_text,y1_text,"");
  pt1->SetTextSize(0.042);
  pt1->SetTextColor(1);
  pt1->SetTextAlign(12);
  pt1->AddText("East PMT 1");
  pt1->SetBorderSize(0);
  pt1->SetFillStyle(0);
  pt1->SetFillColor(0);
  pt1->Draw();

  sprintf(temp,"../residuals/residuals_East_runPeriod_%i_PMTE1.dat",calibrationPeriod);
  ofstream oFileE1(temp);

  // Calculate residuals in [keV]
  Double_t fitEQ_E1[num];
  for (int j=0; j<num; j++) {
    fitEQ_E1[j]    = offsetE1 + slopeE1*ADCE1[j] + quadE1*ADCE1[j]*ADCE1[j];
    if (fitEQ_E1[j] < 200.) {
      resE1[j] = fitEQ_E1[j] - peakCe;
      oFileE1 << "Ce_East" << " " << (int) run[j] << " " << resE1[j] << endl;
      //resE1[j] = (fitEQ - peakCe)/peakCe * 100.;
    }
    else if (fitEQ_E1[j] < 400.) {
      resE1[j] = fitEQ_E1[j] - peakSn;
      //resE1[j] = (fitEQ - peakSn)/peakSn * 100.;
      oFileE1 << "Sn_East" << " " << (int) run[j] << " " << resE1[j] << endl;
    }
    else {
      resE1[j] = fitEQ_E1[j] - peakBiHigh;
      //resE1[j] = (fitEQ - peakBiHigh)/peakBiHigh * 100.;
      oFileE1 << "Bi_East" << " " << (int) run[j] << " " << resE1[j] << endl;
    }
  }

  // East 1 residuals
  c1r = new TCanvas("c1r", "c1r");
  c1r->SetLogy(0);

  TGraphErrors *grE1r = new TGraphErrors(num,ADCE1,resE1,err,err);
  grE1r->SetTitle("");
  grE1r->SetMarkerColor(1);
  grE1r->SetLineColor(1);
  grE1r->SetMarkerStyle(20);
  grE1r->SetMarkerSize(0.75);
  grE1r->GetXaxis()->SetTitle("g(t)#upoint#eta(x,y)#upointADC [channels]");
  grE1r->GetXaxis()->SetTitleOffset(1.2);
  grE1r->GetXaxis()->CenterTitle();
  grE1r->GetYaxis()->SetTitle("Residuals [keV]");
  grE1r->GetYaxis()->SetTitleOffset(1.2);
  grE1r->GetYaxis()->CenterTitle();
  grE1r->GetXaxis()->SetLimits(0.0,2400.0);
  grE1r->SetMinimum(-50.0);
  grE1r->SetMaximum( 50.0);
  grE1r->Draw("AP");

  Int_t n = 2;
  Double_t x[n] = {0, 2400};
  Double_t y[n] = {0.0, 0.0};

  TGraph *gr1 = new TGraph(n,x,y);
  gr1->Draw("Same");
  gr1->SetLineWidth(2);
  gr1->SetLineColor(1);
  gr1->SetLineStyle(2);

  Double_t x1_text = 1200;
  Double_t y1_text = -40;

  TPaveText *pt1 = new TPaveText(x1_text,y1_text,x1_text,y1_text,"");
  pt1->SetTextSize(0.042);
  pt1->SetTextColor(1);
  pt1->SetTextAlign(12);
  pt1->AddText("East PMT 1");
  pt1->SetBorderSize(0);
  pt1->SetFillStyle(0);
  pt1->SetFillColor(0);
  pt1->Draw();

  // East 2
  c2 = new TCanvas("c2", "c2");
  c2->SetLogy(0);

  TGraphErrors *grE2 = new TGraphErrors(num,ADCE2,EQ,err,err);
  grE2->SetTitle("");
  grE2->SetMarkerColor(1);
  grE2->SetLineColor(1);
  grE2->SetMarkerStyle(20);
  grE2->SetMarkerSize(0.75);
  grE2->GetXaxis()->SetTitle("g(t)#upoint#eta(x,y)#upointADC [channels]");
  grE2->GetXaxis()->SetTitleOffset(1.2);
  grE2->GetXaxis()->CenterTitle();
  grE2->GetYaxis()->SetTitle("E_{Q} [keV]");
  grE2->GetYaxis()->SetTitleOffset(1.6);
  grE2->GetYaxis()->CenterTitle();
  grE2->GetXaxis()->SetLimits(0.0,2400.0);
  grE2->SetMinimum(0.0);
  grE2->SetMaximum(1000.0);
  grE2->Draw("AP");

  grE2->Fit("fitADC", "R");
  Double_t offsetE2, slopeE2, quadE2;
  offsetE2 = fitADC->GetParameter(0);
  slopeE2 = fitADC->GetParameter(1);
  quadE2 = fitADC->GetParameter(2);

  Double_t x1_text = 1200;
  Double_t y1_text = 100;

  TPaveText *pt1 = new TPaveText(x1_text,y1_text,x1_text,y1_text,"");
  pt1->SetTextSize(0.042);
  pt1->SetTextColor(1);
  pt1->SetTextAlign(12);
  pt1->AddText("East PMT 2");
  pt1->SetBorderSize(0);
  pt1->SetFillStyle(0);
  pt1->SetFillColor(0);
  pt1->Draw();

  sprintf(temp,"../residuals/residuals_East_runPeriod_%i_PMTE2.dat",calibrationPeriod);
  ofstream oFileE2(temp);

  // Calculate residuals in [keV]
  Double_t fitEQ_E2[num];
  for (int j=0; j<num; j++) {
    fitEQ_E2[j]    = offsetE2 + slopeE2*ADCE2[j] + quadE2*ADCE2[j]*ADCE2[j];
    if (fitEQ_E2[j] < 200.) {
      resE2[j] = fitEQ_E2[j] - peakCe;
      //resE2[j] = (fitEQ - peakCe)/peakCe * 100.;
      oFileE2 << "Ce_East" << " " << (int) run[j] << " " << resE2[j] << endl;
    }
    else if (fitEQ_E2[j] < 400.) {
      resE2[j] = fitEQ_E2[j] - peakSn;
      //resE2[j] = (fitEQ - peakSn)/peakSn * 100.;
      oFileE2 << "Sn_East" << " " << (int) run[j] << " " << resE2[j] << endl;
    }
    else if (fitEQ_E2[j] < 600.) {
      resE2[j] = fitEQ_E2[j] - peakBiLow;
      //resE2[j] = (fitEQ - peakBiLow)/peakBiLow * 100.;
    }
    else {
      resE2[j] = fitEQ_E2[j] - peakBiHigh;
      //resE2[j] = (fitEQ - peakBiHigh)/peakBiHigh * 100.;
      oFileE2 << "Bi_East" << " " << (int) run[j] << " " << resE2[j] << endl;
    }
  }

  // East 2 residuals
  c2r = new TCanvas("c2r", "c2r");
  c2r->SetLogy(0);

  TGraphErrors *grE2r = new TGraphErrors(num,ADCE2,resE2,err,err);
  grE2r->SetTitle("");
  grE2r->SetMarkerColor(1);
  grE2r->SetLineColor(1);
  grE2r->SetMarkerStyle(20);
  grE2r->SetMarkerSize(0.75);
  grE2r->GetXaxis()->SetTitle("g(t)#upoint#eta(x,y)#upointADC [channels]");
  grE2r->GetXaxis()->SetTitleOffset(1.2);
  grE2r->GetXaxis()->CenterTitle();
  grE2r->GetYaxis()->SetTitle("Residuals [keV]");
  grE2r->GetYaxis()->SetTitleOffset(1.2);
  grE2r->GetYaxis()->CenterTitle();
  grE2r->GetXaxis()->SetLimits(0.0,2400.0);
  grE2r->SetMinimum(-50.0);
  grE2r->SetMaximum( 50.0);
  grE2r->Draw("AP");

  Int_t n = 2;
  Double_t x[n] = {0, 2400};
  Double_t y[n] = {0.0, 0.0};

  TGraph *gr1 = new TGraph(n,x,y);
  gr1->Draw("Same");
  gr1->SetLineWidth(2);
  gr1->SetLineColor(1);
  gr1->SetLineStyle(2);

  Double_t x1_text = 1200;
  Double_t y1_text = -40;

  TPaveText *pt1 = new TPaveText(x1_text,y1_text,x1_text,y1_text,"");
  pt1->SetTextSize(0.042);
  pt1->SetTextColor(1);
  pt1->SetTextAlign(12);
  pt1->AddText("East PMT 2");
  pt1->SetBorderSize(0);
  pt1->SetFillStyle(0);
  pt1->SetFillColor(0);
  pt1->Draw();

  // East 3
  c3 = new TCanvas("c3", "c3");
  c3->SetLogy(0);

  TGraphErrors *grE3 = new TGraphErrors(num,ADCE3,EQ,err,err);
  grE3->SetTitle("");
  grE3->SetMarkerColor(1);
  grE3->SetLineColor(1);
  grE3->SetMarkerStyle(20);
  grE3->SetMarkerSize(0.75);
  grE3->GetXaxis()->SetTitle("g(t)#upoint#eta(x,y)#upointADC [channels]");
  grE3->GetXaxis()->SetTitleOffset(1.2);
  grE3->GetXaxis()->CenterTitle();
  grE3->GetYaxis()->SetTitle("E_{Q} [keV]");
  grE3->GetYaxis()->SetTitleOffset(1.6);
  grE3->GetYaxis()->CenterTitle();
  grE3->GetXaxis()->SetLimits(0.0,2000.0);
  grE3->SetMinimum(0.0);
  grE3->SetMaximum(1000.0);
  grE3->Draw("AP");

  grE3->Fit("fitADC", "R");
  Double_t offsetE3, slopeE3, quadE3;
  offsetE3 = fitADC->GetParameter(0);
  slopeE3 = fitADC->GetParameter(1);
  quadE3 = fitADC->GetParameter(2);

  Double_t x1_text = 1200;
  Double_t y1_text = 100;

  TPaveText *pt1 = new TPaveText(x1_text,y1_text,x1_text,y1_text,"");
  pt1->SetTextSize(0.042);
  pt1->SetTextColor(1);
  pt1->SetTextAlign(12);
  pt1->AddText("East PMT 3");
  pt1->SetBorderSize(0);
  pt1->SetFillStyle(0);
  pt1->SetFillColor(0);
  pt1->Draw();

  sprintf(temp,"../residuals/residuals_East_runPeriod_%i_PMTE3.dat",calibrationPeriod);
  ofstream oFileE3(temp);

  // Calculate residuals in [keV]
  Double_t fitEQ_E3[num];
  for (int j=0; j<num; j++) {
    fitEQ_E3[j]    = offsetE3 + slopeE3*ADCE3[j] + quadE3*ADCE3[j]*ADCE3[j];
    if (fitEQ_E3[j] < 200.) {
      resE3[j] = fitEQ_E3[j] - peakCe;
      //resE3[j] = (fitEQ - peakCe)/peakCe * 100.;
      oFileE3 << "Ce_East" << " " << (int) run[j] << " " << resE3[j] << endl;
    }
    else if (fitEQ_E3[j] < 400.) {
      resE3[j] = fitEQ_E3[j] - peakSn;
      //resE3[j] = (fitEQ - peakSn)/peakSn * 100.;
      oFileE3 << "Sn_East" << " " << (int) run[j] << " " << resE3[j] << endl;
    }
    else if (fitEQ_E3[j] < 600.) {
      resE3[j] = fitEQ_E3[j] - peakBiLow;
      //resE3[j] = (fitEQ - peakBiLow)/peakBiLow * 100.;
    }
    else {
      resE3[j] = fitEQ_E3[j] - peakBiHigh;
      //resE3[j] = (fitEQ - peakBiHigh)/peakBiHigh * 100.;
      oFileE3 << "Bi_East" << " " << (int) run[j] << " " << resE3[j] << endl;
    }
  }

  // East 3 residuals
  c3r = new TCanvas("c3r", "c3r");
  c3r->SetLogy(0);

  TGraphErrors *grE3r = new TGraphErrors(num,ADCE3,resE3,err,err);
  grE3r->SetTitle("");
  grE3r->SetMarkerColor(1);
  grE3r->SetLineColor(1);
  grE3r->SetMarkerStyle(20);
  grE3r->SetMarkerSize(0.75);
  grE3r->GetXaxis()->SetTitle("g(t)#upoint#eta(x,y)#upointADC [channels]");
  grE3r->GetXaxis()->SetTitleOffset(1.2);
  grE3r->GetXaxis()->CenterTitle();
  grE3r->GetYaxis()->SetTitle("Residuals [keV]");
  grE3r->GetYaxis()->SetTitleOffset(1.2);
  grE3r->GetYaxis()->CenterTitle();
  grE3r->GetXaxis()->SetLimits(0.0,2400.0);
  grE3r->SetMinimum(-50.0);
  grE3r->SetMaximum( 50.0);
  grE3r->Draw("AP");

  Int_t n = 2;
  Double_t x[n] = {0, 2400};
  Double_t y[n] = {0.0, 0.0};

  TGraph *gr1 = new TGraph(n,x,y);
  gr1->Draw("Same");
  gr1->SetLineWidth(2);
  gr1->SetLineColor(1);
  gr1->SetLineStyle(2);

  Double_t x1_text = 1200;
  Double_t y1_text = -40;

  TPaveText *pt1 = new TPaveText(x1_text,y1_text,x1_text,y1_text,"");
  pt1->SetTextSize(0.042);
  pt1->SetTextColor(1);
  pt1->SetTextAlign(12);
  pt1->AddText("East PMT 3");
  pt1->SetBorderSize(0);
  pt1->SetFillStyle(0);
  pt1->SetFillColor(0);
  pt1->Draw();

  // East 4
  c4 = new TCanvas("c4", "c4");
  c4->SetLogy(0);

  TGraphErrors *grE4 = new TGraphErrors(num,ADCE4,EQ,err,err);
  grE4->SetTitle("");
  grE4->SetMarkerColor(1);
  grE4->SetLineColor(1);
  grE4->SetMarkerStyle(20);
  grE4->SetMarkerSize(0.75);
  grE4->GetXaxis()->SetTitle("g(t)#upoint#eta(x,y)#upointADC [channels]");
  grE4->GetXaxis()->SetTitleOffset(1.2);
  grE4->GetXaxis()->CenterTitle();
  grE4->GetYaxis()->SetTitle("E_{Q} [keV]");
  grE4->GetYaxis()->SetTitleOffset(1.6);
  grE4->GetYaxis()->CenterTitle();
  grE4->GetXaxis()->SetLimits(0.0,2400.0);
  grE4->SetMinimum(0.0);
  grE4->SetMaximum(1000.0);
  grE4->Draw("AP");

  grE4->Fit("fitADC", "R");
  Double_t offsetE4, slopeE4, quadE4;
  offsetE4 = fitADC->GetParameter(0);
  slopeE4 = fitADC->GetParameter(1);
  quadE4 = fitADC->GetParameter(2);

  Double_t x1_text = 1200;
  Double_t y1_text = 100;

  TPaveText *pt1 = new TPaveText(x1_text,y1_text,x1_text,y1_text,"");
  pt1->SetTextSize(0.042);
  pt1->SetTextColor(1);
  pt1->SetTextAlign(12);
  pt1->AddText("East PMT 4");
  pt1->SetBorderSize(0);
  pt1->SetFillStyle(0);
  pt1->SetFillColor(0);
  pt1->Draw();

  sprintf(temp,"../residuals/residuals_East_runPeriod_%i_PMTE4.dat",calibrationPeriod);
  ofstream oFileE4(temp);

  // Calculate residuals in [keV]
  Double_t fitEQ_E4[num];
  for (int j=0; j<num; j++) {
    fitEQ_E4[j]    = offsetE4 + slopeE4*ADCE4[j] + quadE4*ADCE4[j]*ADCE4[j];
    if (fitEQ_E4[j] < 200.) {
      resE4[j] = fitEQ_E4[j] - peakCe;
      //resE4[j] = (fitEQ - peakCe)/peakCe * 100.;
      oFileE4 << "Ce_East" << " " << (int) run[j] << " " << resE4[j] << endl;
    }
    else if (fitEQ_E4[j] < 400.) {
      resE4[j] = fitEQ_E4[j] - peakSn;
      //resE4[j] = (fitEQ - peakSn)/peakSn * 100.;
      oFileE4 << "Sn_East" << " " << (int) run[j] << " " << resE4[j] << endl;
    }
    else if (fitEQ_E4[j] < 600.) {
      resE4[j] = fitEQ_E4[j] - peakBiLow;
      //resE4[j] = (fitEQ - peakBiLow)/peakBiLow * 100.;
    }
    else {
      resE4[j] = fitEQ_E4[j] - peakBiHigh;
      //resE4[j] = (fitEQ - peakBiHigh)/peakBiHigh * 100.;
      oFileE4 << "Bi_East" << " " << (int) run[j] << " " << resE4[j] << endl;
    }
  }

  // East 4 residuals
  c4r = new TCanvas("c4r", "c4r");
  c4r->SetLogy(0);

  TGraphErrors *grE4r = new TGraphErrors(num,ADCE4,resE4,err,err);
  grE4r->SetTitle("");
  grE4r->SetMarkerColor(1);
  grE4r->SetLineColor(1);
  grE4r->SetMarkerStyle(20);
  grE4r->SetMarkerSize(0.75);
  grE4r->GetXaxis()->SetTitle("g(t)#upoint#eta(x,y)#upointADC [channels]");
  grE4r->GetXaxis()->SetTitleOffset(1.2);
  grE4r->GetXaxis()->CenterTitle();
  grE4r->GetYaxis()->SetTitle("Residuals [keV]");
  grE4r->GetYaxis()->SetTitleOffset(1.2);
  grE4r->GetYaxis()->CenterTitle();
  grE4r->GetXaxis()->SetLimits(0.0,2400.0);
  grE4r->SetMinimum(-50.0);
  grE4r->SetMaximum( 50.0);
  grE4r->Draw("AP");

  Int_t n = 2;
  Double_t x[n] = {0, 2400};
  Double_t y[n] = {0.0, 0.0};

  TGraph *gr1 = new TGraph(n,x,y);
  gr1->Draw("Same");
  gr1->SetLineWidth(2);
  gr1->SetLineColor(1);
  gr1->SetLineStyle(2);

  Double_t x1_text = 1200;
  Double_t y1_text = -40;

  TPaveText *pt1 = new TPaveText(x1_text,y1_text,x1_text,y1_text,"");
  pt1->SetTextSize(0.042);
  pt1->SetTextColor(1);
  pt1->SetTextAlign(12);
  pt1->AddText("East PMT 4");
  pt1->SetBorderSize(0);
  pt1->SetFillStyle(0);
  pt1->SetFillColor(0);
  pt1->Draw();

  sprintf(temp,"../residuals/residuals_East_runPeriod_%i.dat",calibrationPeriod);
  ofstream oFileE(temp);

  // For now, simply average the East PMT energies
  cout << "CALCULATING EAST RESIDUALS" << endl;
  double EQ_East[num], res_East[num], x_East[num];
  for (int j=0; j<num; j++) {
    EQ_East[j] = 0.25*(fitEQ_E1[j] + fitEQ_E2[j] + fitEQ_E3[j] + fitEQ_E4[j]);
    //EQ_East[j] = 0.5*(fitEQ_E1[j] + fitEQ_E2[j]);

    if (EQ_East[j] < 200.) {
      res_East[j] = EQ_East[j] - peakCe;
      x_East[j] = peakCe;
      oFileE << "Ce_East" << " " << (int) run[j] << " " << res_East[j] << endl;
    }
    else if (EQ_East[j] < 400.) {
      res_East[j] = EQ_East[j] - peakSn;
      x_East[j] = peakSn;
      oFileE << "Sn_East" << " " << (int) run[j] << " " << res_East[j] << endl;
    }
    else if (EQ_East[j] < 600.) {
      res_East[j] = EQ_East[j] - peakBiLow;
      x_East[j] = peakBiLow;
    }
    else {
      res_East[j] = EQ_East[j] - peakBiHigh;
      x_East[j] = peakBiHigh;
      oFileE << "Bi_East" << " " << (int) run[j] << " " << res_East[j] << endl;
   }

  }
  oFileE.close();

  // East Average residuals
  cEr = new TCanvas("cEr", "cEr");
  cEr->SetLogy(0);

  TGraphErrors *grEr = new TGraphErrors(num,x_East,res_East,err,err);
  grEr->SetTitle("");
  grEr->SetMarkerColor(1);
  grEr->SetLineColor(1);
  grEr->SetMarkerStyle(20);
  grEr->SetMarkerSize(0.75);
  grEr->GetXaxis()->SetTitle("E_{Q} [keV]");
  grEr->GetXaxis()->SetTitleOffset(1.2);
  grEr->GetXaxis()->CenterTitle();
  grEr->GetYaxis()->SetTitle("Residuals [keV]");
  grEr->GetYaxis()->SetTitleOffset(1.2);
  grEr->GetYaxis()->CenterTitle();
  grEr->GetXaxis()->SetLimits(0.0,1200.0);
  grEr->SetMinimum(-30.0);
  grEr->SetMaximum( 30.0);
  grEr->Draw("AP");

  Int_t n = 2;
  Double_t x[n] = {0, 2400};
  Double_t y[n] = {0.0, 0.0};

  TGraph *gr1 = new TGraph(n,x,y);
  gr1->Draw("Same");
  gr1->SetLineWidth(2);
  gr1->SetLineColor(1);
  gr1->SetLineStyle(2);

  Double_t x1_text = 1000;
  Double_t y1_text = -20;

  TPaveText *pt1 = new TPaveText(x1_text,y1_text,x1_text,y1_text,"");
  pt1->SetTextSize(0.042);
  pt1->SetTextColor(1);
  pt1->SetTextAlign(12);
  pt1->AddText("East");
  pt1->SetBorderSize(0);
  pt1->SetFillStyle(0);
  pt1->SetFillColor(0);
  pt1->Draw();

  // West 1
  cW1 = new TCanvas("cW1", "cW1");
  cW1->SetLogy(0);

  TGraphErrors *grW1 = new TGraphErrors(num,ADCW1,EQ,err,err);
  grW1->SetTitle("");
  grW1->SetMarkerColor(1);
  grW1->SetLineColor(1);
  grW1->SetMarkerStyle(20);
  grW1->SetMarkerSize(0.75);
  grW1->GetXaxis()->SetTitle("g(t)#upoint#eta(x,y)#upointADC [channels]");
  grW1->GetXaxis()->SetTitleOffset(1.2);
  grW1->GetXaxis()->CenterTitle();
  grW1->GetYaxis()->SetTitle("E_{Q} [keV]");
  grW1->GetYaxis()->SetTitleOffset(1.6);
  grW1->GetYaxis()->CenterTitle();
  grW1->GetXaxis()->SetLimits(0.0,2400.0);
  grW1->SetMinimum(0.0);
  grW1->SetMaximum(1000.0);
  grW1->Draw("AP");

  grW1->Fit("fitADC", "R");
  Double_t offsetW1, slopeW1, quadW1;
  offsetW1 = fitADC->GetParameter(0);
  slopeW1 = fitADC->GetParameter(1);
  quadW1 = fitADC->GetParameter(2);

  Double_t x1_text = 1200;
  Double_t y1_text = 100;

  TPaveText *pt1 = new TPaveText(x1_text,y1_text,x1_text,y1_text,"");
  pt1->SetTextSize(0.042);
  pt1->SetTextColor(1);
  pt1->SetTextAlign(12);
  pt1->AddText("West PMT 1");
  pt1->SetBorderSize(0);
  pt1->SetFillStyle(0);
  pt1->SetFillColor(0);
  pt1->Draw();

  sprintf(temp,"../residuals/residuals_West_runPeriod_%i_PMTW1.dat",calibrationPeriod);
  ofstream oFileW1(temp);
  
  // Calculate residuals in [keV]
  Double_t fitEQ_W1[num];
  for (int j=0; j<num; j++) {
    fitEQ_W1[j]    = offsetW1 + slopeW1*ADCW1[j] + quadW1*ADCW1[j]*ADCW1[j];
    if (fitEQ_W1[j] < 200.) {
      resW1[j] = fitEQ_W1[j] - peakCe;
      oFileW1 << "Ce_West" << " " << (int) run[j] << " " << resW1[j] << endl;
    }
    else if (fitEQ_W1[j] < 400.) {
      resW1[j] = fitEQ_W1[j] - peakSn;
      oFileW1 << "Sn_West" << " " << (int) run[j] << " " << resW1[j] << endl;
    }
    else {
      resW1[j] = fitEQ_W1[j] - peakBiHigh;
      oFileW1 << "Bi_West" << " " << (int) run[j] << " " << resW1[j] << endl;
    }
  }

  // West 1 residuals
  cW1r = new TCanvas("cW1r", "cW1r");
  cW1r->SetLogy(0);

  TGraphErrors *grW1r = new TGraphErrors(num,ADCW1,resW1,err,err);
  grW1r->SetTitle("");
  grW1r->SetMarkerColor(1);
  grW1r->SetLineColor(1);
  grW1r->SetMarkerStyle(20);
  grW1r->SetMarkerSize(0.75);
  grW1r->GetXaxis()->SetTitle("g(t)#upoint#eta(x,y)#upointADC [channels]");
  grW1r->GetXaxis()->SetTitleOffset(1.2);
  grW1r->GetXaxis()->CenterTitle();
  grW1r->GetYaxis()->SetTitle("Residuals [keV]");
  grW1r->GetYaxis()->SetTitleOffset(1.2);
  grW1r->GetYaxis()->CenterTitle();
  grW1r->GetXaxis()->SetLimits(0.0,2400.0);
  grW1r->SetMinimum(-50.0);
  grW1r->SetMaximum( 50.0);
  grW1r->Draw("AP");

  Int_t n = 2;
  Double_t x[n] = {0, 2400};
  Double_t y[n] = {0.0, 0.0};

  TGraph *gr1 = new TGraph(n,x,y);
  gr1->Draw("Same");
  gr1->SetLineWidth(2);
  gr1->SetLineColor(1);
  gr1->SetLineStyle(2);

  Double_t x1_text = 1200;
  Double_t y1_text = -40;

  TPaveText *pt1 = new TPaveText(x1_text,y1_text,x1_text,y1_text,"");
  pt1->SetTextSize(0.042);
  pt1->SetTextColor(1);
  pt1->SetTextAlign(12);
  pt1->AddText("West PMT 1");
  pt1->SetBorderSize(0);
  pt1->SetFillStyle(0);
  pt1->SetFillColor(0);
  pt1->Draw();

  // West 2
  cW2 = new TCanvas("cW2", "cW2");
  cW2->SetLogy(0);

  TGraphErrors *grW2 = new TGraphErrors(num,ADCW2,EQ,err,err);
  grW2->SetTitle("");
  grW2->SetMarkerColor(1);
  grW2->SetLineColor(1);
  grW2->SetMarkerStyle(20);
  grW2->SetMarkerSize(0.75);
  grW2->GetXaxis()->SetTitle("g(t)#upoint#eta(x,y)#upointADC [channels]");
  grW2->GetXaxis()->SetTitleOffset(1.2);
  grW2->GetXaxis()->CenterTitle();
  grW2->GetYaxis()->SetTitle("E_{Q} [keV]");
  grW2->GetYaxis()->SetTitleOffset(1.6);
  grW2->GetYaxis()->CenterTitle();
  grW2->GetXaxis()->SetLimits(0.0,2400.0);
  grW2->SetMinimum(0.0);
  grW2->SetMaximum(1000.0);
  grW2->Draw("AP");

  grW2->Fit("fitADC", "R");
  Double_t offsetW2, slopeW2, quadW2;
  offsetW2 = fitADC->GetParameter(0);
  slopeW2 = fitADC->GetParameter(1);
  quadW2 = fitADC->GetParameter(2);

  Double_t x1_text = 1200;
  Double_t y1_text = 100;

  TPaveText *pt1 = new TPaveText(x1_text,y1_text,x1_text,y1_text,"");
  pt1->SetTextSize(0.042);
  pt1->SetTextColor(1);
  pt1->SetTextAlign(12);
  pt1->AddText("West PMT 2");
  pt1->SetBorderSize(0);
  pt1->SetFillStyle(0);
  pt1->SetFillColor(0);
  pt1->Draw();

  sprintf(temp,"../residuals/residuals_West_runPeriod_%i_PMTW2.dat",calibrationPeriod);
  ofstream oFileW2(temp);

  // Calculate residuals in [keV]
  Double_t fitEQ_W2[num];
  for (int j=0; j<num; j++) {
    fitEQ_W2[j]    = offsetW2 + slopeW2*ADCW2[j] + quadW2*ADCW2[j]*ADCW2[j];
    if (fitEQ_W2[j] < 200.) {
      resW2[j] = fitEQ_W2[j] - peakCe;
      oFileW2 << "Ce_West" << " " << (int) run[j] << " " << resW2[j] << endl;
    }
    else if (fitEQ_W2[j] < 400.) {
      resW2[j] = fitEQ_W2[j] - peakSn;
      oFileW2 << "Sn_West" << " " << (int) run[j] << " " << resW2[j] << endl;
    }
    else {
      resW2[j] = fitEQ_W2[j] - peakBiHigh;
      oFileW2 << "Bi_West" << " " << (int) run[j] << " " << resW2[j] << endl;
    }
  }

  // West 2 residuals
  cW2r = new TCanvas("cW2r", "cW2r");
  cW2r->SetLogy(0);

  TGraphErrors *grW2r = new TGraphErrors(num,ADCW2,resW2,err,err);
  grW2r->SetTitle("");
  grW2r->SetMarkerColor(1);
  grW2r->SetLineColor(1);
  grW2r->SetMarkerStyle(20);
  grW2r->SetMarkerSize(0.75);
  grW2r->GetXaxis()->SetTitle("g(t)#upoint#eta(x,y)#upointADC [channels]");
  grW2r->GetXaxis()->SetTitleOffset(1.2);
  grW2r->GetXaxis()->CenterTitle();
  grW2r->GetYaxis()->SetTitle("Residuals [keV]");
  grW2r->GetYaxis()->SetTitleOffset(1.2);
  grW2r->GetYaxis()->CenterTitle();
  grW2r->GetXaxis()->SetLimits(0.0,2400.0);
  grW2r->SetMinimum(-50.0);
  grW2r->SetMaximum( 50.0);
  grW2r->Draw("AP");

  Int_t n = 2;
  Double_t x[n] = {0, 2400};
  Double_t y[n] = {0.0, 0.0};

  TGraph *gr1 = new TGraph(n,x,y);
  gr1->Draw("Same");
  gr1->SetLineWidth(2);
  gr1->SetLineColor(1);
  gr1->SetLineStyle(2);

  Double_t x1_text = 1200;
  Double_t y1_text = -40;

  TPaveText *pt1 = new TPaveText(x1_text,y1_text,x1_text,y1_text,"");
  pt1->SetTextSize(0.042);
  pt1->SetTextColor(1);
  pt1->SetTextAlign(12);
  pt1->AddText("West PMT 2");
  pt1->SetBorderSize(0);
  pt1->SetFillStyle(0);
  pt1->SetFillColor(0);
  pt1->Draw();

  // West 3
  cW3 = new TCanvas("cW3", "cW3");
  cW3->SetLogy(0);

  TGraphErrors *grW3 = new TGraphErrors(num,ADCW3,EQ,err,err);
  grW3->SetTitle("");
  grW3->SetMarkerColor(1);
  grW3->SetLineColor(1);
  grW3->SetMarkerStyle(20);
  grW3->SetMarkerSize(0.75);
  grW3->GetXaxis()->SetTitle("g(t)#upoint#eta(x,y)#upointADC [channels]");
  grW3->GetXaxis()->SetTitleOffset(1.2);
  grW3->GetXaxis()->CenterTitle();
  grW3->GetYaxis()->SetTitle("E_{Q} [keV]");
  grW3->GetYaxis()->SetTitleOffset(1.6);
  grW3->GetYaxis()->CenterTitle();
  grW3->GetXaxis()->SetLimits(0.0,2400.0);
  grW3->SetMinimum(0.0);
  grW3->SetMaximum(1000.0);
  grW3->Draw("AP");

  grW3->Fit("fitADC", "R");
  Double_t offsetW3, slopeW3, quadW3;
  offsetW3 = fitADC->GetParameter(0);
  slopeW3 = fitADC->GetParameter(1);
  quadW3 = fitADC->GetParameter(2);

  Double_t x1_text = 1200;
  Double_t y1_text = 100;

  TPaveText *pt1 = new TPaveText(x1_text,y1_text,x1_text,y1_text,"");
  pt1->SetTextSize(0.042);
  pt1->SetTextColor(1);
  pt1->SetTextAlign(12);
  pt1->AddText("West PMT 3");
  pt1->SetBorderSize(0);
  pt1->SetFillStyle(0);
  pt1->SetFillColor(0);
  pt1->Draw();

  sprintf(temp,"../residuals/residuals_West_runPeriod_%i_PMTW3.dat",calibrationPeriod);
  ofstream oFileW3(temp);

  // Calculate residuals in [keV]
  Double_t fitEQ_W3[num];
  for (int j=0; j<num; j++) {
    fitEQ_W3[j]    = offsetW3 + slopeW3*ADCW3[j] + quadW3*ADCW3[j]*ADCW3[j];
    if (fitEQ_W3[j] < 200.) {
      resW3[j] = fitEQ_W3[j] - peakCe;
      oFileW3 << "Ce_West" << " " << (int) run[j] << " " << resW3[j] << endl;
    }
    else if (fitEQ_W3[j] < 400.) {
      resW3[j] = fitEQ_W3[j] - peakSn;
      oFileW3 << "Sn_West" << " " << (int) run[j] << " " << resW3[j] << endl;
    }
    else {
      resW3[j] = fitEQ_W3[j] - peakBiHigh;
      oFileW3 << "Bi_West" << " " << (int) run[j] << " " << resW3[j] << endl;
    }
  }

  // West 3 residuals
  cW3r = new TCanvas("cW3r", "cW3r");
  cW3r->SetLogy(0);

  TGraphErrors *grW3r = new TGraphErrors(num,ADCW3,resW3,err,err);
  grW3r->SetTitle("");
  grW3r->SetMarkerColor(1);
  grW3r->SetLineColor(1);
  grW3r->SetMarkerStyle(20);
  grW3r->SetMarkerSize(0.75);
  grW3r->GetXaxis()->SetTitle("g(t)#upoint#eta(x,y)#upointADC [channels]");
  grW3r->GetXaxis()->SetTitleOffset(1.2);
  grW3r->GetXaxis()->CenterTitle();
  grW3r->GetYaxis()->SetTitle("Residuals [keV]");
  grW3r->GetYaxis()->SetTitleOffset(1.2);
  grW3r->GetYaxis()->CenterTitle();
  grW3r->GetXaxis()->SetLimits(0.0,2400.0);
  grW3r->SetMinimum(-50.0);
  grW3r->SetMaximum( 50.0);
  grW3r->Draw("AP");

  Int_t n = 2;
  Double_t x[n] = {0, 2400};
  Double_t y[n] = {0.0, 0.0};

  TGraph *gr1 = new TGraph(n,x,y);
  gr1->Draw("Same");
  gr1->SetLineWidth(2);
  gr1->SetLineColor(1);
  gr1->SetLineStyle(2);

  Double_t x1_text = 1200;
  Double_t y1_text = -40;

  TPaveText *pt1 = new TPaveText(x1_text,y1_text,x1_text,y1_text,"");
  pt1->SetTextSize(0.042);
  pt1->SetTextColor(1);
  pt1->SetTextAlign(12);
  pt1->AddText("West PMT 3");
  pt1->SetBorderSize(0);
  pt1->SetFillStyle(0);
  pt1->SetFillColor(0);
  pt1->Draw();

  // West 4
  cW4 = new TCanvas("cW4", "cW4");
  cW4->SetLogy(0);

  TGraphErrors *grW4 = new TGraphErrors(num,ADCW4,EQ,err,err);
  grW4->SetTitle("");
  grW4->SetMarkerColor(1);
  grW4->SetLineColor(1);
  grW4->SetMarkerStyle(20);
  grW4->SetMarkerSize(0.75);
  grW4->GetXaxis()->SetTitle("g(t)#upoint#eta(x,y)#upointADC [channels]");
  grW4->GetXaxis()->SetTitleOffset(1.2);
  grW4->GetXaxis()->CenterTitle();
  grW4->GetYaxis()->SetTitle("E_{Q} [keV]");
  grW4->GetYaxis()->SetTitleOffset(1.6);
  grW4->GetYaxis()->CenterTitle();
  grW4->GetXaxis()->SetLimits(0.0,4096.0);
  grW4->SetMinimum(0.0);
  grW4->SetMaximum(1000.0);
  grW4->Draw("AP");

  fitADC->SetRange(0., 3500.);
  grW4->Fit("fitADC", "R");
  Double_t offsetW4, slopeW4, quadW4;
  offsetW4 = fitADC->GetParameter(0);
  slopeW4 = fitADC->GetParameter(1);
  quadW4 = fitADC->GetParameter(2);

  Double_t x1_text = 1200;
  Double_t y1_text = 100;

  TPaveText *pt1 = new TPaveText(x1_text,y1_text,x1_text,y1_text,"");
  pt1->SetTextSize(0.042);
  pt1->SetTextColor(1);
  pt1->SetTextAlign(12);
  pt1->AddText("West PMT 4");
  pt1->SetBorderSize(0);
  pt1->SetFillStyle(0);
  pt1->SetFillColor(0);
  pt1->Draw();

  sprintf(temp,"../residuals/residuals_West_runPeriod_%i_PMTW4.dat",calibrationPeriod);
  ofstream oFileW4(temp);

  // Calculate residuals in [keV]
  Double_t fitEQ_W4[num];
  for (int j=0; j<num; j++) {
    fitEQ_W4[j]    = offsetW4 + slopeW4*ADCW4[j] + quadW4[j]*ADCW4[j]*ADCW4[j];
    if (fitEQ_W4[j] < 200.) {
      resW4[j] = fitEQ_W4[j] - peakCe;
      oFileW4 << "Ce_West" << " " << (int) run[j] << " " << resW4[j] << endl;
    }
    else if (fitEQ_W4[j] < 400.) {
      resW4[j] = fitEQ_W4[j] - peakSn;
      oFileW4 << "Sn_West" << " " << (int) run[j] << " " << resW4[j] << endl;
    }
    else {
      resW4[j] = fitEQ_W4[j] - peakBiHigh;
      oFileW4 << "Bi_West" << " " << (int) run[j] << " " << resW4[j] << endl;
    }
  }

  // West 4 residuals
  cW4r = new TCanvas("cW4r", "cW4r");
  cW4r->SetLogy(0);

  TGraphErrors *grW4r = new TGraphErrors(num,ADCW4,resW4,err,err);
  grW4r->SetTitle("");
  grW4r->SetMarkerColor(1);
  grW4r->SetLineColor(1);
  grW4r->SetMarkerStyle(20);
  grW4r->SetMarkerSize(0.75);
  grW4r->GetXaxis()->SetTitle("g(t)#upoint#eta(x,y)#upointADC [channels]");
  grW4r->GetXaxis()->SetTitleOffset(1.2);
  grW4r->GetXaxis()->CenterTitle();
  grW4r->GetYaxis()->SetTitle("Residuals [keV]");
  grW4r->GetYaxis()->SetTitleOffset(1.2);
  grW4r->GetYaxis()->CenterTitle();
  grW4r->GetXaxis()->SetLimits(0.0,4000.);
  grW4r->SetMinimum(-50.0);
  grW4r->SetMaximum( 50.0);
  grW4r->Draw("AP");

  Int_t n = 2;
  Double_t x[n] = {0, 4000.};
  Double_t y[n] = {0.0, 0.0};

  TGraph *gr1 = new TGraph(n,x,y);
  gr1->Draw("Same");
  gr1->SetLineWidth(2);
  gr1->SetLineColor(1);
  gr1->SetLineStyle(2);

  Double_t x1_text = 1200;
  Double_t y1_text = -40;

  TPaveText *pt1 = new TPaveText(x1_text,y1_text,x1_text,y1_text,"");
  pt1->SetTextSize(0.042);
  pt1->SetTextColor(1);
  pt1->SetTextAlign(12);
  pt1->AddText("West PMT 4");
  pt1->SetBorderSize(0);
  pt1->SetFillStyle(0);
  pt1->SetFillColor(0);
  pt1->Draw();

  sprintf(temp,"../residuals/residuals_West_runPeriod_%i.dat",calibrationPeriod);
  ofstream oFileW(temp);

  // For now, simply average the West PMT energies
  cout << "CALCULATING WEST RESIDUALS" << endl;
  double EQ_West[num], res_West[num], x_West[num];
  for (int j=0; j<num; j++) {
    EQ_West[j] = 0.33*(fitEQ_W1[j] + fitEQ_W2[j] + fitEQ_W3[j]);

    if (EQ_West[j] < 200.) {
      res_West[j] = EQ_West[j] - peakCe;
      x_West[j] = peakCe;
      oFileW << "Ce_West" << " " << (int) run[j] << " " << res_West[j] << endl;
    }
    else if (EQ_West[j] < 400.) {
      res_West[j] = EQ_West[j] - peakSn;
      x_West[j] = peakSn;
      oFileW << "Sn_West" << " " << (int) run[j] << " " << res_West[j] << endl;
    }
    else if (EQ_West[j] < 600.) {
      res_West[j] = EQ_West[j] - peakBiLow;
      x_West[j] = peakBiLow;
    }
    else {
      res_West[j] = EQ_West[j] - peakBiHigh;
      x_West[j] = peakBiHigh;
      oFileW << "Bi_West" << " " << (int) run[j] << " " << res_West[j] << endl;
    }

  }
  oFileW.close();

  // West Average residuals
  cWr = new TCanvas("cWr", "cWr");
  cWr->SetLogy(0);

  TGraphErrors *grWr = new TGraphErrors(num,x_West,res_West,err,err);
  grWr->SetTitle("");
  grWr->SetMarkerColor(1);
  grWr->SetLineColor(1);
  grWr->SetMarkerStyle(20);
  grWr->SetMarkerSize(0.75);
  grWr->GetXaxis()->SetTitle("E_{Q} [keV]");
  grWr->GetXaxis()->SetTitleOffset(1.2);
  grWr->GetXaxis()->CenterTitle();
  grWr->GetYaxis()->SetTitle("Residuals [keV]");
  grWr->GetYaxis()->SetTitleOffset(1.2);
  grWr->GetYaxis()->CenterTitle();
  grWr->GetXaxis()->SetLimits(0.0,1200.0);
  grWr->SetMinimum(-30.0);
  grWr->SetMaximum( 30.0);
  grWr->Draw("AP");

  Int_t n = 2;
  Double_t x[n] = {0, 2400};
  Double_t y[n] = {0.0, 0.0};

  TGraph *gr1 = new TGraph(n,x,y);
  gr1->Draw("Same");
  gr1->SetLineWidth(2);
  gr1->SetLineColor(1);
  gr1->SetLineStyle(2);

  Double_t x1_text = 1000;
  Double_t y1_text = -20;

  TPaveText *pt1 = new TPaveText(x1_text,y1_text,x1_text,y1_text,"");
  pt1->SetTextSize(0.042);
  pt1->SetTextColor(1);
  pt1->SetTextAlign(12);
  pt1->AddText("West");
  pt1->SetBorderSize(0);
  pt1->SetFillStyle(0);
  pt1->SetFillColor(0);
  pt1->Draw();
  
  //outfile->Write();
  //outfile->Close();

}
