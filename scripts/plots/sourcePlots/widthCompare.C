

void widthCompare(Int_t srcPeriod) {


  std::vector <Double_t> dataEastWidth;
  std::vector <Double_t> simEastWidth;

  std::vector <Double_t> dataWestWidth;
  std::vector <Double_t> simWestWidth;


  std::ifstream infile1(TString::Format("%s/residuals/source_runs_Erecon_RunPeriod_%i.dat",
					getenv("ANALYSIS_CODE"),srcPeriod));
  std::ifstream infile2(TString::Format("%s/residuals/SIM_source_runs_Erecon_RunPeriod_%i.dat",
					getenv("ANALYSIS_CODE"),srcPeriod));


  std::string src;
  Int_t run;
  Double_t eastPeak,eastWidth,westPeak,westWidth;
  Double_t SIMeastPeak,SIMeastWidth,SIMwestPeak,SIMwestWidth;

  
  while (infile1 >> run >> src >> eastPeak >> eastWidth >> westPeak >> westWidth) {
    infile2 >> run >> src >> SIMeastPeak >> SIMeastWidth >> SIMwestPeak >> SIMwestWidth;
    //if ( src==std::string("Ce") || src==std::string("Sn") || src.substr(0,2)==std::string("Bi") ) {
    if (src==std::string("In")) continue;  
    std::cout << src << endl;
      dataEastWidth.push_back(eastWidth);
      dataWestWidth.push_back(westWidth);
      simEastWidth.push_back(SIMeastWidth);
      simWestWidth.push_back(SIMwestWidth);
      //}
  }
  std::cout << "made it\n";

  gStyle->SetOptFit(0011);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetTitleYOffset(1.5);

  TCanvas *c1 = new TCanvas("c1","c1",1200,800);
  c1->Divide(2,1);
  c1->cd(1);

  TGraph *gE = new TGraph(dataEastWidth.size(),&dataEastWidth[0],&simEastWidth[0]);
  gE->SetMarkerStyle(kOpenCircle);
  gE->SetMarkerColor(kBlue);
  gE->SetTitle("Source Widths East");
  gE->Draw("AP");
  gE->GetXaxis()->SetLimits(0.,120.);
  gE->SetMaximum(120.);
  gE->SetMinimum(0.);
  gE->GetXaxis()->SetTitle("Data Widths");
  gE->GetYaxis()->SetTitle("Simulated Widths");
  gE->Draw("AP");
  c1->Update();


  TF1 *f = new TF1("f","[0]*x",0.,120.);
  f->SetParameter(0,1.);
  gE->Fit(f);
  gPad->Modified();


  c1->cd(2);

  TGraph *gW = new TGraph(dataWestWidth.size(),&dataWestWidth[0],&simWestWidth[0]);
  gW->SetMarkerStyle(kOpenCircle);
  gW->SetMarkerColor(kBlue);
  gW->SetTitle("Source Widths West");
  gW->Draw("AP");
  gW->GetXaxis()->SetLimits(0.,120.);
  gW->SetMaximum(120.);
  gW->SetMinimum(0.);
  gW->GetXaxis()->SetTitle("Data Widths");
  gW->GetYaxis()->SetTitle("Simulated Widths");
  gW->Draw("AP");
  c1->Update();



  gW->Fit(f);
  gPad->Modified();

  c1->Print(TString::Format("widthsErecon_runPeriod_%i.pdf",srcPeriod));

}
