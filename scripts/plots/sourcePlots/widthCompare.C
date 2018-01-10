

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

  gStyle->SetOptFit(0000);
  gStyle->SetPadLeftMargin(0.17);
  gStyle->SetTitleYOffset(1.2);
  gStyle->SetPadTopMargin(0.05);
  gStyle->SetPadRightMargin(0.05);
  gStyle->SetPadBottomMargin(0.13);
  gStyle->SetTitleXOffset(1.);
  gStyle->SetTitleSize(0.06,"t");
  gStyle->SetLabelSize(0.045,"xy");
  gStyle->SetFillStyle(0000); 


  TCanvas *c1 = new TCanvas("c1","c1",600,800);

  TGraph *gE = new TGraph(dataEastWidth.size(),&dataEastWidth[0],&simEastWidth[0]);
  gE->SetMarkerStyle(kOpenCircle);
  gE->SetMarkerColor(kBlue);
  gE->SetTitle("");//East");
  gE->Draw("AP");
  gE->GetXaxis()->SetLimits(0.,90.);
  gE->SetMaximum(90.);
  gE->SetMinimum(0.);
  gE->GetXaxis()->SetTitle("Data Widths [keV]");
  gE->GetXaxis()->SetTitleSize(0.055);
  gE->GetXaxis()->CenterTitle();
  gE->GetYaxis()->SetTitleSize(0.055);
  gE->GetYaxis()->CenterTitle();
  gE->GetYaxis()->SetTitle("Simulated Widths [keV]");
  gE->Draw("AP");
  c1->Update();


  TF1 *f = new TF1("f","[0]*x",0.,120.);
  f->SetParameter(0,1.);
  gE->Fit(f);
  gPad->Modified();

  TPaveText *pvEast = new TPaveText(0.3,0.8,0.6,0.9,"nbNDC");
  pvEast->SetBorderSize(0);
  pvEast->AddText(TString::Format("slope = %0.3f #pm %0.3f",
				  f->GetParameter(0),f->GetParError(0)));
  pvEast->GetLine(0)->SetTextSize(0.045);
  pvEast->GetLine(0)->SetTextFont(42);
  pvEast->Draw();


  TCanvas *c2 = new TCanvas("c2","c2",600,800);

  TGraph *gW = new TGraph(dataWestWidth.size(),&dataWestWidth[0],&simWestWidth[0]);
  gW->SetMarkerStyle(kOpenCircle);
  gW->SetMarkerColor(kBlue);
  gW->SetTitle("");//West");
  gW->Draw("AP");
  gW->GetXaxis()->SetLimits(0.,90.);
  gW->SetMaximum(90.);
  gW->SetMinimum(0.);
  gW->GetXaxis()->SetTitle("Data Widths [keV]");
  gW->GetYaxis()->SetTitle("Simulated Widths [keV]");
  gW->GetXaxis()->SetTitleSize(0.055);
  gW->GetYaxis()->SetTitleSize(0.055);
  gW->GetXaxis()->CenterTitle();
  gW->GetYaxis()->CenterTitle();
  gW->Draw("AP");
  c2->Update();



  gW->Fit(f);
  gPad->Modified();

  TPaveText *pvWest = new TPaveText(0.3,0.8,0.6,0.9,"nbNDC");
  pvWest->SetBorderSize(0);
  pvWest->AddText(TString::Format("slope = %0.3f #pm %0.3f",
				  f->GetParameter(0),f->GetParError(0)));
  pvWest->GetLine(0)->SetTextSize(0.045);
  pvWest->GetLine(0)->SetTextFont(42);
  pvWest->Draw();

  c1->Print(TString::Format("widthsErecon_runPeriod_%i.pdf(",srcPeriod));
  c2->Print(TString::Format("widthsErecon_runPeriod_%i.pdf)",srcPeriod));

}
