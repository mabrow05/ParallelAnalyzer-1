
#include <string>
#include <vector>

// year = 2011 or 2012
void CombinedTwoDim_ErrorPlotter() {

  gStyle->SetTitleSize(0.08,"t");
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadRightMargin(0.18);
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetPadTopMargin(0.05);
  gStyle->SetTitleOffset(1.0,"xy");
  gStyle->SetTitleSize(0.06,"XYZ");
  gStyle->SetTitleOffset(1.1,"Z");
  gStyle->SetLabelSize(0.05,"xyz");
  gStyle->SetOptFit(1111);
  gStyle->SetStatY(0.85);
  gStyle->SetStatX(0.975);
  gStyle->SetStatW(.09);

  TString drawOpt = "colz"; //"lego20"
  //gStyle->SetPalette(kTemperatureMap);

  
  

  gStyle->SetFillStyle(0000); 
  gStyle->SetStatStyle(0);
  gStyle->SetOptStat(0);
  //gStyle->SetTitleStyle(0); 
  //gStyle->SetCanvasBorderSize(0); 
  //gStyle->SetFrameBorderSize(0); 
  gStyle->SetLegendBorderSize(0); 
  //gStyle->SetStatBorderSize(0); 
  //gStyle->SetTitleBorderSize(0);

  gStyle->SetNumberContours(255);

  std::vector < Double_t > asymm;
  std::vector < Double_t > statistics;
  std::vector < Double_t > sim;
  std::vector < Double_t > total;

  std::vector <Int_t> xBinVal;
  std::vector <Int_t> yBinVal;
 
  std::vector <Double_t> en150;
  std::vector <Double_t> en200;
  std::vector <Double_t> en250;
  std::vector <Double_t> en300;
 
  std::vector <Double_t> a150;
  std::vector <Double_t> a200;
  std::vector <Double_t> a250;
  std::vector <Double_t> a300;

  std::vector <Double_t> a150_err;
  std::vector <Double_t> a200_err;
  std::vector <Double_t> a250_err;
  std::vector <Double_t> a300_err;
  
  Double_t as_minval=1000., as_maxval=-1000.;
  Double_t en_minval=1000., en_maxval=-1000.;
  Double_t stat_minval=1000., stat_maxval=-1000.;
  Double_t syst_minval=1000., syst_maxval=-1000.;
  Double_t tot_minval=1000., tot_maxval=-1000.;

  Double_t finalA = 0.;
  Double_t finalA_err=0.;

  std::vector<Double_t> minPlotEn;
  std::vector<Double_t> minPlotStat;
  std::vector<Double_t> minPlotSyst;
  std::vector<Double_t> minPlotTot;
  
  ifstream infile(TString::Format("CombinedUncertainties.txt").Data());

  std::string stars = "";
  std::string windowStr = "";
  std::string txtHold = "";
  Double_t asymm_hold,stat_hold,syst_hold,tot_hold;
  
  while ( infile >> stars 
	  >> windowStr >> tot_hold
	  >> txtHold >> asymm_hold
	  >> txtHold >> asymm_hold
	  >> txtHold >> stat_hold
	  >> txtHold >> syst_hold )  {

    Double_t upperWindow = atof(windowStr.substr(4,3).c_str());
    Double_t lowerWindow = atof(windowStr.substr(0,3).c_str());


    if ( lowerWindow==190. && upperWindow==740. ) finalA = asymm_hold, finalA_err = stat_hold*100.;
    
    if ( lowerWindow==190. ) {
      minPlotEn.push_back(upperWindow);
      minPlotStat.push_back(stat_hold*100.);
      minPlotSyst.push_back(syst_hold*100.);
      minPlotTot.push_back(tot_hold*100.);
    }
    if ( lowerWindow==150. ) {
      std::cout << lowerWindow << endl;
      
      en150.push_back(upperWindow);
      a150.push_back(asymm_hold);
      a150_err.push_back(TMath::Abs(asymm_hold*stat_hold));
    }
    if ( lowerWindow==200. ) {
      en200.push_back(upperWindow);
      a200.push_back(asymm_hold);
      a200_err.push_back(TMath::Abs(asymm_hold*stat_hold));
    }
    if ( lowerWindow==250. ) {
      en250.push_back(upperWindow);
      a250.push_back(asymm_hold);
      a250_err.push_back(TMath::Abs(asymm_hold*stat_hold));
    }
    if ( lowerWindow==300. ) {
      en300.push_back(upperWindow);
      a300.push_back(asymm_hold);
      a300_err.push_back(TMath::Abs(asymm_hold*stat_hold));
    }
    
    tot_hold*=100.;
    syst_hold*=100.;
    stat_hold*=100.;

    if (as_minval>asymm_hold) as_minval=asymm_hold;
    if (as_maxval<asymm_hold) as_maxval=asymm_hold;
    if (stat_minval>stat_hold) stat_minval=stat_hold;
    if (stat_maxval<stat_hold) stat_maxval=stat_hold;
    if (syst_minval>syst_hold) syst_minval=syst_hold;
    if (syst_maxval<syst_hold) syst_maxval=syst_hold;
    if (tot_minval>tot_hold) tot_minval=tot_hold;
    if (tot_maxval<tot_hold) tot_maxval=tot_hold;
    
    xBinVal.push_back( upperWindow );
    yBinVal.push_back( lowerWindow );

    asymm.push_back(asymm_hold);
    statistics.push_back(stat_hold);
    sim.push_back(syst_hold);
    total.push_back(tot_hold);
    
  }
  Double_t ybinLow = 95.;//495.;
  Double_t ybinHigh = 655.;//805.;
  Int_t nbinsY = (ybinHigh-ybinLow)/10.;

  Double_t xbinLow = 195.;//95.;
  Double_t xbinHigh = 755.;//495;
  Int_t nbinsX = (xbinHigh-xbinLow)/10.;
  
  //Two dimensional plots for each error in minimization
  TH2D *as = new TH2D("as","",//"Asymmetry vs. Analysis Window",
		      nbinsX,xbinLow,xbinHigh,nbinsY,ybinLow,ybinHigh);
  for (int i=0; i<=nbinsX; ++i) {
    for (int j=0; j<=nbinsY; ++j) {
      as->SetBinContent(i,j,-3.);
    }
  }
  TH2D *as2 = new TH2D("as2","",//"Asymmetry vs. Analysis Window",
		      nbinsX,xbinLow,xbinHigh,nbinsY,ybinLow,ybinHigh);
  
  TH2D *stat = new TH2D("stat","",//"Statistical Uncertainty vs. Analysis Window",
		      nbinsX,xbinLow,xbinHigh,nbinsY,ybinLow,ybinHigh);

  TH2D *syst = new TH2D("syst","",//"Systematic Uncertainty vs. Analysis Window",
		      nbinsX,xbinLow,xbinHigh,nbinsY,ybinLow,ybinHigh);

  TH2D *tot = new TH2D("tot","",//"Total Uncertainty vs. Analysis Window",
		      nbinsX,xbinLow,xbinHigh,nbinsY,ybinLow,ybinHigh);

  for (UInt_t bin=0; bin<xBinVal.size(); ++bin) {

    Int_t xbin = tot->GetXaxis()->FindBin(xBinVal[bin]);
    Int_t ybin = tot->GetYaxis()->FindBin(yBinVal[bin]);

    // as->SetBinContent(xbin,ybin,(asymm[bin]/finalA>0.99?asymm[bin]/finalA:0.99));
    Double_t x = (asymm[bin]-finalA)/((total[bin]*asymm[bin])/100.); //+finalA_err*finalA
    as->SetBinContent(xbin,ybin,x<-2.?-2.:x);
    as->SetBinError(xbin,ybin,asymm[bin]*statistics[bin]/100./finalA);
    stat->SetBinContent(xbin,ybin,statistics[bin]);
    syst->SetBinContent(xbin,ybin,sim[bin]);
    tot->SetBinContent(xbin,ybin,total[bin]);

    as2->SetBinContent(xbin,ybin,TMath::Abs(x)<1.?1.:0.);

  }


  const Int_t Number=4;
  Double_t Red[Number]    = { 0., 0.00, 1., 1.00};
  Double_t Green[Number]   = { 0., 1.00, 1., 0.};
  Double_t Blue[Number]  = { 1.00, 0.00, 0.00,0.};
  Int_t nb=50;
  
  Double_t Length[Number] = {0.0, 0.33, 0.66, 1.0}; 
  TColor::CreateGradientColorTable(Number,Length,Red,Green,Blue,nb);

    
  TCanvas *cStat = new TCanvas("cStat","cStat");
  //stat->Smooth(5);
  stat->Draw(drawOpt);
  stat->GetXaxis()->SetTitle("Maximum Energy (keV)");
  stat->GetYaxis()->SetTitle("Minimum Energy (keV)");
  stat->GetZaxis()->SetRangeUser(stat_minval*0.9,stat_maxval*0.7);
  //stat->GetZaxis()->SetRangeUser(0.,tot_maxval*0.9);
  stat->GetZaxis()->SetTitle("#deltaA/A (%)");
  stat->GetXaxis()->CenterTitle();
  stat->GetYaxis()->CenterTitle();
  stat->GetZaxis()->CenterTitle();
  
  TCanvas *cSyst = new TCanvas("cSyst","cSyst");
  syst->Draw(drawOpt);
  syst->GetXaxis()->SetTitle("Maximum Energy (keV)");
  syst->GetYaxis()->SetTitle("Minimum Energy (keV)");
  syst->GetZaxis()->SetRangeUser(syst_minval*0.99,syst_maxval*0.9);
  //syst->GetZaxis()->SetRangeUser(0.,tot_maxval*0.9);
  syst->GetZaxis()->SetTitle("#deltaA/A (%)");
  syst->GetXaxis()->CenterTitle();
  syst->GetYaxis()->CenterTitle();
  syst->GetZaxis()->CenterTitle();


  /*Length[1] = tot_minval;
  Length[2] = 0.50*(tot_maxval+tot_minval);
  Length[3] = tot_maxval; 
  TColor::CreateGradientColorTable(Number,Length,Red,Green,Blue,nb);
  */
  
 
  
  TCanvas *cTot = new TCanvas("cTot","cTot");
  tot->Draw(drawOpt);
  tot->GetXaxis()->SetTitle("Maximum Energy (keV)");
  tot->GetYaxis()->SetTitle("Minimum Energy (keV)");
  tot->GetZaxis()->SetRangeUser(tot_minval*0.99,tot_maxval*0.7);
  //tot->GetZaxis()->SetRangeUser(0.,tot_maxval*0.9);
  tot->GetZaxis()->SetTitle("#deltaA/A (%)");
  tot->GetXaxis()->CenterTitle();
  tot->GetYaxis()->CenterTitle();
  tot->GetZaxis()->CenterTitle();

  cStat->Print("TwoDimUncert.pdf(");
  cSyst->Print("TwoDimUncert.pdf");
  cTot->Print("TwoDimUncert.pdf");

  gStyle->SetErrorX(0);

  TEllipse *ell = new TEllipse(740.,190.,10.,10.);
  ell->SetFillStyle(0);
  ell->SetLineColor(1);
  ell->SetLineStyle(2);
  ell->SetLineWidth(2);
  

  TCanvas *cAs = new TCanvas("cAs","cAs");
  gPad->SetRightMargin(0.24);
  as->Draw(drawOpt);
  as->GetXaxis()->SetTitle("Maximum Energy (keV)");
  as->GetYaxis()->SetTitle("Minimum Energy (keV)");
  //as->SetMarkerStyle(25);
  //as->SetMarkerColor(kBlue);
  //as->GetZaxis()->SetRangeUser(0.99,1.01);
  as->GetZaxis()->SetRangeUser(-2.,2.);
  //as->SetMinimum(0.120);
  //as->SetMaximum(0.125);
  as->GetZaxis()->SetTitle("#frac{A#minusA_{min}}{#deltaA}");
  as->GetZaxis()->SetTitleOffset(1.25);
  as->GetXaxis()->CenterTitle();
  as->GetYaxis()->CenterTitle();
  as->GetZaxis()->CenterTitle();

  ell->Draw("same");
  cAs->Update();
  cAs->Print("TwoDimUncert.pdf");


  gStyle->SetPalette(kDeepSea);
  
  TCanvas *cAs2 = new TCanvas("cAs2","cAs2");
  gPad->SetRightMargin(0.05);
  as2->Draw("col");
  as2->GetXaxis()->SetTitle("Maximum Energy (keV)");
  as2->GetYaxis()->SetTitle("Minimum Energy (keV)");
  //as2->SetMarkerStyle(25);
  //as2->SetMarkerColor(kBlue);
  //as2->GetZaxis()->SetRangeUser(0.99,1.01);
  //as2->GetZaxis()->SetRangeUser(-2.,2.);
  //as2->SetMinimum(0.120);
  //as2->SetMaximum(0.125);
  //as2->GetZaxis()->SetTitle("");
  as2->GetXaxis()->CenterTitle();
  as2->GetYaxis()->CenterTitle();
  as2->GetZaxis()->CenterTitle();

  ell->Draw("SAME");
  
  cAs2->Update();
  cAs2->Print("TwoDimUncert.pdf");
  
  


  TCanvas *c10 = new TCanvas("c10","c10");
  TGraph *st = new TGraph(minPlotEn.size(),&minPlotEn[0],&minPlotStat[0]);
  st->SetLineWidth(2);
  st->SetLineStyle(7);
  st->SetLineColor(1);

  TGraph *sys = new TGraph(minPlotEn.size(),&minPlotEn[0],&minPlotSyst[0]);
  sys->SetLineWidth(2);
  sys->SetLineStyle(8);
  sys->SetLineColor(2);
  
  TGraph *t = new TGraph(minPlotEn.size(),&minPlotEn[0],&minPlotTot[0]);
  t->SetLineWidth(3);
  t->SetLineStyle(1);
  t->SetLineColor(4);

  gPad->SetRightMargin(0.05);

  TMultiGraph *minplot = new TMultiGraph();
  minplot->Add(st,"C");
  minplot->Add(sys,"C");
  minplot->Add(t,"C");
  minplot->SetTitle("");//("Systematic and Statistical Uncertainty vs. Upper Analysis Cut");
  minplot->Draw("A");
  minplot->GetXaxis()->SetTitle("Upper Analysis Cut (keV)");
  minplot->GetYaxis()->SetTitle("#deltaA/A (%)");				

  TLegend *l = new TLegend(0.60,0.65,0.88,0.90);
  l->AddEntry(t,"Total","l");
  l->AddEntry(st,"Statistics","l");
  l->AddEntry(sys,"Systematics","l");
  l->Draw("SAME");
  c10->Update();

  
  c10->Print("TwoDimUncert.pdf)");
}
