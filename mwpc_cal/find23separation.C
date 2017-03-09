

void plot23separation(TString geom) {

  //gStyle->SetOptStat(0);
  //gStyle->SetStatFontSize(0.030);
  //gStyle->SetOptFit(1111);
  //gStyle->SetOptTitle(0);
  //gStyle->SetTitleSize(0.08,"t");
  //gStyle->SetTitleX(0.17);
  //gStyle->SetTitleAlign(13);
  gStyle->SetTitleOffset(0.90, "x");
  //gStyle->SetTitleOffset(1.30, "y");
  //gStyle->SetPadTickX(1);
  //gStyle->SetPadTickY(1);
  //gStyle->SetNdivisions(510,"X");
  //gStyle->SetNdivisions(510,"Y");
  //gStyle->SetNdivisions(9,"Z");
  //gStyle->SetPadLeftMargin(0.13); // 0.13
  //gStyle->SetPadRightMargin(0.15); // 0.04
  //gStyle->SetPadBottomMargin(0.37); // 0.30
  //gStyle->SetLabelSize(0.045, "X");
  //gStyle->SetLabelSize(0.045, "Y");
  //gStyle->SetLabelSize(0.045, "Z");
  //gStyle->SetLabelOffset(0.00, "X");
  //gStyle->SetLabelOffset(0.01, "Y");
  gStyle->SetTitleSize(0.050, "X");
  gStyle->SetTitleSize(0.050, "Y");

  //Borders on Legends and titles and stats
  gStyle->SetFillStyle(0000); 
  //gStyle->SetStatStyle(0); 
  //gStyle->SetTitleStyle(0); 
  //gStyle->SetCanvasBorderSize(0); 
  //gStyle->SetFrameBorderSize(0); 
  gStyle->SetLegendBorderSize(0); 
  //gStyle->SetStatBorderSize(0); 
  //gStyle->SetTitleBorderSize(0);
  
  int octMin = 0;
  int octMax = 0;
  double sepVal = 6.6;
  
  if ( geom==TString("2011-2012") ) octMin = 0, octMax = 59;
  else if ( geom==TString("2012-2013") ) octMin = 80, octMax = 121;
  else if ( geom==TString("2012-2013_isobutane") ) octMin = 67, octMax = 79, sepVal = 5.4 ; 
  else { std::cout << "Bad geometry\n"; exit(0); }

  TH1D *hType2E = new TH1D();
  TH1D *hType3E = new TH1D();
  TH1D *hType2W = new TH1D();
  TH1D *hType3W = new TH1D();
  
  TH1D *hScint2E = new TH1D();
  TH1D *hScint3E = new TH1D();
  TH1D *hScint2W = new TH1D();
  TH1D *hScint3W = new TH1D();

  TH1D::AddDirectory(false);
  
  TFile *f = new TFile(TString::Format("Backscatters_%i-%i.root",octMin,octMax),
		       "READ");

  hType2E = (TH1D*)f->Get("hType2E")->Clone();
  hType3E = (TH1D*)f->Get("hType3E")->Clone();
  hType2W = (TH1D*)f->Get("hType2W")->Clone();
  hType3W = (TH1D*)f->Get("hType3W")->Clone();

  hScint2E = (TH1D*)f->Get("hScint2E")->Clone();
  hScint3E = (TH1D*)f->Get("hScint3E")->Clone();
  hScint2W = (TH1D*)f->Get("hScint2W")->Clone();
  hScint3W = (TH1D*)f->Get("hScint3W")->Clone();

  delete f;

  //
  
  TCanvas *c1 = new TCanvas("c1","c1",1000,600);
  c1->Divide(2,1);

  c1->cd(1);

  hType2E->SetTitle("East Type II/III Separation E_{MWPC}");
  hType2E->GetXaxis()->SetTitle("E_{MWPC} (keV)");
  hType2E->GetXaxis()->CenterTitle();
  hType2E->SetLineColor(kBlue);
  hType2E->SetLineWidth(2);
  hType2E->Draw();

  hType3E->SetLineColor(kBlue);
  hType3E->SetLineWidth(2);
  hType3E->SetLineStyle(7);
  hType3E->Draw("SAME");
  c1->Update();

  TLegend *leg1 = new TLegend(0.55,0.75,0.85,0.85);
  leg1->AddEntry(hType2E,"Type II");
  leg1->AddEntry(hType3E,"Type III");
  leg1->Draw("SAME");
  
  TLine *line = new TLine(sepVal,0.,sepVal,gPad->GetUymax());
  line->SetLineWidth(2);
  line->SetLineStyle(7);
  line->Draw("SAME");
  
  c1->cd(2);

  
  hType2W->SetTitle("West Type II/III Separation E_{MWPC}");
  hType2W->GetXaxis()->SetTitle("E_{MWPC} (keV)");
  hType2W->GetXaxis()->CenterTitle();
  hType2W->SetLineColor(kBlue);
  hType2W->SetLineWidth(2);
  hType2W->Draw();

  hType3W->SetLineColor(kBlue);
  hType3W->SetLineWidth(2);
  hType3W->SetLineStyle(7);
  hType3W->Draw("SAME");
  c1->Update();

  TLegend *leg2 = new TLegend(0.55,0.75,0.85,0.85);
  leg2->AddEntry(hType2W,"Type II");
  leg2->AddEntry(hType3W,"Type III");
  leg2->Draw("SAME");

  line = new TLine(sepVal,0.,sepVal,gPad->GetUymax());
  line->SetLineWidth(2);
  line->SetLineStyle(7);
  line->Draw("SAME");

  
  

}
