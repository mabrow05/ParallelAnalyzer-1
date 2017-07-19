
#include <fstream>
#include <iostream>
#include <vector>
#include <TString.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TAxis.h>
#include <TLegend.h>
#include <TMultiGraph.h>
#include <TMath.h>

std::vector <std::vector<Double_t> > readBSCorr(TString corr,TString year) {

  std::vector <std::vector<Double_t> > corrs;
  std::ifstream infile(TString::Format("../OldMCCorrection/deltaBS%s_anaChC_%s.txt",corr.Data(),year.Data()));
  std::vector<Double_t> enbin;
  std::vector<Double_t> c;
  std::vector<Double_t> cerr;
  Double_t en, corrhold,errhold;
  Int_t i=0;
  Double_t aveEn=0., aveCorr=0.;
  while (infile >> en >> corrhold >> errhold && en<800.) {
    //std::cout << en << " " << corrhold << " " << errhold << "\n";
    if (i<4) {
      aveEn+=en;
      aveCorr+=corrhold;
      i++;
    } else {
      enbin.push_back(aveEn/4.);
      c.push_back(aveCorr/4.*100.);
      cerr.push_back(0.25*aveCorr/4.*100.);
      i=0,aveEn=0.,aveCorr=0.;
    }
  }
  corrs.push_back(enbin);
  corrs.push_back(c);
  corrs.push_back(cerr);
  return corrs;
};

std::vector <std::vector<Double_t> > readAngleCorr(TString year) {

  std::vector <std::vector<Double_t> > corrs;
  std::ifstream infile(TString::Format("../AngleCorrections/%s_DeltaAngle_anaChC.txt",year.Data()));
  std::vector<Double_t> enbin;
  std::vector<Double_t> c;
  std::vector<Double_t> cerr;
  Double_t en, frachold,corrhold;
  Int_t i=0;
  Double_t aveEn=0., aveCorr=0.;
  while (infile >> en >> frachold >> corrhold && en<800.) {
    //std::cout << en << " " << corrhold << " " << errhold << "\n";
    if (i<4) {
      aveEn+=en;
      aveCorr+=corrhold;
      i++;
    } else {
      enbin.push_back(aveEn/4.);
      c.push_back(aveCorr/4.*100.);
      cerr.push_back(0.25*aveCorr/4.*100.);
      i=0,aveEn=0.,aveCorr=0.;
    }
  }
  corrs.push_back(enbin);
  corrs.push_back(c);
  corrs.push_back(cerr);
  return corrs;
};

void EnergyDependentCorrections() {

  gStyle->SetOptStat(0);
  gStyle->SetTitleSize(0.07,"t");
  gStyle->SetPadLeftMargin(0.12);
  //gStyle->SetPadRightMargin(0.15);
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetTitleYSize(0.05);
  gStyle->SetTitleYOffset(1.);
  gStyle->SetTitleXSize(0.05);
  gStyle->SetTitleXOffset(1.2);
  gStyle->SetTitleSize(0.04,"Z");
  gStyle->SetTitleOffset(1.6,"Z");
  gStyle->SetLabelSize(0.04,"xyz");
  gStyle->SetLegendBorderSize(0);
  gStyle->SetFillStyle(0);
  gStyle->SetGridStyle(2);
  //  gStyle->SetHatchesSpacing(0.6);
  //gStyle->SetHatchesLineWidth(2);
  // read in the corrections and plot their deltaA/A for each bin
  // This can be done for both theory and MC

  //Int_t octmin = year==2011?0:60;
  //Int_t octmax = year==2011?59:121;

  int fill2011 = 3001;
  int fill2012 = 3001;

  int col2011 = 4;
  int col2012 = 3;
  
  TString year = "2011-2012";
  std::vector<std::vector<Double_t> > bs0_2011 = readBSCorr("0",year);
  std::vector<std::vector<Double_t> > bs1_2011 = readBSCorr("1",year);
  std::vector<std::vector<Double_t> > bs2_2011 = readBSCorr("2",year);
  std::vector<std::vector<Double_t> > bs3_2011 = readBSCorr("3",year);
  std::vector<std::vector<Double_t> > bsALL_2011 = readBSCorr("ALL",year);
  std::vector<std::vector<Double_t> > cosTheta_2011 = readAngleCorr(year);

  year = "2012-2013";
  std::vector<std::vector<Double_t> > bs0_2012 = readBSCorr("0",year);
  std::vector<std::vector<Double_t> > bs1_2012 = readBSCorr("1",year);
  std::vector<std::vector<Double_t> > bs2_2012 = readBSCorr("2",year);
  std::vector<std::vector<Double_t> > bs3_2012 = readBSCorr("3",year);
  std::vector<std::vector<Double_t> > bsALL_2012 = readBSCorr("ALL",year);
  std::vector<std::vector<Double_t> > cosTheta_2012 = readAngleCorr(year);

  std::vector<std::vector<Double_t> > totalCorr_2011(3,std::vector<Double_t>(0));
  std::vector<std::vector<Double_t> > totalCorr_2012(3,std::vector<Double_t>(0));

  for (UInt_t i=0;i<bsALL_2011[0].size();++i) {
    Double_t c2011 = ((1.+bsALL_2011[1][i]/100.)*(1.+cosTheta_2011[1][i]/100.)-1.)*100.;
    Double_t c2011err = TMath::Sqrt(TMath::Power((1.+cosTheta_2011[1][i]/100.)*bsALL_2011[2][i]/100.,2.)
				    +TMath::Power((1.+bsALL_2011[1][i]/100.)*cosTheta_2011[2][i]/100.,2.))*100.;
    std::cout << bsALL_2011[0][i] << " " << c2011 << " " << c2011err << std::endl;
    totalCorr_2011[0].push_back(bsALL_2011[0][i]);
    totalCorr_2011[1].push_back(c2011);
    totalCorr_2011[2].push_back(c2011err);

    Double_t c2012 = ((1.+bsALL_2012[1][i]/100.)*(1.+cosTheta_2012[1][i]/100.)-1.)*100.;
    Double_t c2012err = TMath::Sqrt(TMath::Power((1.+cosTheta_2012[1][i]/100.)*bsALL_2012[2][i]/100.,2.)
				    +TMath::Power((1.+bsALL_2012[1][i]/100.)*cosTheta_2012[2][i]/100.,2.))*100.;
    totalCorr_2012[0].push_back(bsALL_2012[0][i]);
    totalCorr_2012[1].push_back(c2012);
    totalCorr_2012[2].push_back(c2012err);    
  }
  
  TCanvas *c0 = new TCanvas("c0","c0");
  gPad->SetGrid(0,1);

  TMultiGraph *mg0 = new TMultiGraph();
  mg0->SetTitle("#Delta_{BS,0} vs. Energy");  
  
  TGraphErrors *g_bs0_2011 = new TGraphErrors(bs0_2011[0].size(),&bs0_2011[0][0],&bs0_2011[1][0],0,&bs0_2011[2][0]);
  g_bs0_2011->SetMarkerStyle(0);
  g_bs0_2011->SetLineWidth(3);
  g_bs0_2011->SetLineColor(col2011);
  g_bs0_2011->SetFillColor(col2011);
  g_bs0_2011->SetFillStyle(fill2011);
  
  TGraphErrors *g_bs0_2012 = new TGraphErrors(bs0_2012[0].size(),&bs0_2012[0][0],&bs0_2012[1][0],0,&bs0_2012[2][0]);
  g_bs0_2012->SetMarkerStyle(0);
  g_bs0_2012->SetLineWidth(3);
  g_bs0_2012->SetLineColor(col2012);
  g_bs0_2012->SetFillColor(col2012);
  g_bs0_2012->SetFillStyle(fill2012);
  
  mg0->Add(g_bs0_2011);
  mg0->Add(g_bs0_2012);
  mg0->SetMinimum(-5.);
  mg0->SetMaximum(5.);
  
  mg0->Draw("ALP3");
  mg0->GetYaxis()->SetTitle("#DeltaA/A (%)");
  mg0->GetYaxis()->CenterTitle();
  mg0->GetXaxis()->SetTitle("Energy (keV)");
  mg0->GetXaxis()->CenterTitle();
  gPad->Modified();

  TLegend *leg0 = new TLegend(0.57,0.7,0.87,0.8);
  leg0->AddEntry(g_bs0_2011,"2011-2012","lf");
  leg0->AddEntry(g_bs0_2012,"2012-2013","lf");
  leg0->SetTextSize(0.05);
  leg0->Draw("SAME");


  //////////////////// BS1///////////////////
  TCanvas *c1 = new TCanvas("c1","c1");
  gPad->SetGrid(0,1);

  TMultiGraph *mg1 = new TMultiGraph();
  mg1->SetTitle("#Delta_{BS,1} vs. Energy");  
  
  TGraphErrors *g_bs1_2011 = new TGraphErrors(bs1_2011[0].size(),&bs1_2011[0][0],&bs1_2011[1][0],0,&bs1_2011[2][0]);
  g_bs1_2011->SetMarkerStyle(0);
  g_bs1_2011->SetLineWidth(3);
  g_bs1_2011->SetLineColor(col2011);
  g_bs1_2011->SetFillColor(col2011);
  g_bs1_2011->SetFillStyle(fill2011);
  
  TGraphErrors *g_bs1_2012 = new TGraphErrors(bs1_2012[0].size(),&bs1_2012[0][0],&bs1_2012[1][0],0,&bs1_2012[2][0]);
  g_bs1_2012->SetMarkerStyle(0);
  g_bs1_2012->SetLineWidth(3);
  g_bs1_2012->SetLineColor(col2012);
  g_bs1_2012->SetFillColor(col2012);
  g_bs1_2012->SetFillStyle(fill2012);
  
  mg1->Add(g_bs1_2011);
  mg1->Add(g_bs1_2012);
  mg1->SetMinimum(-5.);
  mg1->SetMaximum(5.);
  
  mg1->Draw("ALP3");
  mg1->GetYaxis()->SetTitle("#DeltaA/A (%)");
  mg1->GetYaxis()->CenterTitle();
  mg1->GetXaxis()->SetTitle("Energy (keV)");
  mg1->GetXaxis()->CenterTitle();
  gPad->Modified();

  TLegend *leg1 = new TLegend(0.57,0.7,0.87,0.8);
  leg1->AddEntry(g_bs1_2011,"2011-2012","lf");
  leg1->AddEntry(g_bs1_2012,"2012-2013","lf");
  leg1->SetTextSize(0.05);
  leg1->Draw("SAME");

  //////////////////// BS2///////////////////
  TCanvas *c2 = new TCanvas("c2","c2");
  gPad->SetGrid(0,1);

  TMultiGraph *mg2 = new TMultiGraph();
  mg2->SetTitle("#Delta_{BS,2} vs. Energy");  
  
  TGraphErrors *g_bs2_2011 = new TGraphErrors(bs2_2011[0].size(),&bs2_2011[0][0],&bs2_2011[1][0],0,&bs2_2011[2][0]);
  g_bs2_2011->SetMarkerStyle(0);
  g_bs2_2011->SetLineWidth(3);
  g_bs2_2011->SetLineColor(col2011);
  g_bs2_2011->SetFillColor(col2011);
  g_bs2_2011->SetFillStyle(fill2011);
  
  TGraphErrors *g_bs2_2012 = new TGraphErrors(bs2_2012[0].size(),&bs2_2012[0][0],&bs2_2012[1][0],0,&bs2_2012[2][0]);
  g_bs2_2012->SetMarkerStyle(0);
  g_bs2_2012->SetLineWidth(3);
  g_bs2_2012->SetLineColor(col2012);
  g_bs2_2012->SetFillColor(col2012);
  g_bs2_2012->SetFillStyle(fill2012);
  
  mg2->Add(g_bs2_2011);
  mg2->Add(g_bs2_2012);
  mg2->SetMinimum(-5.);
  mg2->SetMaximum(5.);
  
  mg2->Draw("ALP3");
  mg2->GetYaxis()->SetTitle("#DeltaA/A (%)");
  mg2->GetYaxis()->CenterTitle();
  mg2->GetXaxis()->SetTitle("Energy (keV)");
  mg2->GetXaxis()->CenterTitle();
  gPad->Modified();

  TLegend *leg2 = new TLegend(0.57,0.7,0.87,0.8);
  leg2->AddEntry(g_bs2_2011,"2011-2012","lf");
  leg2->AddEntry(g_bs2_2012,"2012-2013","lf");
  leg2->SetTextSize(0.05);
  leg2->Draw("SAME");

  //////////////////// BS3///////////////////
  TCanvas *c3 = new TCanvas("c3","c3");
  gPad->SetGrid(0,1);

  TMultiGraph *mg3 = new TMultiGraph();
  mg3->SetTitle("#Delta_{BS,3} vs. Energy");  
  
  TGraphErrors *g_bs3_2011 = new TGraphErrors(bs3_2011[0].size(),&bs3_2011[0][0],&bs3_2011[1][0],0,&bs3_2011[2][0]);
  g_bs3_2011->SetMarkerStyle(0);
  g_bs3_2011->SetLineWidth(3);
  g_bs3_2011->SetLineColor(col2011);
  g_bs3_2011->SetFillColor(col2011);
  g_bs3_2011->SetFillStyle(fill2011);
  
  TGraphErrors *g_bs3_2012 = new TGraphErrors(bs3_2012[0].size(),&bs3_2012[0][0],&bs3_2012[1][0],0,&bs3_2012[2][0]);
  g_bs3_2012->SetMarkerStyle(0);
  g_bs3_2012->SetLineWidth(3);
  g_bs3_2012->SetLineColor(col2012);
  g_bs3_2012->SetFillColor(col2012);
  g_bs3_2012->SetFillStyle(fill2012);
  
  mg3->Add(g_bs3_2011);
  mg3->Add(g_bs3_2012);
  mg3->SetMinimum(-5.);
  mg3->SetMaximum(5.);
  
  mg3->Draw("ALP3");
  mg3->GetYaxis()->SetTitle("#DeltaA/A (%)");
  mg3->GetYaxis()->CenterTitle();
  mg3->GetXaxis()->SetTitle("Energy (keV)");
  mg3->GetXaxis()->CenterTitle();
  gPad->Modified();

  TLegend *leg3 = new TLegend(0.57,0.7,0.87,0.8);
  leg3->AddEntry(g_bs3_2011,"2011-2012","lf");
  leg3->AddEntry(g_bs3_2012,"2012-2013","lf");
  leg3->SetTextSize(0.05);
  leg3->Draw("SAME");
  

  //////////////////// BSALL///////////////////
  TCanvas *cALLBS = new TCanvas("cALLBS","cALLBS");
  gPad->SetGrid(0,1);

  TMultiGraph *mgALLBS = new TMultiGraph();
  mgALLBS->SetTitle("#Delta_{BS} vs. Energy");  
  
  TGraphErrors *g_bsALL_2011 = new TGraphErrors(bsALL_2011[0].size(),&bsALL_2011[0][0],&bsALL_2011[1][0],0,&bsALL_2011[2][0]);
  g_bsALL_2011->SetMarkerStyle(0);
  g_bsALL_2011->SetLineWidth(3);
  g_bsALL_2011->SetLineColor(col2011);
  g_bsALL_2011->SetFillColor(col2011);
  g_bsALL_2011->SetFillStyle(fill2011);
  
  TGraphErrors *g_bsALL_2012 = new TGraphErrors(bsALL_2012[0].size(),&bsALL_2012[0][0],&bsALL_2012[1][0],0,&bsALL_2012[2][0]);
  g_bsALL_2012->SetMarkerStyle(0);
  g_bsALL_2012->SetLineWidth(3);
  g_bsALL_2012->SetLineColor(col2012);
  g_bsALL_2012->SetFillColor(col2012);
  g_bsALL_2012->SetFillStyle(fill2012);
  
  mgALLBS->Add(g_bsALL_2011);
  mgALLBS->Add(g_bsALL_2012);
  mgALLBS->SetMinimum(-5.);
  mgALLBS->SetMaximum(5.);
  
  mgALLBS->Draw("ALP3");
  mgALLBS->GetYaxis()->SetTitle("#DeltaA/A (%)");
  mgALLBS->GetYaxis()->CenterTitle();
  mgALLBS->GetXaxis()->SetTitle("Energy (keV)");
  mgALLBS->GetXaxis()->CenterTitle();
  gPad->Modified();

  TLegend *legALLBS = new TLegend(0.57,0.7,0.87,0.8);
  legALLBS->AddEntry(g_bsALL_2011,"2011-2012","lf");
  legALLBS->AddEntry(g_bsALL_2012,"2012-2013","lf");
  legALLBS->SetTextSize(0.05);
  legALLBS->Draw("SAME");

  //////////////////// Angle ///////////////////
  TCanvas *cAngle = new TCanvas("cAngle","cAngle");
  gPad->SetGrid(0,1);

  TMultiGraph *mgAngle = new TMultiGraph();
  mgAngle->SetTitle("#Delta_{cos#theta} vs. Energy");  
  
  TGraphErrors *g_cosTheta_2011 = new TGraphErrors(cosTheta_2011[0].size(),&cosTheta_2011[0][0],&cosTheta_2011[1][0],0,&cosTheta_2011[2][0]);
  g_cosTheta_2011->SetMarkerStyle(0);
  g_cosTheta_2011->SetLineWidth(3);
  g_cosTheta_2011->SetLineColor(col2011);
  g_cosTheta_2011->SetFillColor(col2011);
  g_cosTheta_2011->SetFillStyle(fill2011);
  
  TGraphErrors *g_cosTheta_2012 = new TGraphErrors(cosTheta_2012[0].size(),&cosTheta_2012[0][0],&cosTheta_2012[1][0],0,&cosTheta_2012[2][0]);
  g_cosTheta_2012->SetMarkerStyle(0);
  g_cosTheta_2012->SetLineWidth(3);
  g_cosTheta_2012->SetLineColor(col2012);
  g_cosTheta_2012->SetFillColor(col2012);
  g_cosTheta_2012->SetFillStyle(fill2012);
  
  mgAngle->Add(g_cosTheta_2011);
  mgAngle->Add(g_cosTheta_2012);
  mgAngle->SetMinimum(-5.);
  mgAngle->SetMaximum(5.);
  
  mgAngle->Draw("ALP3");
  mgAngle->GetYaxis()->SetTitle("#DeltaA/A (%)");
  mgAngle->GetYaxis()->CenterTitle();
  mgAngle->GetXaxis()->SetTitle("Energy (keV)");
  mgAngle->GetXaxis()->CenterTitle();
  gPad->Modified();

  TLegend *legAngle = new TLegend(0.57,0.7,0.87,0.8);
  legAngle->AddEntry(g_cosTheta_2011,"2011-2012","lf");
  legAngle->AddEntry(g_cosTheta_2012,"2012-2013","lf");
  legAngle->SetTextSize(0.05);
  legAngle->Draw("SAME");

  //////////////////// TotalCorr ///////////////////
  TCanvas *cTotalCorr = new TCanvas("cTotalCorr","cTotalCorr");
  gPad->SetGrid(0,1);

  TMultiGraph *mgTotalCorr = new TMultiGraph();
  mgTotalCorr->SetTitle("#Delta_{MC} vs. Energy");  
  
  TGraphErrors *g_totalCorr_2011 = new TGraphErrors(totalCorr_2011[0].size(),&totalCorr_2011[0][0],&totalCorr_2011[1][0],0,&totalCorr_2011[2][0]);
  g_totalCorr_2011->SetMarkerStyle(0);
  g_totalCorr_2011->SetLineWidth(3);
  g_totalCorr_2011->SetLineColor(col2011);
  g_totalCorr_2011->SetFillColor(col2011);
  g_totalCorr_2011->SetFillStyle(fill2011);
  
  TGraphErrors *g_totalCorr_2012 = new TGraphErrors(totalCorr_2012[0].size(),&totalCorr_2012[0][0],&totalCorr_2012[1][0],0,&totalCorr_2012[2][0]);
  g_totalCorr_2012->SetMarkerStyle(0);
  g_totalCorr_2012->SetLineWidth(3);
  g_totalCorr_2012->SetLineColor(col2012);
  g_totalCorr_2012->SetFillColor(col2012);
  g_totalCorr_2012->SetFillStyle(fill2012);
  
  mgTotalCorr->Add(g_totalCorr_2011);
  mgTotalCorr->Add(g_totalCorr_2012);
  mgTotalCorr->SetMinimum(-5.);
  mgTotalCorr->SetMaximum(5.);
  
  mgTotalCorr->Draw("ALP3");
  mgTotalCorr->GetYaxis()->SetTitle("#DeltaA/A (%)");
  mgTotalCorr->GetYaxis()->CenterTitle();
  mgTotalCorr->GetXaxis()->SetTitle("Energy (keV)");
  mgTotalCorr->GetXaxis()->CenterTitle();
  gPad->Modified();

  TLegend *legTotalCorr = new TLegend(0.57,0.7,0.87,0.8);
  legTotalCorr->AddEntry(g_totalCorr_2011,"2011-2012","lf");
  legTotalCorr->AddEntry(g_totalCorr_2012,"2012-2013","lf");
  legTotalCorr->SetTextSize(0.05);
  legTotalCorr->Draw("SAME");
  

  TString pdffile = "MC_Corrections.pdf";

  c0->Print(TString::Format("%s(",pdffile.Data()));
  c1->Print(pdffile);
  c2->Print(pdffile);
  c3->Print(pdffile);
  cALLBS->Print(pdffile);
  cAngle->Print(pdffile);
  cTotalCorr->Print(TString::Format("%s)",pdffile.Data()));
  
}
