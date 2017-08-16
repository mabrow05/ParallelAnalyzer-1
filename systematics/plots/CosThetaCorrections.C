
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

Int_t groupBin=1;
Double_t enStart=40.;//10.*groupBin;

Double_t BSlimitLow = -1., BSlimitHigh = 3.;
Double_t AnglelimitLow = -5., AnglelimitHigh = 3.;

TString drawOpt = "0AL3";//"0AC4" "0AL3"

std::vector <std::vector<Double_t> > readAngleCorr(TString year, TString type) {

  std::vector <std::vector<Double_t> > corrs;
  std::ifstream infile(TString::Format("../AngleCorrections/%s_delta_3%s.txt",year.Data(),(type==TString("ALL")?"":type.Data())));
  std::vector<Double_t> enbin;
  std::vector<Double_t> c;
  std::vector<Double_t> cerr;
  TString hold1,hold2;
  Double_t en;
  Int_t i=0;
  Double_t aveEn=0., aveCorr=0., aveErr=0.;
  while (infile >> en >> hold1 >> hold2 && en<800.) {
    if (en<enStart) continue;
    //std::cout << en << " " << hold2 << " " << errhold << "\n";
    if (i<groupBin) {
      aveEn+=en;
      aveCorr+=( hold1==TString("nan")||hold1==TString("-nan"))?0.:atof(hold1.Data());
      aveErr+=(hold1==TString("nan")||hold1==TString("-nan"))?0.:atof(hold2.Data());
      i++;
    } else {
      enbin.push_back(aveEn/groupBin);
      c.push_back(aveCorr/groupBin*100.);
      cerr.push_back(aveErr/groupBin*100.);
      i=0,aveEn=0.,aveCorr=0., aveErr=0.;
      aveEn+=en;
      aveCorr+=( hold1==TString("nan")||hold1==TString("-nan"))?0.:atof(hold1.Data());
      aveErr+=(hold1==TString("nan")||hold1==TString("-nan"))?0.:atof(hold2.Data());
      i++;
    }
  }
  corrs.push_back(enbin);
  corrs.push_back(c);
  corrs.push_back(cerr);
  return corrs;
};

void CosThetaCorrections() {

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

  int startPoint=0;
  
  TString year = "2011-2012";
  std::vector<std::vector<Double_t> > delta30_2011 = readAngleCorr(year,"0");
  std::vector<std::vector<Double_t> > delta31_2011 = readAngleCorr(year,"1");
  std::vector<std::vector<Double_t> > delta32_2011 = readAngleCorr(year,"2");
  std::vector<std::vector<Double_t> > delta33_2011 = readAngleCorr(year,"3");
  std::vector<std::vector<Double_t> > delta3_2011 = readAngleCorr(year,"ALL");
    
  year = "2012-2013";
  std::vector<std::vector<Double_t> > delta30_2012 = readAngleCorr(year,"0");
  std::vector<std::vector<Double_t> > delta31_2012 = readAngleCorr(year,"1");
  std::vector<std::vector<Double_t> > delta32_2012 = readAngleCorr(year,"2");
  std::vector<std::vector<Double_t> > delta33_2012 = readAngleCorr(year,"3");
  std::vector<std::vector<Double_t> > delta3_2012 = readAngleCorr(year,"ALL");
  
  
  TCanvas *c0 = new TCanvas("c0","c0");
  gPad->SetGrid(0,1);

  TMultiGraph *mg0 = new TMultiGraph();
  mg0->SetTitle("#Delta_{3,0} vs. Energy");  
  
  TGraphErrors *g_delta30_2011 = new TGraphErrors(delta30_2011[0].size()-startPoint,&delta30_2011[0][startPoint],&delta30_2011[1][startPoint],0,&delta30_2011[2][startPoint]);
  g_delta30_2011->SetMarkerStyle(0);
  g_delta30_2011->SetLineWidth(3);
  g_delta30_2011->SetLineColor(col2011);
  g_delta30_2011->SetFillColor(col2011);
  g_delta30_2011->SetFillStyle(fill2011);
  
  TGraphErrors *g_delta30_2012 = new TGraphErrors(delta30_2012[0].size()-startPoint,&delta30_2012[0][startPoint],&delta30_2012[1][startPoint],0,&delta30_2012[2][startPoint]);
  g_delta30_2012->SetMarkerStyle(0);
  g_delta30_2012->SetLineWidth(3);
  g_delta30_2012->SetLineColor(col2012);
  g_delta30_2012->SetFillColor(col2012);
  g_delta30_2012->SetFillStyle(fill2012);
  
  mg0->Add(g_delta30_2011);
  mg0->Add(g_delta30_2012);
  mg0->SetMinimum(AnglelimitLow);
  mg0->SetMaximum(AnglelimitHigh);
  
  mg0->Draw("ALP3");
  mg0->GetYaxis()->SetTitle("#DeltaA/A (%)");
  mg0->GetYaxis()->CenterTitle();
  mg0->GetXaxis()->SetTitle("Energy (keV)");
  mg0->GetXaxis()->CenterTitle();
  gPad->Modified();

  TLegend *leg0 = new TLegend(0.57,0.7,0.87,0.8);
  leg0->AddEntry(g_delta30_2011,"2011-2012","lf");
  leg0->AddEntry(g_delta30_2012,"2012-2013","lf");
  leg0->SetTextSize(0.05);
  leg0->Draw("SAME");


  //////////////////// DELTA31///////////////////
  TCanvas *c1 = new TCanvas("c1","c1");
  gPad->SetGrid(0,1);

  TMultiGraph *mg1 = new TMultiGraph();
  mg1->SetTitle("#Delta_{3,1} vs. Energy");  
  
  TGraphErrors *g_delta31_2011 = new TGraphErrors(delta31_2011[0].size()-startPoint,&delta31_2011[0][startPoint],&delta31_2011[1][startPoint],0,&delta31_2011[2][startPoint]);
  g_delta31_2011->SetMarkerStyle(0);
  g_delta31_2011->SetLineWidth(3);
  g_delta31_2011->SetLineColor(col2011);
  g_delta31_2011->SetFillColor(col2011);
  g_delta31_2011->SetFillStyle(fill2011);
  
  TGraphErrors *g_delta31_2012 = new TGraphErrors(delta31_2012[0].size()-startPoint,&delta31_2012[0][startPoint],&delta31_2012[1][startPoint],0,&delta31_2012[2][startPoint]);
  g_delta31_2012->SetMarkerStyle(0);
  g_delta31_2012->SetLineWidth(3);
  g_delta31_2012->SetLineColor(col2012);
  g_delta31_2012->SetFillColor(col2012);
  g_delta31_2012->SetFillStyle(fill2012);
  
  mg1->Add(g_delta31_2011);
  mg1->Add(g_delta31_2012);
  mg1->SetMinimum(AnglelimitLow);
  mg1->SetMaximum(AnglelimitHigh);
  
  mg1->Draw("ALP3");
  mg1->GetYaxis()->SetTitle("#DeltaA/A (%)");
  mg1->GetYaxis()->CenterTitle();
  mg1->GetXaxis()->SetTitle("Energy (keV)");
  mg1->GetXaxis()->CenterTitle();
  gPad->Modified();

  TLegend *leg1 = new TLegend(0.57,0.7,0.87,0.8);
  leg1->AddEntry(g_delta31_2011,"2011-2012","lf");
  leg1->AddEntry(g_delta31_2012,"2012-2013","lf");
  leg1->SetTextSize(0.05);
  leg1->Draw("SAME");

  //////////////////// DELTA32///////////////////
  TCanvas *c2 = new TCanvas("c2","c2");
  gPad->SetGrid(0,1);

  TMultiGraph *mg2 = new TMultiGraph();
  mg2->SetTitle("#Delta_{3,2} vs. Energy");  
  
  TGraphErrors *g_delta32_2011 = new TGraphErrors(delta32_2011[0].size()-startPoint,&delta32_2011[0][startPoint],&delta32_2011[1][startPoint],0,&delta32_2011[2][startPoint]);
  g_delta32_2011->SetMarkerStyle(0);
  g_delta32_2011->SetLineWidth(3);
  g_delta32_2011->SetLineColor(col2011);
  g_delta32_2011->SetFillColor(col2011);
  g_delta32_2011->SetFillStyle(fill2011);
  
  TGraphErrors *g_delta32_2012 = new TGraphErrors(delta32_2012[0].size()-startPoint,&delta32_2012[0][startPoint],&delta32_2012[1][startPoint],0,&delta32_2012[2][startPoint]);
  g_delta32_2012->SetMarkerStyle(0);
  g_delta32_2012->SetLineWidth(3);
  g_delta32_2012->SetLineColor(col2012);
  g_delta32_2012->SetFillColor(col2012);
  g_delta32_2012->SetFillStyle(fill2012);
  
  mg2->Add(g_delta32_2011);
  mg2->Add(g_delta32_2012);
  mg2->SetMinimum(AnglelimitLow);
  mg2->SetMaximum(AnglelimitHigh);
  
  mg2->Draw("ALP3");
  mg2->GetYaxis()->SetTitle("#DeltaA/A (%)");
  mg2->GetYaxis()->CenterTitle();
  mg2->GetXaxis()->SetTitle("Energy (keV)");
  mg2->GetXaxis()->CenterTitle();
  gPad->Modified();

  TLegend *leg2 = new TLegend(0.57,0.7,0.87,0.8);
  leg2->AddEntry(g_delta32_2011,"2011-2012","lf");
  leg2->AddEntry(g_delta32_2012,"2012-2013","lf");
  leg2->SetTextSize(0.05);
  leg2->Draw("SAME");

  //////////////////// DELTA33///////////////////
  TCanvas *c3 = new TCanvas("c3","c3");
  gPad->SetGrid(0,1);

  TMultiGraph *mg3 = new TMultiGraph();
  mg3->SetTitle("#Delta_{3,3} vs. Energy");  
  
  TGraphErrors *g_delta33_2011 = new TGraphErrors(delta33_2011[0].size()-startPoint,&delta33_2011[0][startPoint],&delta33_2011[1][startPoint],0,&delta33_2011[2][startPoint]);
  g_delta33_2011->SetMarkerStyle(0);
  g_delta33_2011->SetLineWidth(3);
  g_delta33_2011->SetLineColor(col2011);
  g_delta33_2011->SetFillColor(col2011);
  g_delta33_2011->SetFillStyle(fill2011);
  
  TGraphErrors *g_delta33_2012 = new TGraphErrors(delta33_2012[0].size()-startPoint,&delta33_2012[0][startPoint],&delta33_2012[1][startPoint],0,&delta33_2012[2][startPoint]);
  g_delta33_2012->SetMarkerStyle(0);
  g_delta33_2012->SetLineWidth(3);
  g_delta33_2012->SetLineColor(col2012);
  g_delta33_2012->SetFillColor(col2012);
  g_delta33_2012->SetFillStyle(fill2012);
  
  mg3->Add(g_delta33_2011);
  mg3->Add(g_delta33_2012);
  mg3->SetMinimum(AnglelimitLow);
  mg3->SetMaximum(AnglelimitHigh);
  
  mg3->Draw("ALP3");
  mg3->GetYaxis()->SetTitle("#DeltaA/A (%)");
  mg3->GetYaxis()->CenterTitle();
  mg3->GetXaxis()->SetTitle("Energy (keV)");
  mg3->GetXaxis()->CenterTitle();
  gPad->Modified();

  TLegend *leg3 = new TLegend(0.57,0.7,0.87,0.8);
  leg3->AddEntry(g_delta33_2011,"2011-2012","lf");
  leg3->AddEntry(g_delta33_2012,"2012-2013","lf");
  leg3->SetTextSize(0.05);
  leg3->Draw("SAME");
  

  //////////////////// DELTA3///////////////////
  TCanvas *cDELTA3 = new TCanvas("cDELTA3","cDELTA3");
  gPad->SetGrid(0,1);

  TMultiGraph *mgDELTA3 = new TMultiGraph();
  mgDELTA3->SetTitle("#Delta_{3} vs. Energy");  
  
  TGraphErrors *g_delta3_2011 = new TGraphErrors(delta3_2011[0].size()-startPoint,&delta3_2011[0][startPoint],&delta3_2011[1][startPoint],0,&delta3_2011[2][startPoint]);
  g_delta3_2011->SetMarkerStyle(0);
  g_delta3_2011->SetLineWidth(3);
  g_delta3_2011->SetLineColor(col2011);
  g_delta3_2011->SetFillColor(col2011);
  g_delta3_2011->SetFillStyle(fill2011);
  
  TGraphErrors *g_delta3_2012 = new TGraphErrors(delta3_2012[0].size()-startPoint,&delta3_2012[0][startPoint],&delta3_2012[1][startPoint],0,&delta3_2012[2][startPoint]);
  g_delta3_2012->SetMarkerStyle(0);
  g_delta3_2012->SetLineWidth(3);
  g_delta3_2012->SetLineColor(col2012);
  g_delta3_2012->SetFillColor(col2012);
  g_delta3_2012->SetFillStyle(fill2012);
  
  mgDELTA3->Add(g_delta3_2011);
  mgDELTA3->Add(g_delta3_2012);
  mgDELTA3->SetMinimum(AnglelimitLow);
  mgDELTA3->SetMaximum(AnglelimitHigh);
  
  mgDELTA3->Draw("ALP3");
  mgDELTA3->GetYaxis()->SetTitle("#DeltaA/A (%)");
  mgDELTA3->GetYaxis()->CenterTitle();
  mgDELTA3->GetXaxis()->SetTitle("Energy (keV)");
  mgDELTA3->GetXaxis()->CenterTitle();
  gPad->Modified();

  TLegend *legDELTA3 = new TLegend(0.57,0.7,0.87,0.8);
  legDELTA3->AddEntry(g_delta3_2011,"2011-2012","lf");
  legDELTA3->AddEntry(g_delta3_2012,"2012-2013","lf");
  legDELTA3->SetTextSize(0.05);
  legDELTA3->Draw("SAME");

  

  TString pdffile = TString::Format("Delta_3_byType_anaChC.pdf");

  c0->Print(TString::Format("%s(",pdffile.Data()));
  c1->Print(pdffile);
  c2->Print(pdffile);
  c3->Print(pdffile);
  cDELTA3->Print(TString::Format("%s)",pdffile.Data()));
  
  
}
