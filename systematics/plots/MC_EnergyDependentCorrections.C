
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

// Remember to print whichever one you want at the bottom to a pdf
bool color = true;

Int_t groupBin=5;
Double_t enStart=50.;//10.*groupBin;

Double_t BSlimitLow = -1., BSlimitHigh = 3.;
Double_t AnglelimitLow = -5., AnglelimitHigh = 3.;

TString drawOpt = "0AL3";//"0AC4" "0AL3"

std::vector <std::vector<Double_t> > readBSCorr(TString corr,TString year,TString anaCh) {

  std::vector <std::vector<Double_t> > corrs;
  //  std::ifstream infile(TString::Format("../OldMCCorrection/deltaBS%s_anaCh%s_%s.txt",corr.Data(),anaCh.Data(),year.Data()));
  std::ifstream infile(TString::Format("../../systematics/OldMCCorrection/%s_delta_2%s_anaCh%s.txt",year.Data(),corr!=std::string("ALL")?corr.Data():"",anaCh.Data()));
  std::vector<Double_t> enbin;
  std::vector<Double_t> c;
  std::vector<Double_t> cerr;
  Double_t en;
  TString hold1,hold2;
  Int_t i=0;
  Double_t aveEn=0., aveCorr=0.,aveErr=0.;
  while (infile >> en >> hold1 >> hold2 && en<800) {
    if (en<enStart) continue;
    std::cout << en << " " << hold1 << " " << hold2 << "\n";
    if (i<groupBin) {
      aveEn+=en;
      aveCorr+=( hold1==TString("nan")||hold1==TString("-nan"))?0.:atof(hold1.Data());
      aveErr+=(hold1==TString("nan")||hold1==TString("-nan"))?0.:atof(hold2.Data());
      i++;
    } else {
      enbin.push_back(aveEn/groupBin);
      c.push_back(aveCorr/groupBin*100.);
      cerr.push_back(aveErr/groupBin*100.);
      std::cout << aveErr/groupBin*100. << std::endl;
      i=0,aveEn=0.,aveCorr=0.,aveErr=0.;
      aveEn+=en;
      aveCorr+=( hold1==TString("nan")||hold1==TString("-nan"))?0.:atof(hold1.Data());
      aveErr+=(hold1==TString("nan")||hold1==TString("-nan"))?0.:atof(hold2.Data());
      i++;
    }
  }
  enbin.push_back(enbin[enbin.size()-1]+10.*groupBin);
  c.push_back(c[c.size()-1]);
  cerr.push_back(cerr[cerr.size()-1]);
  corrs.push_back(enbin);
  corrs.push_back(c);
  corrs.push_back(cerr);
  return corrs;
};

std::vector <std::vector<Double_t> > readAngleCorr(TString year, TString type,TString anaCh) {

  std::vector <std::vector<Double_t> > corrs;
  std::ifstream infile(TString::Format("../AngleCorrections/%s_delta_3%s_anaCh%s.txt",year.Data(),(type==TString("ALL")?"":type.Data()),anaCh.Data()));
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
  enbin.push_back(enbin[enbin.size()-1]+10.*groupBin);
  c.push_back(c[c.size()-1]);
  cerr.push_back(cerr[cerr.size()-1]);
  corrs.push_back(enbin);
  corrs.push_back(c);
  corrs.push_back(cerr);
  return corrs;
};


void EnergyDependentCorrections(TString anaCh) {

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

  Int_t col2011 = color?9:14;
  Int_t col2012 = color?8:17;
  Int_t fill2011 = 1001;
  Int_t fill2012 = 1001;
  Int_t marker2011 = 20;
  Int_t marker2012 = 21;
  Int_t lineStyle2011 = 1;
  Int_t lineStyle2012 = color?1:2;

  int startPoint=0;
  
  TString year = "2011-2012";
  std::vector<std::vector<Double_t> > bs0_2011 = readBSCorr("0",year,anaCh);
  std::vector<std::vector<Double_t> > bs1_2011 = readBSCorr("1",year,anaCh);
  std::vector<std::vector<Double_t> > bs2_2011 = readBSCorr("2",year,anaCh);
  std::vector<std::vector<Double_t> > bs3_2011 = readBSCorr("3",year,anaCh);
  std::vector<std::vector<Double_t> > bsALL_2011 = readBSCorr("ALL",year,anaCh);
  
  //std::vector<std::vector<Double_t> > delta30_2011 = readAngleCorr(year,"0");
  //std::vector<std::vector<Double_t> > delta31_2011 = readAngleCorr(year,"1");
  //std::vector<std::vector<Double_t> > delta32_2011 = readAngleCorr(year,"2");
  //std::vector<std::vector<Double_t> > delta33_2011 = readAngleCorr(year,"3");
  std::vector<std::vector<Double_t> > delta3_2011 = readAngleCorr(year,"ALL",anaCh);

  year = "2012-2013";
  std::vector<std::vector<Double_t> > bs0_2012 = readBSCorr("0",year,anaCh);
  std::vector<std::vector<Double_t> > bs1_2012 = readBSCorr("1",year,anaCh);
  std::vector<std::vector<Double_t> > bs2_2012 = readBSCorr("2",year,anaCh);
  std::vector<std::vector<Double_t> > bs3_2012 = readBSCorr("3",year,anaCh);
  std::vector<std::vector<Double_t> > bsALL_2012 = readBSCorr("ALL",year,anaCh);
  // std::vector<std::vector<Double_t> > delta30_2012 = readAngleCorr(year,"0");
  //std::vector<std::vector<Double_t> > delta31_2012 = readAngleCorr(year,"1");
  //std::vector<std::vector<Double_t> > delta32_2012 = readAngleCorr(year,"2");
  //std::vector<std::vector<Double_t> > delta33_2012 = readAngleCorr(year,"3");
  std::vector<std::vector<Double_t> > delta3_2012 = readAngleCorr(year,"ALL",anaCh);

  std::vector<std::vector<Double_t> > totalCorr_2011(3,std::vector<Double_t>(0));
  std::vector<std::vector<Double_t> > totalCorr_2012(3,std::vector<Double_t>(0));

  for (UInt_t i=0;i<bsALL_2011[0].size();++i) {
    Double_t c2011 = ((1.+bsALL_2011[1][i]/100.)*(1.+delta3_2011[1][i]/100.)-1.)*100.;
    Double_t c2011err = TMath::Sqrt(TMath::Power((1.+delta3_2011[1][i]/100.)*bsALL_2011[2][i]/100.,2.)
				    +TMath::Power((1.+bsALL_2011[1][i]/100.)*delta3_2011[2][i]/100.,2.))*100.;
    std::cout << bsALL_2011[0][i] << " " << c2011 << " " << c2011err << std::endl;
    totalCorr_2011[0].push_back(bsALL_2011[0][i]);
    totalCorr_2011[1].push_back(c2011);
    totalCorr_2011[2].push_back(c2011err);

    Double_t c2012 = ((1.+bsALL_2012[1][i]/100.)*(1.+delta3_2012[1][i]/100.)-1.)*100.;
    Double_t c2012err = TMath::Sqrt(TMath::Power((1.+delta3_2012[1][i]/100.)*bsALL_2012[2][i]/100.,2.)
				    +TMath::Power((1.+bsALL_2012[1][i]/100.)*delta3_2012[2][i]/100.,2.))*100.;
    totalCorr_2012[0].push_back(bsALL_2012[0][i]);
    totalCorr_2012[1].push_back(c2012);
    totalCorr_2012[2].push_back(c2012err);    
  }
  TCanvas *c0 = new TCanvas("c0","c0");
  gPad->SetGrid(0,1);

  TMultiGraph *mg0 = new TMultiGraph();
  mg0->SetTitle("#Delta_{2,0} vs. Energy");  
  
  TGraphErrors *g_bs0_2011 = new TGraphErrors(bs0_2011[0].size()-startPoint,&bs0_2011[0][startPoint],&bs0_2011[1][startPoint],0,&bs0_2011[2][startPoint]);
  g_bs0_2011->SetMarkerStyle(0);
  g_bs0_2011->SetLineWidth(3);
  g_bs0_2011->SetLineColor(!color?1:col2011);
  g_bs0_2011->SetFillColorAlpha(col2011,0.7);
  g_bs0_2011->SetFillStyle(fill2011);
  g_bs0_2011->SetLineStyle(lineStyle2011);
  

  TGraph *g_bs0_2011_noErr = new TGraph(bs0_2011[0].size()-startPoint,&bs0_2011[0][startPoint],&bs0_2011[1][startPoint]);
 
  TGraphErrors *g_bs0_2012 = new TGraphErrors(bs0_2012[0].size()-startPoint,&bs0_2012[0][startPoint],&bs0_2012[1][startPoint],0,&bs0_2012[2][startPoint]);
  g_bs0_2012->SetMarkerStyle(0);
  g_bs0_2012->SetLineWidth(3);
  g_bs0_2012->SetLineColor(!color?1:col2012);
  g_bs0_2012->SetFillColorAlpha(col2012,0.7);
  g_bs0_2012->SetFillStyle(fill2012);
  g_bs0_2012->SetLineStyle(lineStyle2012);

  mg0->Add(g_bs0_2011);
  mg0->Add(g_bs0_2012);
  mg0->SetMinimum(BSlimitLow);
  mg0->SetMaximum(BSlimitHigh);
  
  mg0->Draw(drawOpt);
  mg0->GetYaxis()->SetTitle("#DeltaA/A (%)");
  mg0->GetYaxis()->CenterTitle();
  mg0->GetXaxis()->SetTitle("Energy (keV)");
  mg0->GetXaxis()->CenterTitle();
  mg0->GetXaxis()->SetLimits(0.,780.);
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
  mg1->SetTitle("#Delta_{2,1} vs. Energy");  
  
  TGraphErrors *g_bs1_2011 = new TGraphErrors(bs1_2011[0].size()-startPoint,&bs1_2011[0][startPoint],&bs1_2011[1][startPoint],0,&bs1_2011[2][startPoint]);
  g_bs1_2011->SetMarkerStyle(0);
  g_bs1_2011->SetLineWidth(3);
  g_bs1_2011->SetLineColor(!color?1:col2011);
  g_bs1_2011->SetFillColorAlpha(col2011,0.7);
  g_bs1_2011->SetFillStyle(fill2011);
  g_bs1_2011->SetLineStyle(lineStyle2011);

  
  TGraphErrors *g_bs1_2012 = new TGraphErrors(bs1_2012[0].size()-startPoint,&bs1_2012[0][startPoint],&bs1_2012[1][startPoint],0,&bs1_2012[2][startPoint]);
  g_bs1_2012->SetMarkerStyle(0);
  g_bs1_2012->SetLineWidth(3);
  g_bs1_2012->SetLineColor(!color?1:col2012);
  g_bs1_2012->SetFillColorAlpha(col2012,0.7);
  g_bs1_2012->SetFillStyle(fill2012);
  g_bs1_2012->SetLineStyle(lineStyle2012);
  
  mg1->Add(g_bs1_2011);
  mg1->Add(g_bs1_2012);
  mg1->SetMinimum(BSlimitLow);
  mg1->SetMaximum(BSlimitHigh);
  
  mg1->Draw(drawOpt);
  mg1->GetYaxis()->SetTitle("#DeltaA/A (%)");
  mg1->GetYaxis()->CenterTitle();
  mg1->GetXaxis()->SetTitle("Energy (keV)");
  mg1->GetXaxis()->CenterTitle();
  mg1->GetXaxis()->SetLimits(0.,780.);

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
  mg2->SetTitle("#Delta_{2,2} vs. Energy");  
  
  TGraphErrors *g_bs2_2011 = new TGraphErrors(bs2_2011[0].size()-startPoint,&bs2_2011[0][startPoint],&bs2_2011[1][startPoint],0,&bs2_2011[2][startPoint]);
  g_bs2_2011->SetMarkerStyle(0);
  g_bs2_2011->SetLineWidth(3);
  g_bs2_2011->SetLineColor(!color?1:col2011);
  g_bs2_2011->SetFillColorAlpha(col2011,0.7);
  g_bs2_2011->SetFillStyle(fill2011);
  g_bs2_2011->SetLineStyle(lineStyle2011);

  
  TGraphErrors *g_bs2_2012 = new TGraphErrors(bs2_2012[0].size()-startPoint,&bs2_2012[0][startPoint],&bs2_2012[1][startPoint],0,&bs2_2012[2][startPoint]);
  g_bs2_2012->SetMarkerStyle(0);
  g_bs2_2012->SetLineWidth(3);
  g_bs2_2012->SetLineColor(!color?1:col2012);
  g_bs2_2012->SetFillColorAlpha(col2012,0.7);
  g_bs2_2012->SetFillStyle(fill2012);
  g_bs2_2012->SetLineStyle(lineStyle2012);

  
  mg2->Add(g_bs2_2011);
  mg2->Add(g_bs2_2012);
  mg2->SetMinimum(BSlimitLow);
  mg2->SetMaximum(BSlimitHigh);
  
  mg2->Draw(drawOpt);
  mg2->GetYaxis()->SetTitle("#DeltaA/A (%)");
  mg2->GetYaxis()->CenterTitle();
  mg2->GetXaxis()->SetTitle("Energy (keV)");
  mg2->GetXaxis()->CenterTitle();
  mg2->GetXaxis()->SetLimits(0.,780.);
    
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
  mg3->SetTitle("#Delta_{2,3} vs. Energy");  
  
  TGraphErrors *g_bs3_2011 = new TGraphErrors(bs3_2011[0].size()-startPoint,&bs3_2011[0][startPoint],&bs3_2011[1][startPoint],0,&bs3_2011[2][startPoint]);
  g_bs3_2011->SetMarkerStyle(0);
  g_bs3_2011->SetLineWidth(3);
  g_bs3_2011->SetLineColor(!color?1:col2011);
  g_bs3_2011->SetFillColorAlpha(col2011,0.7);
  g_bs3_2011->SetFillStyle(fill2011);
  g_bs3_2011->SetLineStyle(lineStyle2011);

  TGraphErrors *g_bs3_2012 = new TGraphErrors(bs3_2012[0].size()-startPoint,&bs3_2012[0][startPoint],&bs3_2012[1][startPoint],0,&bs3_2012[2][startPoint]);
  g_bs3_2012->SetMarkerStyle(0);
  g_bs3_2012->SetLineWidth(3);
  g_bs3_2012->SetLineColor(!color?1:col2012);
  g_bs3_2012->SetFillColorAlpha(col2012,0.7);
  g_bs3_2012->SetFillStyle(fill2012);
  g_bs3_2012->SetLineStyle(lineStyle2012);

  mg3->Add(g_bs3_2011);
  mg3->Add(g_bs3_2012);
  mg3->SetMinimum(BSlimitLow);
  mg3->SetMaximum(BSlimitHigh);
  
  mg3->Draw(drawOpt);
  mg3->GetYaxis()->SetTitle("#DeltaA/A (%)");
  mg3->GetYaxis()->CenterTitle();
  mg3->GetXaxis()->SetTitle("Energy (keV)");
  mg3->GetXaxis()->CenterTitle();
  mg3->GetXaxis()->SetLimits(0.,780.);

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
  mgALLBS->SetTitle("#Delta_{2} vs. Energy");  
  
  TGraphErrors *g_bsALL_2011 = new TGraphErrors(bsALL_2011[0].size()-startPoint,&bsALL_2011[0][startPoint],&bsALL_2011[1][startPoint],0,&bsALL_2011[2][startPoint]);
  g_bsALL_2011->SetMarkerStyle(0);
  g_bsALL_2011->SetLineWidth(3);
  g_bsALL_2011->SetLineColor(!color?1:col2011);
  g_bsALL_2011->SetFillColorAlpha(col2011,0.7);
  g_bsALL_2011->SetFillStyle(fill2011);
  g_bsALL_2011->SetLineStyle(lineStyle2011);

  
  TGraphErrors *g_bsALL_2012 = new TGraphErrors(bsALL_2012[0].size()-startPoint,&bsALL_2012[0][startPoint],&bsALL_2012[1][startPoint],0,&bsALL_2012[2][startPoint]);
  g_bsALL_2012->SetMarkerStyle(0);
  g_bsALL_2012->SetLineWidth(3);
  g_bsALL_2012->SetLineColor(!color?1:col2012);
  g_bsALL_2012->SetFillColorAlpha(col2012,0.7);
  g_bsALL_2012->SetFillStyle(fill2012);
  g_bsALL_2012->SetLineStyle(lineStyle2012);
    
  mgALLBS->Add(g_bsALL_2011);
  mgALLBS->Add(g_bsALL_2012);
  mgALLBS->SetMinimum(BSlimitLow);
  mgALLBS->SetMaximum(BSlimitHigh);
  
  mgALLBS->Draw(drawOpt);
  mgALLBS->GetYaxis()->SetTitle("#DeltaA/A (%)");
  mgALLBS->GetYaxis()->CenterTitle();
  mgALLBS->GetXaxis()->SetTitle("Energy (keV)");
  mgALLBS->GetXaxis()->CenterTitle();
  mgALLBS->GetXaxis()->SetLimits(0.,780.);

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
  mgAngle->SetTitle("#Delta_{3} vs. Energy");  
  
  TGraphErrors *g_delta3_2011 = new TGraphErrors(delta3_2011[0].size()-startPoint,&delta3_2011[0][startPoint],&delta3_2011[1][startPoint],0,&delta3_2011[2][startPoint]);
  g_delta3_2011->SetMarkerStyle(0);
  g_delta3_2011->SetLineWidth(3);
  g_delta3_2011->SetLineColor(!color?1:col2011);
  g_delta3_2011->SetFillColorAlpha(col2011,0.7);
  g_delta3_2011->SetFillStyle(fill2011);
  g_delta3_2011->SetLineStyle(lineStyle2011);

  
  TGraphErrors *g_delta3_2012 = new TGraphErrors(delta3_2012[0].size()-startPoint,&delta3_2012[0][startPoint],&delta3_2012[1][startPoint],0,&delta3_2012[2][startPoint]);
  g_delta3_2012->SetMarkerStyle(0);
  g_delta3_2012->SetLineWidth(3);
  g_delta3_2012->SetLineColor(!color?1:col2012);
  g_delta3_2012->SetFillColorAlpha(col2012,0.7);
  g_delta3_2012->SetFillStyle(fill2012);
  g_delta3_2012->SetLineStyle(lineStyle2012);

  mgAngle->Add(g_delta3_2011);
  mgAngle->Add(g_delta3_2012);
  mgAngle->SetMinimum(AnglelimitLow);
  mgAngle->SetMaximum(AnglelimitHigh);
  
  mgAngle->Draw(drawOpt);
  mgAngle->GetYaxis()->SetTitle("#DeltaA/A (%)");
  mgAngle->GetYaxis()->CenterTitle();
  mgAngle->GetXaxis()->SetTitle("Energy (keV)");
  mgAngle->GetXaxis()->CenterTitle();
  mgAngle->GetXaxis()->SetLimits(0.,780.);

  gPad->Modified();

  TLegend *legAngle = new TLegend(0.57,0.7,0.87,0.8);
  legAngle->AddEntry(g_delta3_2011,"2011-2012","lf");
  legAngle->AddEntry(g_delta3_2012,"2012-2013","lf");
  legAngle->SetTextSize(0.05);
  legAngle->Draw("SAME");

  //////////////////// TotalCorr ///////////////////
  TCanvas *cTotalCorr = new TCanvas("cTotalCorr","cTotalCorr");
  gPad->SetGrid(0,1);

  TMultiGraph *mgTotalCorr = new TMultiGraph();
  mgTotalCorr->SetTitle("#Delta_{MC} vs. Energy");  
  
  TGraphErrors *g_totalCorr_2011 = new TGraphErrors(totalCorr_2011[0].size()-startPoint,&totalCorr_2011[0][startPoint],&totalCorr_2011[1][startPoint],0,&totalCorr_2011[2][startPoint]);
  g_totalCorr_2011->SetMarkerStyle(0);
  g_totalCorr_2011->SetLineWidth(3);
  g_totalCorr_2011->SetLineColor(!color?1:col2011);
  g_totalCorr_2011->SetFillColorAlpha(col2011,0.7);
  g_totalCorr_2011->SetFillStyle(fill2011);
  g_totalCorr_2011->SetLineStyle(lineStyle2011);

  
  TGraphErrors *g_totalCorr_2012 = new TGraphErrors(totalCorr_2012[0].size()-startPoint,&totalCorr_2012[0][startPoint],&totalCorr_2012[1][startPoint],0,&totalCorr_2012[2][startPoint]);
  g_totalCorr_2012->SetMarkerStyle(0);
  g_totalCorr_2012->SetLineWidth(3);
  g_totalCorr_2012->SetLineColor(!color?1:col2012);
  g_totalCorr_2012->SetFillColorAlpha(col2012,0.7);
  g_totalCorr_2012->SetFillStyle(fill2012);
  g_totalCorr_2012->SetLineStyle(lineStyle2012);

  mgTotalCorr->Add(g_totalCorr_2011);
  mgTotalCorr->Add(g_totalCorr_2012);
  mgTotalCorr->SetMinimum(AnglelimitLow);
  mgTotalCorr->SetMaximum(AnglelimitHigh);
  
  mgTotalCorr->Draw(drawOpt);
  mgTotalCorr->GetYaxis()->SetTitle("#DeltaA/A (%)");
  mgTotalCorr->GetYaxis()->CenterTitle();
  mgTotalCorr->GetXaxis()->SetTitle("Energy (keV)");
  mgTotalCorr->GetXaxis()->CenterTitle();
  mgTotalCorr->GetXaxis()->SetLimits(0.,780.);

  gPad->Modified();

  TLegend *legTotalCorr = new TLegend(0.57,0.7,0.87,0.8);
  legTotalCorr->AddEntry(g_totalCorr_2011,"2011-2012","lf");
  legTotalCorr->AddEntry(g_totalCorr_2012,"2012-2013","lf");
  legTotalCorr->SetTextSize(0.05);
  legTotalCorr->Draw("SAME");
  

  TString pdffile = TString::Format("Delta_2_byType_anaCh%s_%iBinAve%s.pdf",anaCh.Data(),groupBin,(color?"_color":""));

  c0->Print(TString::Format("%s(",pdffile.Data()));
  c1->Print(pdffile);
  c2->Print(pdffile);
  c3->Print(pdffile);
  cALLBS->Print(pdffile);
  cAngle->Print(pdffile);
  cTotalCorr->Print(TString::Format("%s)",pdffile.Data()));
  
}
