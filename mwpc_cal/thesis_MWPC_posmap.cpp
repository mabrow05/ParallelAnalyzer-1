

#include <string>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <ios>
#include "positionMapHandler.hh"

#include <TH2D.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TString.h>
#include <TPaveText.h>
#include <TMath.h>

using namespace std;

void plot_position_map(int iRunPeriod, double binWidth);


int main(int argc, char* argv[]) {
  int runPeriod = atoi(argv[1]);
  double binWidth = atof(argv[2]);

  plot_position_map(runPeriod, binWidth);
  

}

void plot_position_map(int XePeriod, double binWidth)
{
  cout.setf(ios::fixed, ios::floatfield);
  cout.precision(12);

  // Style options
  //gROOT->SetStyle("Plain");
  Int_t palette[10];
  palette[0] = 0;
  palette[1] = 19;
  palette[2] = 18;
  palette[3] = 17;
  palette[4] = 16;
  palette[5] = 15;
  palette[6] = 14;
  palette[7] = 13;
  palette[8] = 12;
  palette[9] = 1;
  //gStyle->SetPalette(10,palette);  // z-axis color scale for 2D histograms
  gStyle->SetPalette(57);  // z-axis color scale for 2D histograms
  //gStyle->SetOptStat(11);
  gStyle->SetOptStat(0);
  gStyle->SetStatFontSize(0.030);
  gStyle->SetOptFit(1111);
  //gStyle->SetOptTitle(0);
  gStyle->SetTitleSize(0.08,"t");
  //gStyle->SetTitleX(0.17);
  //gStyle->SetTitleAlign(13);
  //gStyle->SetTitleOffset(0.80, "x");
  gStyle->SetTitleOffset(0.9, "y");
  //gStyle->SetPadTickX(1);
  //gStyle->SetPadTickY(1);
  //gStyle->SetNdivisions(510,"X");
  //gStyle->SetNdivisions(510,"Y");
  //gStyle->SetNdivisions(9,"Z");
  gStyle->SetPadLeftMargin(0.12); // 0.13
  gStyle->SetPadRightMargin(0.12); // 0.04
  gStyle->SetPadTopMargin(0.05); // 0.30
  gStyle->SetPadBottomMargin(0.12); // 0.30
  gStyle->SetLabelSize(0.045, "XYZ");
  //gStyle->SetLabelSize(0.045, "Y");
  //gStyle->SetLabelSize(0.045, "Z");
  //gStyle->SetLabelOffset(0.00, "X");
  //gStyle->SetLabelOffset(0.01, "Y");
  gStyle->SetTitleSize(0.050, "X");
  gStyle->SetTitleSize(0.050, "Y");
  //gROOT->ForceStyle();

  gStyle->SetNumberContours(255);

  //This is the position map object which will be used to construct all position values for plotting
  MWPCPositionMap plotMap(0.5,50.);

  int nPosBinsXY = plotMap.getNbinsXY();
  double xyMin = plotMap.getBinLower(0);
  double xyMax = plotMap.getBinUpper(nPosBinsXY-1);  

  cout << nPosBinsXY << " " << xyMin << " " << xyMax << endl;

  double EnergyBinWidth = 50.;
  int nEnergyBins = 14;
  double enBinStart = 100.;
  
  // Histograms
  TH2D *hisE[nEnergyBins];
  TH2D *hisW[nEnergyBins];

  for ( int i=0; i<nEnergyBins; ++i ) {
    
    double elow = enBinStart + i*EnergyBinWidth;
    double ehigh = elow + EnergyBinWidth;

    hisE[i] = new TH2D(TString::Format("East_%.0f-%.0f",elow,ehigh),
		       TString::Format("East %.0f-%.0f keV Reconstructed Energy",elow,ehigh)
		       ,nPosBinsXY,xyMin,xyMax, nPosBinsXY,xyMin,xyMax);
    hisW[i] = new TH2D(TString::Format("West_%.0f-%.0f",elow,ehigh),
		       TString::Format("West %.0f-%.0f keV Reconstructed Energy",elow,ehigh)
		       ,nPosBinsXY,xyMin,xyMax, nPosBinsXY,xyMin,xyMax);
	
  }

  //Fill all the bins

  for ( int i=0; i<nEnergyBins; ++i ) {

    double elow = enBinStart + i*EnergyBinWidth;
    double ehigh = elow + EnergyBinWidth;

    MWPCPositionMap posmap(binWidth, 50.);
    posmap.readMWPCPositionMap(XePeriod,elow,ehigh);
    
    for (int xb=0; xb<nPosBinsXY; xb++) {

      for (int yb=0; yb<nPosBinsXY; yb++) {

	if (TMath::Sqrt(plotMap.getBinCenter(xb)*plotMap.getBinCenter(xb)+plotMap.getBinCenter(yb)*plotMap.getBinCenter(yb)) <= 50.) {
	  
	  std::vector <Double_t> eta = posmap.getInterpolatedEta(plotMap.getBinCenter(xb),plotMap.getBinCenter(yb),
								 plotMap.getBinCenter(xb),plotMap.getBinCenter(yb));
	   
	  hisE[i]->Fill(plotMap.getBinCenter(xb),plotMap.getBinCenter(yb),eta[0]);
	  hisW[i]->Fill(plotMap.getBinCenter(xb),plotMap.getBinCenter(yb),eta[1]);
	  //std::cout << eta[0] << "\t" << eta[1] << std::endl;
	}
      }
    }
  }
  

  // Output PDF file
  TString filenameOut;
  Char_t runPeriodString[10];
  std::cout << "Made it...\n";
  sprintf(runPeriodString,"%i",XePeriod);
  filenameOut = "MWPC_position_map_"+TString(runPeriodString)+TString::Format("_%0.1fmm.pdf",binWidth);

  TString filenameOutFirst;
  filenameOutFirst  = filenameOut;
  filenameOutFirst += "(";

  TString filenameOutLast;
  filenameOutLast  = filenameOut;
  filenameOutLast += ")";

  TCanvas * c0 = new TCanvas("c0", "canvas", 900, 800);
  hisE[4]->Draw("colz");
  hisE[4]->SetXTitle("x [mm]");
  hisE[4]->SetYTitle("y [mm]");
  hisE[4]->GetXaxis()->CenterTitle();
  hisE[4]->GetYaxis()->CenterTitle();
  hisE[4]->GetYaxis()->CenterTitle();
  hisE[4]->SetAxisRange(0.,2.,"Z");
  hisE[4]->SetTitle("");
    
  c0->Print(filenameOutFirst);

  TCanvas * c1 = new TCanvas("c1", "canvas", 900, 800);
  hisW[4]->Draw("colz");
  hisW[4]->SetXTitle("x [mm]");
  hisW[4]->SetYTitle("y [mm]");
  hisW[4]->GetXaxis()->CenterTitle();
  hisW[4]->GetYaxis()->CenterTitle();
  hisW[4]->GetYaxis()->CenterTitle();
  hisW[4]->SetAxisRange(0.,2.,"Z");
  hisW[4]->SetTitle("");

  c1->Print(filenameOutLast);

  /*
  TEllipse *ell = new TEllipse(0,0,45,45);
  ell->SetFillStyle(0);
  ell->SetLineColor(1);
  ell->SetLineStyle(2);
  ell->SetLineWidth(3);
  ell->Draw("same");
  */

  

  delete c0;  delete c1; //delete c3;
  for ( int i=0; i<nEnergyBins; ++i ) { delete hisE[i]; delete hisW[i]; }
  
}
