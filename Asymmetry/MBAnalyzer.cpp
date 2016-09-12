/*
Compilable code which houses the heart of the analysis. 
Here are classes which handle calculating asymmetries of either data 
or simulation. There will be flags to be passed dictating what type of 
asymmetry, what type of data, and what octet/runs to use. 

The code should be able to do BG subtraction, construct asymmetries, 
or do both depending on the flags passed.

Maybe even add in writing the final answer to the database if the user wants to

*/

#include "SQLinterface.hh"
#include "EvtRateHandler.hh"
#include "Asymmetries.hh"
#include "SystematicCorrections.hh"
#include "MBUtils.hh"
#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <cmath>
#include <fstream>
#include <algorithm>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <TStyle.h>


bool BLINDED = true;

std::vector <Int_t> badOct = {7,9,59}; // Octet 7 had W anode dead for part of run
                                  // Either need to discard, or apply the 
                                  // Charge cloud method to determine a 
                                  // coincidence

// The Process functions will calculate all the raw asymmetries on a bin-by-bin basis and write them to file
void ProcessOctets(Int_t octBegin, Int_t octEnd, Int_t anaChoice=1, Double_t enBinWidth=10., bool UKdata=true, bool simulation=false, bool AsymmOn=true);
void ProcessQuartets(Int_t octBegin, Int_t octEnd, Int_t anaChoice=1, Double_t enBinWidth=10., bool UKdata=true, bool simulation=false, bool AsymmOn=true);
void ProcessPairs(Int_t octBegin, Int_t octEnd, Int_t anaChoice=1, Double_t enBinWidth=10., bool UKdata=true, bool simulation=false, bool AsymmOn=true); 

// makes plots of 2*A/Beta for each octet (pair, quartet, and octet as a whole) and fits over the range specified, for whatever grouping provided ("Octet", "Quartet", "Pair")
void PlotAsymmetriesByGrouping(std::string groupType, Int_t octBegin, Int_t octEnd, Int_t anaChoice, Double_t Elow=220., Double_t Ehigh=680., Double_t enBinWidth=10., bool UKdata=true, bool simulation=false);

// Collects all the asymmetries, bin-by-bin, and produces a final asymmetry plot, both A_SR and 2*A/Beta, for whatever grouping provided ("Octet", "Quartet", "Pair"). Also makes a plot of the 
// integrated asymmetry vs octet/quartet/pair number
void PlotFinalAsymmetries(std::string groupType, Int_t octBegin, Int_t octEnd, Int_t anaChoice, Double_t Elow=220., Double_t Ehigh=680., Double_t enBinWidth=10., bool UKdata=true, bool simulation=false, bool AsymmOn=true);

// Function to return beta when given the kinetic energy of an electron
Double_t returnBeta(Double_t En) { 
  Double_t me = 510.998928; //rest mass energy of electron in keV
  return sqrt(En*En+2.*En*me)/(En+me);
};



int main()
{

  Int_t analysisChoice = 1;
  Int_t octBegin = 0;
  Int_t octEnd = 59;
  Double_t enBinWidth = 10.;
  Double_t Elow = 180.;
  Double_t Ehigh = 780.;
  bool UKdata = true;
  bool simulation = false;
  bool applyAsymm = true;
  
  try {
    
    ProcessOctets(octBegin, octEnd, analysisChoice, enBinWidth, UKdata, simulation, applyAsymm);
    PlotAsymmetriesByGrouping("Octet",octBegin, octEnd, analysisChoice, Elow, Ehigh, enBinWidth, UKdata, simulation);
    PlotFinalAsymmetries("Octet",octBegin, octEnd, analysisChoice, Elow, Ehigh, enBinWidth, UKdata, simulation, applyAsymm);
    
    ProcessQuartets(octBegin, octEnd, analysisChoice, enBinWidth, UKdata, simulation, applyAsymm);
    PlotAsymmetriesByGrouping("Quartet",octBegin, octEnd, analysisChoice, Elow, Ehigh, enBinWidth, UKdata, simulation);
    PlotFinalAsymmetries("Quartet",octBegin, octEnd, analysisChoice, Elow, Ehigh, enBinWidth, UKdata, simulation);
     
    ProcessPairs(octBegin, octEnd, analysisChoice, enBinWidth, UKdata, simulation, applyAsymm);
    PlotAsymmetriesByGrouping("Pair",octBegin, octEnd, analysisChoice, Elow, Ehigh, enBinWidth, UKdata, simulation);
    PlotFinalAsymmetries("Pair",octBegin, octEnd, analysisChoice, Elow, Ehigh, enBinWidth, UKdata, simulation);
    
  }
  catch(const char* ex){
    std::cerr << "Error: " << ex << std::endl;
  }
  
  
  return 0;
}


void ProcessOctets(Int_t octBegin, Int_t octEnd, Int_t anaChoice, Double_t enBinWidth, bool UKdata, bool simulation, bool AsymmOn) {

  for (Int_t octet=octBegin; octet<=octEnd; octet++) {
    try {
      OctetAsymmetry oct(octet,enBinWidth, 50., UKdata, simulation, AsymmOn);
      oct.calcAsymmetryBinByBin(anaChoice);     
      oct.calcSuperSum(anaChoice);
      oct.writeAsymToFile();
      oct.writeSuperSumToFile();
    }
    catch(const char* ex){
      std::cerr << "Error: " << ex << std::endl;
    }
  }
};

void ProcessQuartets(Int_t octBegin, Int_t octEnd, Int_t anaChoice, Double_t enBinWidth, bool UKdata, bool simulation, bool AsymmOn) {
  
  for (Int_t octet=octBegin; octet<=octEnd; octet++) {
    try {
      QuartetAsymmetry quart(octet,enBinWidth, 50., UKdata, simulation, AsymmOn);
      quart.calcAsymmetryBinByBin(anaChoice);
      quart.calcSuperSum(anaChoice);
      quart.writeAsymToFile();
      quart.writeSuperSumToFile();
    }
    catch(const char* ex){
      std::cerr << "Error: " << ex << std::endl;
    }
  }
  
};

void ProcessPairs(Int_t octBegin, Int_t octEnd, Int_t anaChoice, Double_t enBinWidth, bool UKdata, bool simulation, bool AsymmOn) {
  
  for (Int_t octet=octBegin; octet<=octEnd; octet++) {
    try {
      PairAsymmetry pair(octet,enBinWidth, 50., UKdata, simulation, AsymmOn);
      pair.calcAsymmetryBinByBin(anaChoice);
      pair.calcSuperSum(anaChoice);
      pair.writeAsymToFile();
      pair.writeSuperSumToFile();
    }
    catch(const char* ex){
      std::cerr << "Error: " << ex << std::endl;
    }
  }
  
};


void PlotFinalAsymmetries(std::string groupType, Int_t octBegin, Int_t octEnd, Int_t anaChoice, Double_t Elow, Double_t Ehigh, Double_t enBinWidth, bool UKdata, bool simulation, bool AsymmOn) {

  if (groupType!="Quartet" && groupType!="Octet" && groupType!="Pair") throw "Bad group type given to PlotFinalAsymmetries. Options are \"Octet\", \"Quartet\", or \"Pair\""; 

  Int_t numBins = 1200/enBinWidth;

  std::vector < Double_t > enBinMedian; //Holds the center of the energy bins
  for (Int_t i=0; i<numBins; i++) {
    Double_t En = i*enBinWidth+enBinWidth/2.;
    enBinMedian.push_back(En);
  }

  std::vector < std::vector <Double_t> > rawAsymByGroup(3,std::vector <Double_t> (0));
  std::vector < std::vector <Double_t> > AsymByGroup(3,std::vector <Double_t> (0));
  
  std::vector < std::vector <Double_t> > groupRawAsymByBin(3,std::vector <Double_t> (numBins,0.));
  std::vector < std::vector <Double_t> > groupAsymByBin(3,std::vector <Double_t> (numBins,0.));

  
  //bool isOctet = groupType=="Octet" ? true : false;
  bool isQuartet = groupType=="Quartet" ? true : false;
  bool isPair = groupType=="Pair" ? true : false;
  std::map < std::string, Int_t > numQuartetLoops = {{"Octet",1},{"Quartet",2},{"Pair",2}};
  std::map < std::string, Int_t > numPairLoops = {{"Octet",1},{"Quartet",1},{"Pair",2}};
  std::string quartetName[2] = {"A","B"};
  
  std::string basePath = simulation ? getenv("SIM_ANALYSIS_RESULTS") : UKdata ? getenv("ANALYSIS_RESULTS") : getenv("MPM_ANALYSIS_RESULTS");
  
  std::string txt;
  Double_t binEdge, Asym, AsymError;

  Int_t pair = 0, quartet = 0;

  for (Int_t octet=octBegin; octet<=octEnd; octet++) {

    if (std::find(badOct.begin(), badOct.end(),octet) != badOct.end()) { pair+=4; quartet+=2; continue; } //Checking if octet should be ignored for data quality reasons
    
    for (Int_t q = 0; q<numQuartetLoops[groupType]; q++) {

      for (Int_t p=0; p<numPairLoops[groupType]; p++) {

	std::string path = basePath+"Octet_"+itos(octet)+"/"+groupType+"Asymmetry/"+"FittedAsymmetry_Octet"+itos(octet)+"_AnaCh"+itos(anaChoice)+"_"+
	  (isQuartet?("Quartet_"+quartetName[q]+"_"):isPair?("Pair_"+quartetName[q]+itos(p)+"_"):"")+itos((int)Elow)+"-"+itos((int)Ehigh) + ".dat";
	//std::cout << path << std::endl;

	ifstream infile(path.c_str());  
   
	if (infile.is_open()) {
	  infile >> txt >> Asym >> AsymError;
	  rawAsymByGroup[0].push_back(isQuartet?quartet:isPair?pair:octet);
	  rawAsymByGroup[1].push_back(Asym);
	  rawAsymByGroup[2].push_back(AsymError);

	  //std::cout << octet << " " << Asym << " " << AsymError << std::endl;
      
	  infile >> txt >> Asym >> AsymError;
	  AsymByGroup[0].push_back(isQuartet?quartet:isPair?pair:octet);
	  AsymByGroup[1].push_back(Asym);
	  AsymByGroup[2].push_back(AsymError);
	  infile.close();
	
	  path = basePath + "Octet_" + itos(octet) + "/" + groupType + "Asymmetry/" + "rawAsymmetry_Octet" + itos(octet) + "_AnaCh" + itos(anaChoice)  +
	    (isQuartet?("_Quartet_"+quartetName[q]):isPair?("_Pair_"+quartetName[q]+itos(p)):"")+".dat"; 
	  infile.open(path.c_str());
	  //std::cout << path << std::endl;
	  
	  Int_t i = 0;
	  while (infile >> binEdge >> Asym >> AsymError) {
	    groupRawAsymByBin[0][i] = binEdge;
	    groupRawAsymByBin[1][i] += AsymError>0. ? 1./power(AsymError,2)*Asym : 0.;
	    groupRawAsymByBin[2][i] += AsymError>0. ? 1/power(AsymError,2) : 0.;
	    std::cout << binEdge << " " << groupRawAsymByBin[1][i] << " " << groupRawAsymByBin[2][i] << std::endl;
	    i++;
	  }
	  infile.close();
	}
	pair++;
      }
      quartet++;
    }
  }
  //Do final calculations of the rates in each bin and their associated errors
  for (unsigned int i=0; i<groupRawAsymByBin[1].size(); i++) {
    groupRawAsymByBin[1][i] = groupRawAsymByBin[2][i]>0. ? -groupRawAsymByBin[1][i] / groupRawAsymByBin[2][i] : 0.; // This finishes the calculation of the weighted average and errors 
    groupRawAsymByBin[2][i] = groupRawAsymByBin[2][i]>0. ? 1./sqrt(groupRawAsymByBin[2][i]) : 0.;
    
    groupAsymByBin[0][i] = groupRawAsymByBin[0][i];
    groupAsymByBin[1][i] = 2*groupRawAsymByBin[1][i]/returnBeta(enBinMedian[i]);
    groupAsymByBin[2][i] = 2*groupRawAsymByBin[2][i]/returnBeta(enBinMedian[i]);
    
    //std::cout << enBinMedian[i] << " " << groupAsymByBin[1][i] << " " << groupAsymByBin[2][i] << std::endl;
  }
  
  // Plotting stuff
  std::string pdfFile = basePath + "Asymmetries/"+groupType+"RawAsymmetries_AnaCh" + itos(anaChoice) 
    + std::string("_") + itos((int)Elow) + std::string("-") + itos((int)Ehigh) +std::string(".pdf");
  
  TCanvas *c1 = new TCanvas("c1", "c1", 1000., 1400.);
  c1->Divide(1,4);
  
  gStyle->SetOptFit(1111);
  gStyle->SetTitleX(0.25);
  gStyle->SetStatX(0.75);
  gStyle->SetStatY(0.85);
  
  std::vector <double> errorX(1000,0.);
  std::string title;
  
  c1->cd(1);
  
  TGraphErrors *g = new TGraphErrors(rawAsymByGroup[0].size(), &rawAsymByGroup[0][0],&rawAsymByGroup[1][0],&errorX[0], &rawAsymByGroup[2][0]);
  title = groupType+"-By-"+groupType+" Raw Measured Asymmetry " + itos((int)Elow) + std::string("-") +itos((int)Ehigh) + std::string(" keV Window");
  g->SetTitle(title.c_str());
  g->SetMarkerStyle(20);
  g->SetLineWidth(2);
  g->GetXaxis()->SetLimits(-2., rawAsymByGroup[0][rawAsymByGroup[0].size()-1]+2.);
  g->GetXaxis()->SetTitle("Number");
  g->GetYaxis()->SetTitle("Raw Asymmetry");
  g->GetXaxis()->CenterTitle();
  g->GetYaxis()->CenterTitle();
  
  TF1 *fit = new TF1("fit","[0]",rawAsymByGroup[0][0]-1., rawAsymByGroup[0][rawAsymByGroup[0].size()-1]+1.);
  fit->SetLineColor(kRed);
  fit->SetLineWidth(3);
  fit->SetParameter(0,0.05);
  
  g->Fit("fit","R");
  
  g->Draw("AP");
  g->SetMinimum((simulation && !AsymmOn) ? -0.05 : 0.03);
  g->SetMaximum((simulation && !AsymmOn) ? 0.05 : 0.07);
  c1->Update();
  
  c1->cd(2);
  
  TGraphErrors *gBeta = new TGraphErrors(AsymByGroup[0].size(), &AsymByGroup[0][0],&AsymByGroup[1][0],&errorX[0], &AsymByGroup[2][0]);
  title = groupType+"-By-"+groupType+" 2A_{SR}/^{}#beta " + itos((int)Elow) + std::string("-") +itos((int)Ehigh) + std::string(" keV Window");
  gBeta->SetTitle(title.c_str());
  gBeta->SetMarkerStyle(20);
  gBeta->SetLineWidth(2);
  gBeta->GetXaxis()->SetLimits(-2., rawAsymByGroup[0][rawAsymByGroup[0].size()-1]+2.);
  gBeta->GetXaxis()->SetTitle("Number");
  gBeta->GetYaxis()->SetTitle("Asymmetry");
  gBeta->GetXaxis()->CenterTitle();
  gBeta->GetYaxis()->CenterTitle();
  
  TF1 *fitBeta = new TF1("fitBeta","[0]",rawAsymByGroup[0][0]-1., rawAsymByGroup[0][rawAsymByGroup[0].size()-1]+1.);
  fitBeta->SetLineColor(kRed);
  fitBeta->SetLineWidth(3);
  fitBeta->SetParameter(0,0.05);
  
  gBeta->Fit("fitBeta","R");
  
  gBeta->Draw("AP");
  gBeta->SetMinimum((simulation && !AsymmOn) ? -0.5 : -0.14);
  gBeta->SetMaximum((simulation && !AsymmOn) ? 0.5 : -0.09);
  c1->Update();
  
  
  //c1->Print(pdfFile.c_str());

  c1->cd(3);
  Int_t offset = 0;
  TGraphErrors *g2 = new TGraphErrors(enBinMedian.size()-offset, &enBinMedian[offset], &groupRawAsymByBin[1][offset], 0, &groupRawAsymByBin[2][offset]);
  title = groupType + std::string(" Raw Asymmetry A_{SR} ") + itos((int)Elow) + std::string("-") +itos((int)Ehigh) + std::string(" keV Window");
  g2->SetTitle(title.c_str());
  g2->SetMarkerStyle(20);
  g2->SetLineWidth(2);
  g2->GetXaxis()->SetLimits(0., 800.);
  g2->GetXaxis()->SetTitle("Energy (keV)");
  g2->GetYaxis()->SetTitle("Uncorrected Asymmetry");
  g2->GetXaxis()->CenterTitle();
  g2->GetYaxis()->CenterTitle();
  
  g2->Draw("AP");
  g2->SetMinimum((simulation && !AsymmOn) ? -0.05 : -0.08);
  g2->SetMaximum((simulation && !AsymmOn) ? 0.05 : -0.01);
  c1->Update();
	
  c1->cd(4);
	
  TGraphErrors *gBeta2 = new TGraphErrors(enBinMedian.size()-offset, &enBinMedian[offset], &groupAsymByBin[1][offset], 0, &groupAsymByBin[2][offset]);
  title = groupType + std::string(" 2A_{SR}/^{}#beta ") + itos((int)Elow) + std::string("-") +itos((int)Ehigh) + std::string(" keV Window");
  gBeta2->SetTitle(title.c_str());
  gBeta2->SetMarkerStyle(20);
  gBeta2->SetLineWidth(2);
  gBeta2->GetXaxis()->SetLimits(0., 800.);
  gBeta2->GetXaxis()->SetTitle("Energy (keV)");
  gBeta2->GetYaxis()->SetTitle("Asymmetry");
  gBeta2->GetXaxis()->CenterTitle();
  gBeta2->GetYaxis()->CenterTitle();
  
  TF1 *fitBeta2 = new TF1("fitBeta2","[0]",Elow, Ehigh);
  fitBeta2->SetLineColor(kRed);
  fitBeta2->SetLineWidth(3);
  fitBeta2->SetParameter(0,-0.12);
  
  gBeta2->Fit("fitBeta2","R");
  
  gBeta2->Draw("AP");
  gBeta2->SetMinimum((simulation && !AsymmOn) ? -0.5 : -0.16);
  gBeta2->SetMaximum((simulation && !AsymmOn) ? 0.5 : -0.07);
  c1->Update();
  c1->Print(pdfFile.c_str());
  
  
  delete c1; delete g; delete fit; delete gBeta; delete fitBeta; 
  delete g2; delete gBeta2; delete fitBeta2;

};


//This will create asymmetry plots for each pair, quartet, and octet, and also write the raw integrated asymmetry and Beta corrected integrated asymmetry
// to file for each grouping to be used by PlotFinalAsymmetries

void PlotAsymmetriesByGrouping(std::string groupType, Int_t octBegin, Int_t octEnd, Int_t anaChoice, Double_t Elow, Double_t Ehigh, Double_t enBinWidth, bool UKdata, bool simulation) {
  if (groupType!="Quartet" && groupType!="Octet" && groupType!="Pair") throw "Bad group type given to PlotAsymmetriesByGrouping. Options are \"Octet\", \"Quartet\", or \"Pair\""; 

  std::string basePath = simulation ? getenv("SIM_ANALYSIS_RESULTS") : UKdata ? getenv("ANALYSIS_RESULTS") : getenv("MPM_ANALYSIS_RESULTS");

  for (Int_t octet=octBegin; octet<=octEnd; octet++) {

    ifstream infile;

    if (groupType=="Octet")
    {
      std::string basePath2 =basePath+ "Octet_"+itos(octet)+"/" + groupType + "Asymmetry/"; 
      std::string infilePath = basePath2 + "rawAsymmetry_Octet" + itos(octet) + "_AnaCh" + itos(anaChoice) + ".dat";
      std::string outfilePath = basePath2 + "FittedAsymmetry_Octet" + itos(octet) + "_AnaCh" + itos(anaChoice) + "_" + itos((int)Elow) + std::string("-") +itos((int)Ehigh) + ".dat";
      
      //Remove old output files in case they were created on accident and aren't filled with good values
      std::string command = "rm " + outfilePath; 
      system(command.c_str());

      //First check that Octet was good
      std::string checkStatus;
      infile.open(infilePath.c_str());

      if ( infile.is_open() ) {
	infile >> checkStatus; 
	infile.close();
	if (checkStatus=="BAD") continue;
      }
      else {
	std::cout << "Could not open binned Asymmetries for Octet " << octet << " so skipping...\n";
	continue; 
      }

      std::vector < std::vector <Double_t > > AsymAndError;
      std::vector < std::vector <Double_t > > RawAsymAndError;
      std::vector < Double_t > enBinMedian;
	
      AsymAndError.resize(2,std::vector <Double_t> (0));
      RawAsymAndError.resize(2,std::vector <Double_t> (0));
	
      Double_t eBinLow, Asym, AsymError;
	
      infile.open(infilePath.c_str());
	
      while (infile >> eBinLow >> Asym >> AsymError) {
	Double_t En = eBinLow+enBinWidth/2.;
	enBinMedian.push_back(En);
	Double_t Beta = returnBeta(En);
	AsymAndError[0].push_back(-Asym*2./Beta);
	AsymAndError[1].push_back(AsymError*2./Beta);
	RawAsymAndError[0].push_back(-Asym);
	RawAsymAndError[1].push_back(AsymError);
      }
      
      infile.close();
	
      std::string pdfPath = basePath2 + "Asymmetry_Octet" + itos(octet) + "_AnaCh" + itos(anaChoice) + "_" + itos((int)Elow) + std::string("-") +itos((int)Ehigh) + ".pdf";
      TCanvas *c1 = new TCanvas("c1", "c1",800, 400);
      gStyle->SetOptFit(1111);
      gStyle->SetTitleX(0.25);
	
      TGraphErrors *gOct = new TGraphErrors(enBinMedian.size()-5, &enBinMedian[5], &RawAsymAndError[0][5], 0, &RawAsymAndError[1][5]);
      std::string title = std::string("Raw Asymmetry A_{SR} ") + itos((int)Elow) + std::string("-") +itos((int)Ehigh) + std::string(" keV Window");
      gOct->SetTitle(title.c_str());
      gOct->SetMarkerStyle(20);
      gOct->SetLineWidth(2);
      gOct->GetXaxis()->SetLimits(0., 800.);
      gOct->GetXaxis()->SetTitle("Energy (keV)");
      gOct->GetYaxis()->SetTitle("Uncorrected Asymmetry");
      gOct->GetXaxis()->CenterTitle();
      gOct->GetYaxis()->CenterTitle();
	
      TF1 *fitOct = new TF1("fitOct","[0]",Elow, Ehigh);
      fitOct->SetLineColor(kRed);
      fitOct->SetLineWidth(3);
      fitOct->SetParameter(0,-0.05);
	
      gOct->Fit("fitOct","R");
	
      ofstream ofile(outfilePath.c_str());
      ofile << "RawA_SR " << -(fitOct->GetParameter(0)) << " " << fitOct->GetParError(0) << std::endl;
	
      gOct->Draw("AP");
      gOct->SetMinimum(-0.1);
      gOct->SetMaximum(0.);
      c1->Update();
      c1->Print(pdfPath.c_str());
	
	
      delete c1; delete gOct; delete fitOct;
	
      pdfPath = basePath2 + "BetaCorrectedAsymmetry_Octet" + itos(octet) + "_AnaCh" + itos(anaChoice) + "_" + itos((int)Elow) + std::string("-") +itos((int)Ehigh) + ".pdf";
      c1 = new TCanvas("c1", "c1",800, 400);
	
      gOct = new TGraphErrors(enBinMedian.size()-5, &enBinMedian[5], &AsymAndError[0][5], 0, &AsymAndError[1][5]);
      title = std::string("#frac{2A_{SR}}{#beta} ") + itos((int)Elow) + std::string("-") +itos((int)Ehigh) + std::string(" keV Window");
      gOct->SetTitle(title.c_str());
      gOct->SetMarkerStyle(20);
      gOct->SetLineWidth(2);
      gOct->GetXaxis()->SetLimits(0., 800.);
      gOct->GetXaxis()->SetTitle("Energy (keV)");
      gOct->GetYaxis()->SetTitle("Uncorrected Asymmetry");
      gOct->GetXaxis()->CenterTitle();
      gOct->GetYaxis()->CenterTitle();
	
      fitOct = new TF1("fitOct","[0]",Elow, Ehigh);
      fitOct->SetLineColor(kRed);
      fitOct->SetLineWidth(3);
      fitOct->SetParameter(0,-0.12);
	
      gOct->Fit("fitOct","R");
	
      gOct->Draw("AP");
      gOct->SetMinimum(-0.6);
      gOct->SetMaximum(0.4);
      c1->Update();
      c1->Print(pdfPath.c_str());
	
	
      ofile << "2*A/Beta " << fitOct->GetParameter(0) << " " << fitOct->GetParError(0);
      ofile.close();
	
      delete c1; delete gOct; delete fitOct;
      
    }

    //Now the Quartets
    else if (groupType=="Quartet")
    {	  

      std::string basePath2 =basePath+ "Octet_"+itos(octet)+"/" + groupType + "Asymmetry/"; 

      std::string infilePath[2];
      infilePath[0] = basePath2 + "rawAsymmetry_Octet" + itos(octet) + "_AnaCh" + itos(anaChoice) + "_Quartet_A" + ".dat";
      infilePath[1] = basePath2 + "rawAsymmetry_Octet" + itos(octet) + "_AnaCh" + itos(anaChoice) + "_Quartet_B" + ".dat";

      std::string outfilePath[2];
      outfilePath[0] = basePath2 + "FittedAsymmetry_Octet" + itos(octet) + "_AnaCh" + itos(anaChoice) + "_Quartet_A" + "_" + itos((int)Elow) + std::string("-") +itos((int)Ehigh) +".dat";
      outfilePath[1] = basePath2 + "FittedAsymmetry_Octet" + itos(octet) + "_AnaCh" + itos(anaChoice) + "_Quartet_B" + "_" + itos((int)Elow) + std::string("-") +itos((int)Ehigh) +".dat";

      std::string pdfPathRaw[2];
      pdfPathRaw[0] = basePath2 + "Asymmetry_Octet" + itos(octet) + "_AnaCh" + itos(anaChoice) + "_Quartet_A" + "_" + itos((int)Elow) + std::string("-") +itos((int)Ehigh) + ".pdf";
      pdfPathRaw[1] = basePath2 + "Asymmetry_Octet" + itos(octet) + "_AnaCh" + itos(anaChoice) + "_Quartet_B" + "_" + itos((int)Elow) + std::string("-") +itos((int)Ehigh) + ".pdf";

      std::string pdfPathCorr[2];
      pdfPathCorr[0] = basePath2 + "BetaCorrectedAsymmetry_Octet" + itos(octet) + "_AnaCh" + itos(anaChoice) + "_Quartet_A" + "_" + itos((int)Elow) + std::string("-") +itos((int)Ehigh) + ".pdf";
      pdfPathCorr[1] = basePath2 + "BetaCorrectedAsymmetry_Octet" + itos(octet) + "_AnaCh" + itos(anaChoice) + "_Quartet_B" + "_" + itos((int)Elow) + std::string("-") +itos((int)Ehigh) + ".pdf";
      
      for (int quart=0; quart<2; quart++) {

	std::string currentQuart = quart==0?"A":"B";

	//Remove old output files in case they were created on accident and aren't filled with good values
	std::string command = "rm " + outfilePath[quart]; 
	system(command.c_str());

	//First check that Quartet was good
	std::string checkStatus;
	infile.open(infilePath[quart].c_str());
	
	if ( infile.is_open() ) {
	  infile >> checkStatus; 
	  infile.close();
	  if (checkStatus=="BAD") continue;
	}
	else {
	  std::cout << "Could not open binned Asymmetries for Quartet " << currentQuart << " in Octet " << octet << " so skipping...\n";
	  continue; 
	}
	    
	infile.open(infilePath[quart].c_str());
      
	std::vector < std::vector <Double_t > > AsymAndError;
	std::vector < std::vector <Double_t > > RawAsymAndError;
	std::vector < Double_t > enBinMedian;
	
	AsymAndError.resize(2,std::vector <Double_t> (0));
	RawAsymAndError.resize(2,std::vector <Double_t> (0));
	
	Double_t eBinLow, Asym, AsymError;
	
	while (infile >> eBinLow >> Asym >> AsymError) {
	  Double_t En = eBinLow+enBinWidth/2.;
	  enBinMedian.push_back(En);
	  Double_t Beta = returnBeta(En);
	  AsymAndError[0].push_back(-Asym*2./Beta);
	  AsymAndError[1].push_back(AsymError*2./Beta);
	  RawAsymAndError[0].push_back(-Asym);
	  RawAsymAndError[1].push_back(AsymError);
	}

	infile.close();
	
	TCanvas *c1 = new TCanvas("c1", "c1",800, 400);
	gStyle->SetOptFit(1111);
	gStyle->SetTitleX(0.25);
      
	TGraphErrors *gOct = new TGraphErrors(enBinMedian.size()-5, &enBinMedian[5], &RawAsymAndError[0][5], 0, &RawAsymAndError[1][5]);
	std::string title = std::string("Raw Asymmetry A_{SR} ") + itos((int)Elow) + std::string("-") +itos((int)Ehigh) + std::string(" keV Window");
	gOct->SetTitle(title.c_str());
	gOct->SetMarkerStyle(20);
	gOct->SetLineWidth(2);
	gOct->GetXaxis()->SetLimits(0., 800.);
	gOct->GetXaxis()->SetTitle("Energy (keV)");
	gOct->GetYaxis()->SetTitle("Uncorrected Asymmetry");
	gOct->GetXaxis()->CenterTitle();
	gOct->GetYaxis()->CenterTitle();

	TF1 *fitOct = new TF1("fitOct","[0]",Elow, Ehigh);
	fitOct->SetLineColor(kRed);
	fitOct->SetLineWidth(3);
	fitOct->SetParameter(0,-0.05);
	
	gOct->Fit("fitOct","R");
	
	ofstream ofile(outfilePath[quart].c_str());
	ofile << "RawA_SR " << -(fitOct->GetParameter(0)) << " " << fitOct->GetParError(0) << std::endl;
      
	gOct->Draw("AP");
	gOct->SetMinimum(-0.1);
	gOct->SetMaximum(0.);
	c1->Update();
	c1->Print(pdfPathRaw[quart].c_str());
	
	delete c1; delete gOct; delete fitOct;
	
	c1 = new TCanvas("c1", "c1",800, 400);
	
	gOct = new TGraphErrors(enBinMedian.size()-5, &enBinMedian[5], &AsymAndError[0][5], 0, &AsymAndError[1][5]);
	title = std::string("#frac{2A_{SR}}{#beta} ") + itos((int)Elow) + std::string("-") +itos((int)Ehigh) + std::string(" keV Window");
	gOct->SetTitle(title.c_str());
	gOct->SetMarkerStyle(20);
	gOct->SetLineWidth(2);
	gOct->GetXaxis()->SetLimits(0., 800.);
	gOct->GetXaxis()->SetTitle("Energy (keV)");
	gOct->GetYaxis()->SetTitle("Uncorrected Asymmetry");
	gOct->GetXaxis()->CenterTitle();
	gOct->GetYaxis()->CenterTitle();
	
        fitOct = new TF1("fitOct","[0]",Elow, Ehigh);
	fitOct->SetLineColor(kRed);
	fitOct->SetLineWidth(3);
	fitOct->SetParameter(0,-0.12);
	
	gOct->Fit("fitOct","R");
	
	gOct->Draw("AP");
	gOct->SetMinimum(-0.6);
	gOct->SetMaximum(0.4);
	c1->Update();
	c1->Print(pdfPathCorr[quart].c_str());

	ofile << "2*A/Beta " << fitOct->GetParameter(0) << " " << fitOct->GetParError(0);
	ofile.close();

	delete c1; delete gOct; delete fitOct;
      }
    }

    //And last of all the Pairs
    else if (groupType=="Pair")
    {
      std::string basePath2 =basePath+ "Octet_"+itos(octet)+"/" + groupType + "Asymmetry/"; 

      for (int pair=0; pair<2; pair++) {
	std::string infilePath[2];
	infilePath[0] = basePath2 + "rawAsymmetry_Octet" + itos(octet) + "_AnaCh" + itos(anaChoice) + "_Pair_A" + itos(pair) + ".dat";
	infilePath[1] = basePath2 + "rawAsymmetry_Octet" + itos(octet) + "_AnaCh" + itos(anaChoice) + "_Pair_B" + itos(pair) + ".dat";

	std::string outfilePath[2];
	outfilePath[0] = basePath2 + "FittedAsymmetry_Octet" + itos(octet) + "_AnaCh" + itos(anaChoice) + "_Pair_A" + itos(pair) + "_" + itos((int)Elow) + std::string("-") +itos((int)Ehigh) + ".dat";
	outfilePath[1] = basePath2 + "FittedAsymmetry_Octet" + itos(octet) + "_AnaCh" + itos(anaChoice) + "_Pair_B" + itos(pair) + "_" + itos((int)Elow) + std::string("-") +itos((int)Ehigh) + ".dat";

	std::string pdfPathRaw[2];
	pdfPathRaw[0] = basePath2 + "Asymmetry_Octet" + itos(octet) + "_AnaCh" + itos(anaChoice) + "_Pair_A" + itos(pair) + "_" + itos((int)Elow) + std::string("-") +itos((int)Ehigh) + ".pdf";
	pdfPathRaw[1] = basePath2 + "Asymmetry_Octet" + itos(octet) + "_AnaCh" + itos(anaChoice) + "_Pair_B" + itos(pair) + "_" + itos((int)Elow) + std::string("-") +itos((int)Ehigh) + ".pdf";
	
	std::string pdfPathCorr[2];
	pdfPathCorr[0] = basePath2 + "BetaCorrectedAsymmetry_Octet" + itos(octet) + "_AnaCh" + itos(anaChoice) + "_Pair_A" + itos(pair) + "_" + itos((int)Elow) + std::string("-") +itos((int)Ehigh) + ".pdf";
	pdfPathCorr[1] = basePath2 + "BetaCorrectedAsymmetry_Octet" + itos(octet) + "_AnaCh" + itos(anaChoice) + "_Pair_B" + itos(pair) + "_" + itos((int)Elow) + std::string("-") +itos((int)Ehigh) + ".pdf";
	

	//Remove old output files in case they were created on accident and aren't filled with good values
	std::string command = "rm " + outfilePath[0]; 
	system(command.c_str());
	command = "rm " + outfilePath[1]; 
	system(command.c_str());
	
	for (int quart=0; quart<2; quart++) {
	  
	  std::string currentPair = (quart==0?"A":"B") + itos(pair);	  

	  //First check that pair was good
	  std::string checkStatus;
	  infile.open(infilePath[quart].c_str());
	  
	  if ( infile.is_open() ) {
	    infile >> checkStatus; 
	    infile.close();
	    if (checkStatus=="BAD") continue;
	  }
	  else {
	    std::cout << "Could not open binned Asymmetries for Pair " << currentPair << " in Octet " << octet << " so skipping...\n";
	    continue; 
	  }


	 
	  
	  infile.open(infilePath[quart].c_str());
	  
	  std::vector < std::vector <Double_t > > AsymAndError;
	  std::vector < std::vector <Double_t > > RawAsymAndError;
	  std::vector < Double_t > enBinMedian;
	  
	  AsymAndError.resize(2,std::vector <Double_t> (0));
	  RawAsymAndError.resize(2,std::vector <Double_t> (0));
	  
	  Double_t eBinLow, Asym, AsymError;
	
	  while (infile >> eBinLow >> Asym >> AsymError) {
	    Double_t En = eBinLow+enBinWidth/2.;
	    enBinMedian.push_back(En);
	    Double_t Beta = returnBeta(En);
	    AsymAndError[0].push_back(-Asym*2./Beta);
	    AsymAndError[1].push_back(AsymError*2./Beta);
	    RawAsymAndError[0].push_back(-Asym);
	    RawAsymAndError[1].push_back(AsymError);
	  }

	  infile.close();
	  
	  
	  TCanvas *c1 = new TCanvas("c1", "c1",800, 400);
	  gStyle->SetOptFit(1111);
	  gStyle->SetTitleX(0.25);
      
	  TGraphErrors *gOct = new TGraphErrors(enBinMedian.size()-5, &enBinMedian[5], &RawAsymAndError[0][5], 0, &RawAsymAndError[1][5]);
	  std::string title = std::string("Raw Asymmetry A_{SR} ") + itos((int)Elow) + std::string("-") +itos((int)Ehigh) + std::string(" keV Window");
	  gOct->SetTitle(title.c_str());
	  gOct->SetMarkerStyle(20);
	  gOct->SetLineWidth(2);
	  gOct->GetXaxis()->SetLimits(0., 800.);
	  gOct->GetXaxis()->SetTitle("Energy (keV)");
	  gOct->GetYaxis()->SetTitle("Uncorrected Asymmetry");
	  gOct->GetXaxis()->CenterTitle();
	  gOct->GetYaxis()->CenterTitle();

	  TF1 *fitOct = new TF1("fitOct","[0]",Elow, Ehigh);
	  fitOct->SetLineColor(kRed);
	  fitOct->SetLineWidth(3);
	  fitOct->SetParameter(0,-0.05);
	  
	  gOct->Fit("fitOct","R");
	  
	  ofstream ofile(outfilePath[quart].c_str());
	  ofile << "RawA_SR " << -(fitOct->GetParameter(0)) << " " << fitOct->GetParError(0) << std::endl;
	  
	  gOct->Draw("AP");
	  gOct->SetMinimum(-0.1);
	  gOct->SetMaximum(0.);
	  c1->Update();
	  c1->Print(pdfPathRaw[quart].c_str());
	  
	  delete c1; delete gOct; delete fitOct;
	  
	  c1 = new TCanvas("c1", "c1",800, 400);
	  
	  gOct = new TGraphErrors(enBinMedian.size()-5, &enBinMedian[5], &AsymAndError[0][5], 0, &AsymAndError[1][5]);
	  title = std::string("#frac{2A_{SR}}{#beta} ") + itos((int)Elow) + std::string("-") +itos((int)Ehigh) + std::string(" keV Window");
	  gOct->SetTitle(title.c_str());
	  gOct->SetMarkerStyle(20);
	  gOct->SetLineWidth(2);
	  gOct->GetXaxis()->SetLimits(0., 800.);
	  gOct->GetXaxis()->SetTitle("Energy (keV)");
	  gOct->GetYaxis()->SetTitle("Uncorrected Asymmetry");
	  gOct->GetXaxis()->CenterTitle();
	  gOct->GetYaxis()->CenterTitle();
	  
	  fitOct = new TF1("fitOct","[0]",Elow, Ehigh);
	  fitOct->SetLineColor(kRed);
	  fitOct->SetLineWidth(3);
	  fitOct->SetParameter(0,-0.12);
	  
	  gOct->Fit("fitOct","R");
	  
	  gOct->Draw("AP");
	  gOct->SetMinimum(-0.6);
	  gOct->SetMaximum(0.4);
	  c1->Update();
	  c1->Print(pdfPathCorr[quart].c_str());

	  ofile << "2*A/Beta " << fitOct->GetParameter(0) << " " << fitOct->GetParError(0);
	  ofile.close();
	  
	  delete c1; delete gOct; delete fitOct;
	}
      }
    }
  }
};
      
  

//These are summed over the energy range given as Elow-Ehigh

void ProduceRawAsymmetries(Int_t octBegin, Int_t octEnd, Int_t anaChoice, Double_t Elow, Double_t Ehigh, Double_t enBinWidth, bool UKdata, bool simulation, bool AsymmOn) {
 
    //Looking at octet evts and asymmetries
  std::string basePath = simulation ? getenv("SIM_ANALYSIS_RESULTS") : UKdata ? getenv("ANALYSIS_RESULTS") : getenv("MPM_ANALYSIS_RESULTS");
  
  std::string pairFile = basePath + std::string("Asymmetries/Pair_RawAsymmetries_AnaCh") + itos(anaChoice) 
    + std::string("_") + itos((int)Elow) + std::string("-") + itos((int)Ehigh) +std::string(".dat");
  std::string quartetFile = basePath + std::string("Asymmetries/Quartet_RawAsymmetries_AnaCh") + itos(anaChoice)
    + std::string("_") + itos((int)Elow) + std::string("-") + itos((int)Ehigh) +std::string(".dat");
  std::string octetFile = basePath + std::string("Asymmetries/Octet_RawAsymmetries_AnaCh") + itos(anaChoice) 
    + std::string("_") + itos((int)Elow) + std::string("-") + itos((int)Ehigh) +std::string(".dat");

    
    //unsigned int octetNum = 1;
  ofstream octAsym(octetFile.c_str());
  ofstream quartAsym(quartetFile.c_str());
  ofstream pairAsym(pairFile.c_str());
 
  Int_t numPair=0, numQuart=0;
  std::vector < std::vector <Double_t > > octetRawAsymAndError;
  std::vector < std::vector <Double_t > > quartetRawAsymAndError;
  std::vector < std::vector <Double_t > > pairRawAsymAndError;
  octetRawAsymAndError.resize(3,std::vector <Double_t> (0));
  quartetRawAsymAndError.resize(3,std::vector <Double_t> (0));
  pairRawAsymAndError.resize(3,std::vector <Double_t> (0)); 
  
  Double_t Asym=0., AsymError=0.;

  for (Int_t octet=octBegin; octet<=octEnd; octet++) {

    if (std::find(badOct.begin(), badOct.end(),octet) != badOct.end()) { numPair+=4; numQuart+=2; continue; } //Checking if octet should be ignored for data quality reasons

    try {
      //OCTETS
      OctetAsymmetry oct(octet,enBinWidth, 50., UKdata, simulation, AsymmOn);     
      oct.calcTotalAsymmetry(Elow,Ehigh,anaChoice);
      //oct.calcAsymmetryBinByBin(1);
      //oct.calcTotalAsymmetry(170.,630.,1);
      
      Asym = oct.returnTotalAsymmetry();
      AsymError = oct.returnTotalAsymmetryError();
      octetRawAsymAndError[0].push_back(octet);
      octetRawAsymAndError[1].push_back(Asym);
      octetRawAsymAndError[2].push_back(AsymError);

      std::cout << "Asymmetry for Octet " << octet << ":\n";
      std::cout << Asym << " +- " << AsymError << std::endl;
      octAsym << octet << " " << Asym << " " << AsymError << "\n";
     
    }
    catch(const char* ex){
      std::cerr << "Error: " << ex << std::endl;
    }
    
    try {
      // QUARTETS
      QuartetAsymmetry quart(octet,enBinWidth, 50., UKdata, simulation, AsymmOn); 
      quart.calcTotalAsymmetry(Elow,Ehigh,anaChoice);
	//quart.calcAsymmetryBinByBin(1);
	//quart.calcTotalAsymmetry(170.,630.,1);
      
      if (quart.boolGoodQuartet(0)) {
	
	Asym = quart.returnTotalAsymmetry_QuartetA();
	AsymError = quart.returnTotalAsymmetryError_QuartetA();
	quartetRawAsymAndError[0].push_back(numQuart);
	quartetRawAsymAndError[1].push_back(Asym);
	quartetRawAsymAndError[2].push_back(AsymError);
  
	std::cout << "Asymmetry for Quartet " << numQuart << ":\n";
	std::cout << Asym << " +- " << AsymError << std::endl;
	quartAsym << numQuart << " " << Asym << " " << AsymError << "\n";
      }
      numQuart++;

      if (quart.boolGoodQuartet(1)) {
	Asym = quart.returnTotalAsymmetry_QuartetB();
	AsymError = quart.returnTotalAsymmetryError_QuartetB();
	quartetRawAsymAndError[0].push_back(numQuart);
	quartetRawAsymAndError[1].push_back(Asym);
	quartetRawAsymAndError[2].push_back(AsymError);

	std::cout << "Asymmetry for Quartet " << numQuart << ":\n";
	std::cout << Asym << " +- " << AsymError << std::endl;
	quartAsym << numQuart << " " << Asym << " " << AsymError << "\n";
      }
      numQuart++;
    }
    catch(const char* ex){
      numQuart+=2;
      std::cerr << "Error: " << ex << std::endl;
    }

    try {
      // PAIRS
      PairAsymmetry pair(octet,enBinWidth, 50., UKdata, simulation, AsymmOn); 
      pair.calcTotalAsymmetry(Elow,Ehigh,anaChoice);
	//pair.calcAsymmetryBinByBin(1);
	//pair.calcTotalAsymmetry(170.,630.,1);
      
      if (pair.boolGoodPair(0,0)) {
	
	Asym = pair.returnTotalAsymmetry_PairA0();
	AsymError = pair.returnTotalAsymmetryError_PairA0();
	pairRawAsymAndError[0].push_back(numPair);
	pairRawAsymAndError[1].push_back(Asym);
	pairRawAsymAndError[2].push_back(AsymError);
	
	std::cout << "Asymmetry for Pair " << numPair << ":\n";
	std::cout << Asym << " +- " << AsymError << std::endl;
	pairAsym << numPair << " " << Asym << " " << AsymError << "\n";
      }
      numPair++;

      if (pair.boolGoodPair(0,1)) {
	Asym = pair.returnTotalAsymmetry_PairA1();
	AsymError = pair.returnTotalAsymmetryError_PairA1();
	pairRawAsymAndError[0].push_back(numPair);
	pairRawAsymAndError[1].push_back(Asym);
	pairRawAsymAndError[2].push_back(AsymError);

	std::cout << "Asymmetry for Pair " << numPair << ":\n";
	std::cout << Asym << " +- " << AsymError << std::endl;
	pairAsym << numPair << " " << Asym << " " << AsymError << "\n";
      }
      numPair++;

      if (pair.boolGoodPair(1,0)) {
	
	Asym = pair.returnTotalAsymmetry_PairB0();
	AsymError = pair.returnTotalAsymmetryError_PairB0();
	pairRawAsymAndError[0].push_back(numPair);
	pairRawAsymAndError[1].push_back(Asym);
	pairRawAsymAndError[2].push_back(AsymError);

	std::cout << "Asymmetry for Pair " << numPair << ":\n";
	std::cout << Asym << " +- " << AsymError << std::endl;
	pairAsym << numPair << " " << Asym << " " << AsymError << "\n";
      }
      numPair++;

      if (pair.boolGoodPair(1,1)) {
	Asym = pair.returnTotalAsymmetry_PairB1();
	AsymError = pair.returnTotalAsymmetryError_PairB1();
	pairRawAsymAndError[0].push_back(numPair);
	pairRawAsymAndError[1].push_back(Asym);
	pairRawAsymAndError[2].push_back(AsymError);

	std::cout << "Asymmetry for Pair " << numPair << ":\n";
	std::cout << Asym << " +- " << AsymError << std::endl;
	pairAsym << numPair << " " << Asym << " " << AsymError << "\n";
      }
      numPair++;
    }
    catch(const char* ex){
      numPair+=4;
      std::cerr << "Error: " << ex << std::endl;
    }
  }

  octAsym.close();
  quartAsym.close();
  pairAsym.close();

  std::string pdfFile = basePath + std::string("Asymmetries/RawAsymmetries_AnaCh") + itos(anaChoice) 
    + std::string("_") + itos((int)Elow) + std::string("-") + itos((int)Ehigh) +std::string(".pdf");

  TCanvas *c1 = new TCanvas("c1", "c1", 1000., 1000.);
  c1->Divide(1,3);

  gStyle->SetOptFit(1111);
  gStyle->SetTitleX(0.25);
 
  std::vector <double> errorX(pairRawAsymAndError[0].size(),0.);
  std::string title;

  c1->cd(1);

  TGraphErrors *gOct = new TGraphErrors(octetRawAsymAndError[0].size(), &octetRawAsymAndError[0][0],&octetRawAsymAndError[1][0],&errorX[0], &octetRawAsymAndError[2][0]);
  title = std::string("Octet-By-Octet Raw Measured Asymmetry ") + itos((int)Elow) + std::string("-") +itos((int)Ehigh) + std::string(" keV Window");
  gOct->SetTitle(title.c_str());
  gOct->SetMarkerStyle(20);
  gOct->SetLineWidth(2);
  gOct->GetXaxis()->SetLimits(-2., octetRawAsymAndError[0][octetRawAsymAndError[0].size()-1]+2.);
  gOct->GetXaxis()->SetTitle("Octet Number");
  gOct->GetYaxis()->SetTitle("Raw Asymmetry");
  gOct->GetXaxis()->CenterTitle();
  gOct->GetYaxis()->CenterTitle();
  
  TF1 *fitOct = new TF1("fitOct","[0]",octetRawAsymAndError[0][0]-1., octetRawAsymAndError[0][octetRawAsymAndError[0].size()-1]+1.);
  fitOct->SetLineColor(kRed);
  fitOct->SetLineWidth(3);
  fitOct->SetParameter(0,0.05);

  gOct->Fit("fitOct","R");
  
  gOct->Draw("AP");
  gOct->SetMinimum(0.03);
  gOct->SetMaximum(0.07);
  c1->Update();

  c1->cd(2);

  TGraphErrors *gQuart = new TGraphErrors(quartetRawAsymAndError[0].size(), &quartetRawAsymAndError[0][0],&quartetRawAsymAndError[1][0],&errorX[0], &quartetRawAsymAndError[2][0]);
  title = std::string("Quartet-By-Quartet Raw Measured Asymmetry ") + itos((int)Elow) + std::string("-") +itos((int)Ehigh) + std::string(" keV Window");
  gQuart->SetTitle(title.c_str());
  gQuart->SetMarkerStyle(20);
  gQuart->SetLineWidth(2);
  gQuart->GetXaxis()->SetLimits(-2., quartetRawAsymAndError[0][quartetRawAsymAndError[0].size()-1]+2.);
  gQuart->GetXaxis()->SetTitle("Quartet Number");
  gQuart->GetYaxis()->SetTitle("Raw Asymmetry");
  gQuart->GetXaxis()->CenterTitle();
  gQuart->GetYaxis()->CenterTitle();
  
  TF1 *fitQuart = new TF1("fitQuart","[0]",quartetRawAsymAndError[0][0]-1., quartetRawAsymAndError[0][quartetRawAsymAndError[0].size()-1]+1.);
  fitQuart->SetLineColor(kRed);
  fitQuart->SetLineWidth(3);
  fitQuart->SetParameter(0,0.05);


  gQuart->Fit("fitQuart","R");
  
  gQuart->Draw("AP");
  gQuart->SetMinimum(0.03);
  gQuart->SetMaximum(0.07);
  c1->Update();

  c1->cd(3);

  TGraphErrors *gPair = new TGraphErrors(pairRawAsymAndError[0].size(), &pairRawAsymAndError[0][0],&pairRawAsymAndError[1][0],&errorX[0], &pairRawAsymAndError[2][0]);
  title = std::string("Pair-By-Pair Raw Measured Asymmetry ") + itos((int)Elow) + std::string("-") +itos((int)Ehigh) + std::string(" keV Window");
  gPair->SetTitle(title.c_str());
  gPair->SetMarkerStyle(20);
  gPair->SetLineWidth(2);
  gPair->GetXaxis()->SetLimits(-2., pairRawAsymAndError[0][pairRawAsymAndError[0].size()-1]+2.);
  gPair->GetXaxis()->SetTitle("Pair Number");
  gPair->GetYaxis()->SetTitle("Raw Asymmetry");
  gPair->GetXaxis()->CenterTitle();
  gPair->GetYaxis()->CenterTitle();
  
  TF1 *fitPair = new TF1("fitPair","[0]",pairRawAsymAndError[0][0]-1., pairRawAsymAndError[0][pairRawAsymAndError[0].size()-1]+1.);
  fitPair->SetLineColor(kRed);
  fitPair->SetLineWidth(3);
  fitPair->SetParameter(0,0.05);


  gPair->Fit("fitPair","R");
  
  gPair->Draw("AP");
  gPair->SetMinimum(0.03);
  gPair->SetMaximum(0.07);
  c1->Update();
  
  c1->Print(pdfFile.c_str());

  delete gPair; delete gQuart; delete gOct; delete fitPair; delete fitQuart; delete fitOct; delete c1;
																     

};
  
  
