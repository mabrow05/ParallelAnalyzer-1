/*
 Compares the BG subtracted data spectra to simulation data which has 
 detector effects already taken into account
*/

#include "Asymmetries.hh"
#include "MBUtils.hh"
#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>

#include <TMath.h>
#include <TString.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TH1F.h>
#include <TH1I.h>
#include <TGraphErrors.h>
#include <TPDF.h>
#include <TLegend.h>
#include <TPaveText.h>

int main(int argc, char *argv[])
{
  if (argc!=6) {
    std::cout << "Usage: ./UK2SIMcomparison [octet start] [octet end] [energy window low] [energy window high] [energy bin width]\n";
    std::cout << "The code will produce comparisons for every octet in the range given, on an octet-by-octet basis\n";
    exit(0);
  }

  gStyle->SetOptStat(0);

  int octetNumStart = atoi(argv[1]);
  int octetNumEnd = atoi(argv[2]);
  double enWinLow = (double) atof(argv[3]);
  double enWinHigh = (double) atof(argv[4]);
  int enBinWidth = atoi(argv[5]);

  int numBins = 1000/enBinWidth;

  //set up energy binning vector
  std::vector <double> enBins(numBins,0.);
  
  for (int bin=0; bin<numBins; bin++) {
    enBins[bin] = ((double)bin+.5)*(double)enBinWidth;
  }
  
  std::vector <std::vector <double> > superSumTotal_uk(4, std::vector <double>(numBins,0.)); //type0=0, type1=1, type23=2, ALL=3
  std::vector <std::vector <double> > superSumTotalError_uk(4, std::vector <double>(numBins,0.));
  std::vector <std::vector <double> > superSumTotal_sim(4, std::vector <double>(numBins,0.));
  std::vector <std::vector <double> > superSumTotalError_sim(4, std::vector <double>(numBins,0.));
  std::vector <double> xErr(numBins,0.);

  for (int octetNum=octetNumStart; octetNum<octetNumEnd+1; octetNum++) {
 
    try {
      //std::string pdfFileBase = "comparison_plots/SimulationComp_Octet"+itos(octetNum)+"_evis_"+itos((int)enWinLow)+"-"+itos((int)enWinHigh)+".pdf";
      std::string file1 = "superSumPlots/octet_"+itos(octetNum)+".root";
      TFile *f = new TFile(file1.c_str(),"RECREATE");
      
      

      //int uk0E=0, sim0E=0, uk1E=0, sim1E=0, uk23E=0, sim23E=0, ukTotE=0, simTotE=0, uk0W=0, sim0W=0, uk1W=0, sim1W=0, uk23W=0, sim23W=0, ukTotW=0, simTotW=0; //Event totals
      //double ukAsym = 0., simAsym=0., ukAsymError = 0., simAsymError=0.;
      
      double normFactor = 0., ukIntegral=0., simIntegral=0.;


      { //Setting scope for first OctetAsymmetry so that it will be deleted to clear memory
	OctetAsymmetry UK(octetNum, enBinWidth, 50., true, false);

	//std::vector <double> evts;

	std::vector <double> superSum;
	std::vector <double> superSumError;
	

	f->cd();
	
	//Type 0 events
	UK.calcSuperSum(4);
	superSum = UK.returnSuperSum();
	superSumError = UK.returnSuperSumError();

	for (int n=0; n<numBins;n++) {
	  //std::cout << enBins[n] << " " << superSum[n] << " " << superSumError[n] << std::endl; 
	  superSumTotal_uk[0][n]+=superSumError[n]>0.?(1./power(superSumError[n],2))*superSum[n]:0.;
	  superSumTotalError_uk[0][n]+=superSumError[n]>0.?(1./power(superSumError[n],2)):0.;
	  if (enBins[n]>=enWinLow && enBins[n]<=enWinHigh) ukIntegral+=superSum[n];
	}
	
	TGraphErrors *Erecon0_uk = new TGraphErrors(numBins,&enBins[0],&superSum[0],&xErr[0],&superSumError[0]);
	Erecon0_uk->SetMarkerColor(kBlue);
	Erecon0_uk->SetMarkerStyle(21);
	Erecon0_uk->GetXaxis()->SetTitle("Energy (keV)");
	Erecon0_uk->SetTitle("Type 0 Super-Sum");
	Erecon0_uk->Draw("AP");
	
	Erecon0_uk->Write("Erecon0_uk");
	
	delete Erecon0_uk;

	//Type 1 events
	UK.calcSuperSum(6);
	superSum = UK.returnSuperSum();
	superSumError = UK.returnSuperSumError();

	for (int n=0; n<numBins;n++) {
	  superSumTotal_uk[1][n]+=superSumError[n]>0.?(1./power(superSumError[n],2))*superSum[n]:0.;
	  superSumTotalError_uk[1][n]+=superSumError[n]>0.?(1./power(superSumError[n],2)):0.;
	}
	
	TGraphErrors *Erecon1_uk = new TGraphErrors(numBins,&enBins[0],&superSum[0],&xErr[0],&superSumError[0]);
	Erecon1_uk->SetMarkerColor(kBlue);
	Erecon1_uk->SetMarkerStyle(21);
	Erecon1_uk->GetXaxis()->SetTitle("Energy (keV)");
	Erecon1_uk->SetTitle("Type 1 Super-Sum");
	Erecon1_uk->Draw("AP");
	
	Erecon1_uk->Write("Erecon1_uk");
	
	delete Erecon1_uk;

	//Type 2/3 events
	UK.calcSuperSum(7);
	superSum = UK.returnSuperSum();
	superSumError = UK.returnSuperSumError();

	for (int n=0; n<numBins;n++) {
	  superSumTotal_uk[2][n]+=superSumError[n]>0.?(1./power(superSumError[n],2))*superSum[n]:0.;
	  superSumTotalError_uk[2][n]+=superSumError[n]>0.?(1./power(superSumError[n],2)):0.;
	}
	
	TGraphErrors *Erecon23_uk = new TGraphErrors(numBins,&enBins[0],&superSum[0],&xErr[0],&superSumError[0]);
	Erecon23_uk->SetMarkerColor(kBlue);
	Erecon23_uk->SetMarkerStyle(21);
	Erecon23_uk->GetXaxis()->SetTitle("Energy (keV)");
	Erecon23_uk->SetTitle("All Types Super-Sum");
	Erecon23_uk->Draw("AP");
	
	Erecon23_uk->Write("Erecon23_uk");
	
	delete Erecon23_uk;

	//All event types
	UK.calcSuperSum(1);
	superSum = UK.returnSuperSum();
	superSumError = UK.returnSuperSumError();

	for (int n=0; n<numBins;n++) {
	  superSumTotal_uk[3][n]+=superSumError[n]>0.?(1./power(superSumError[n],2))*superSum[n]:0.;
	  superSumTotalError_uk[3][n]+=superSumError[n]>0.?(1./power(superSumError[n],2)):0.;
	}
	
	TGraphErrors *EreconALL_uk = new TGraphErrors(numBins,&enBins[0],&superSum[0],&xErr[0],&superSumError[0]);
	EreconALL_uk->SetMarkerColor(kBlue);
	EreconALL_uk->SetMarkerStyle(21);
	EreconALL_uk->GetXaxis()->SetTitle("Energy (keV)");
	EreconALL_uk->SetTitle("All Types Super-Sum");
	EreconALL_uk->Draw("AP");
	
	EreconALL_uk->Write("EreconALL_uk");
	
	delete EreconALL_uk;
      }

      // SIMULATION

      { //Setting scope for second OctetAsymmetry so that it will be deleted to clear memory
	OctetAsymmetry SIM(octetNum, enBinWidth, 50., true, true, false);
	
	std::vector <double> superSum;
	std::vector <double> superSumError;
	std::vector <double> xErr(numBins,0.);

	f->cd();

	//Type 0 events
	SIM.calcSuperSum(4);
	superSum = SIM.returnSuperSum();
	superSumError = SIM.returnSuperSumError();

	for (int n=0; n<numBins;n++) {
	  if (enBins[n]>=enWinLow && enBins[n]<=enWinHigh) simIntegral+=superSum[n];
	  //std::cout << enBins[n] << " " << superSum[n] << " " << superSumError[n] << std::endl; 	  
	}
	normFactor = ukIntegral/simIntegral;
	for (int n=0; n<numBins;n++) {
	  superSum[n] = normFactor*superSum[n];
	  superSumError[n] = normFactor*superSumError[n];
	  superSumTotal_sim[0][n]+=superSumError[n]>0.?(1./power(superSumError[n],2))*superSum[n]:0.;
	  superSumTotalError_sim[0][n]+=superSumError[n]>0.?(1./power(superSumError[n],2)):0.;
	  //std::cout << enBins[n] << " " << superSum[n] << " " << superSumError[n] << std::endl; 	  
	}
	
	TGraphErrors *Erecon0_sim = new TGraphErrors(numBins,&enBins[0],&superSum[0],&xErr[0],&superSumError[0]);
	Erecon0_sim->SetMarkerColor(kBlue);
	Erecon0_sim->SetMarkerStyle(21);
	Erecon0_sim->GetXaxis()->SetTitle("Energy (keV)");
	Erecon0_sim->SetTitle("All Types Super-Sum");
	Erecon0_sim->Draw("AP");
	//Erecon0_uk->SetDirectory(f);
	
	Erecon0_sim->Write("Erecon0_sim");
	//f->Write();

	delete Erecon0_sim;
	
       
	//Type 1 events
	SIM.calcSuperSum(6);
	superSum = SIM.returnSuperSum();
	superSumError = SIM.returnSuperSumError();

	for (int n=0; n<numBins;n++) {
	  superSum[n] = normFactor*superSum[n];
	  superSumError[n] = normFactor*superSumError[n];
	  superSumTotal_sim[1][n]+=superSumError[n]>0.?(1./power(superSumError[n],2))*superSum[n]:0.;
	  superSumTotalError_sim[1][n]+=superSumError[n]>0.?(1./power(superSumError[n],2)):0.;
	  //std::cout << enBins[n] << " " << superSum[n] << " " << superSumError[n] << std::endl; 	  
	}
	
	TGraphErrors *Erecon1_sim = new TGraphErrors(numBins,&enBins[0],&superSum[0],&xErr[0],&superSumError[0]);
	Erecon1_sim->SetMarkerColor(kBlue);
	Erecon1_sim->SetMarkerStyle(21);
	Erecon1_sim->GetXaxis()->SetTitle("Energy (keV)");
	Erecon1_sim->SetTitle("All Types Super-Sum");
	Erecon1_sim->Draw("AP");
	
	Erecon1_sim->Write("Erecon1_sim");
	
	delete Erecon1_sim;

	//Type 2/3 events
	SIM.calcSuperSum(7);
	superSum = SIM.returnSuperSum();
	superSumError = SIM.returnSuperSumError();
	for (int n=0; n<numBins;n++) {
	  superSum[n] = normFactor*superSum[n];
	  superSumError[n] = normFactor*superSumError[n];
	  superSumTotal_sim[2][n]+=superSumError[n]>0.?(1./power(superSumError[n],2))*superSum[n]:0.;
	  superSumTotalError_sim[2][n]+=superSumError[n]>0.?(1./power(superSumError[n],2)):0.;
	  //std::cout << enBins[n] << " " << superSum[n] << " " << superSumError[n] << std::endl; 	  
	}
	
	TGraphErrors *Erecon23_sim = new TGraphErrors(numBins,&enBins[0],&superSum[0],&xErr[0],&superSumError[0]);
	Erecon23_sim->SetMarkerColor(kBlue);
	Erecon23_sim->SetMarkerStyle(21);
	Erecon23_sim->GetXaxis()->SetTitle("Energy (keV)");
	Erecon23_sim->SetTitle("All Types Super-Sum");
	Erecon23_sim->Draw("AP");
	
	Erecon23_sim->Write("Erecon23_sim");
	
	delete Erecon23_sim;

	//All event types
	SIM.calcSuperSum(1);
	superSum = SIM.returnSuperSum();
	superSumError = SIM.returnSuperSumError();

	for (int n=0; n<numBins;n++) {
	  superSum[n] = normFactor*superSum[n];
	  superSumError[n] = normFactor*superSumError[n];
	  superSumTotal_sim[3][n]+=superSumError[n]>0.?(1./power(superSumError[n],2))*superSum[n]:0.;
	  superSumTotalError_sim[3][n]+=superSumError[n]>0.?(1./power(superSumError[n],2)):0.;
	  //std::cout << enBins[n] << " " << superSum[n] << " " << superSumError[n] << std::endl; 	  
	}
	
	TGraphErrors *EreconALL_sim = new TGraphErrors(numBins,&enBins[0],&superSum[0],&xErr[0],&superSumError[0]);
	EreconALL_sim->SetMarkerColor(kBlue);
	EreconALL_sim->SetMarkerStyle(21);
	EreconALL_sim->GetXaxis()->SetTitle("Energy (keV)");
	EreconALL_sim->SetTitle("All Types Super-Sum");
	EreconALL_sim->Draw("AP");
	
	EreconALL_sim->Write("EreconALL_sim");
	
	delete EreconALL_sim;

	f->Close();

      }

    }

    catch(const char* ex){
      std::cerr << "Error: " << ex << std::endl;
    }
  }

  //Creating file which is summed over all octets in range
  TString fileName = "superSumPlots/SuperSum_octets_";
  fileName+= octetNumStart;
  fileName+= "-";
  fileName+= octetNumEnd;
  fileName+= ".root";
  TFile *f2 = new TFile(fileName,"RECREATE");
  f2->cd();
  for (int t=0; t<4; t++) {
    for (int bin=0; bin<numBins; bin++) {
      superSumTotal_uk[t][bin] = superSumTotalError_uk[t][bin]>0.? superSumTotal_uk[t][bin]/superSumTotalError_uk[t][bin] : 0.;
      superSumTotalError_uk[t][bin] = superSumTotalError_uk[t][bin]>0.? (1./TMath::Sqrt(superSumTotalError_uk[t][bin])) : 0.;
      //std::cout << enBins[bin] << " " << superSumTotal_uk[0][bin] << " " << superSumTotalError_uk[0][bin] << std::endl;
      
      superSumTotal_sim[t][bin] = superSumTotalError_sim[t][bin]>0.? superSumTotal_sim[t][bin]/superSumTotalError_sim[t][bin] : 0.;
      superSumTotalError_sim[t][bin] = superSumTotalError_sim[t][bin]>0.? (1./TMath::Sqrt(superSumTotalError_sim[t][bin])) : 0.;
    }
  }

  TGraphErrors *Erecon_uk, *Erecon_sim;
  
  for (int t=0;t<4;t++) {
    Erecon_uk = new TGraphErrors(numBins,&enBins[0],&superSumTotal_uk[t][0],&xErr[0],&superSumTotalError_uk[t][0]);
    Erecon_uk->SetMarkerColor(kBlue);
    Erecon_uk->SetMarkerStyle(21);
    Erecon_uk->GetXaxis()->SetTitle("Energy (keV)");
    //std::string title = std::string("Type ") + (t<2)?itos(t):(t==2?std::string("2/3"):std::string("ALL")) + std::string(" Super-Sum");
    TString title = TString("Type ") + ((t<2) ? itos(t) : (t==2 ? TString("23") : TString("ALL"))) + TString(" Super-Sum");
    Erecon_uk->SetTitle(title);
    //Erecon0_uk->Draw("AP");
    title = TString("Erecon") + ((t<2)? itos(t) :(t==2?TString("23"):TString("ALL"))) + TString("_uk");
    Erecon_uk->Write(title);
	
    delete Erecon_uk;

    Erecon_sim = new TGraphErrors(numBins,&enBins[0],&superSumTotal_sim[t][0],&xErr[0],&superSumTotalError_sim[t][0]);
    Erecon_sim->SetMarkerColor(kBlue);
    Erecon_sim->SetMarkerStyle(21);
    Erecon_sim->GetXaxis()->SetTitle("Energy (keV)");
    title = TString("Type ") + ((t<2) ? itos(t) : (t==2 ? TString("23") : TString("ALL"))) + TString(" Super-Sum");
    Erecon_sim->SetTitle(title);
    //Erecon0_sim->Draw("AP");
    title = TString("Erecon") + ((t<2)? itos(t) : (t==2?TString("23"):TString("ALL"))) + TString("_sim");
    Erecon_sim->Write(title);
	
    delete Erecon_sim;
  }

  f2->Close();
  
  
  return 0;

}

    
    
      /* std::cout << "uk Asym " << ukAsym << std::endl;
      std::cout << "sim Asym " << simAsym << std::endl;
      
      TCanvas cEvisALL("cEvisALL"," ", 1400., 600.);
      cEvisALL.Divide(2,1);
      cEvisALL.cd(1);
  
      EreconEALL_uk.Draw();
      EreconEALL_sim.Draw("SAME");
  
      TLegend legALL(0.55,0.75,0.85,0.875);
      legALL.AddEntry(&EreconEALL_sim,"  SIM","l");
      legALL.AddEntry(&EreconEALL_uk,"  UK","l");
      legALL.Draw();
  
      cEvisALL.cd(2);
  
      EreconWALL_uk.Draw();
      EreconWALL_sim.Draw("SAME");
  
      std::string pdfstart = pdfFileBase + "[";
      cEvisALL.Print(pdfstart.c_str());
      cEvisALL.Print(pdfFileBase.c_str());
  
  
      TCanvas cEvis0("cEvis0"," ", 1400., 600.);
      cEvis0.Divide(2,1);
      cEvis0.cd(1);
  
      EreconE0_uk.Draw();
      EreconE0_sim.Draw("SAME");
  
      TLegend leg0(0.55,0.75,0.85,0.875);
      leg0.AddEntry(&EreconE0_sim,"  SIM","l");
      leg0.AddEntry(&EreconE0_uk,"  UK","l");
      leg0.Draw();
  
      cEvis0.cd(2);
  
      EreconW0_uk.Draw();
      EreconW0_sim.Draw("SAME");
  
      cEvis0.Print(pdfFileBase.c_str());
  
  
      TCanvas cEvis1("cEvis1"," ", 1400., 600.);
      cEvis1.Divide(2,1);
      cEvis1.cd(1);
  
      EreconE1_uk.Draw();
      EreconE1_sim.Draw("SAME");
  
      TLegend leg1(0.55,0.75,0.85,0.875);
      leg1.AddEntry(&EreconE1_sim,"  SIM","l");
      leg1.AddEntry(&EreconE1_uk,"  UK","l");
      leg1.Draw();
  
      cEvis1.cd(2);
  
      EreconW1_uk.Draw();
      EreconW1_sim.Draw("SAME");
  
      cEvis1.Print(pdfFileBase.c_str());
    
  
      TCanvas cEvis23("cEvis23"," ", 1400., 600.);
      cEvis23.Divide(2,1);
      cEvis23.cd(1);
  
      EreconE23_uk.Draw();
      EreconE23_sim.Draw("SAME");
  
      TLegend leg23(0.55,0.75,0.85,0.875);
      leg23.AddEntry(&EreconE23_sim,"  SIM","l");
      leg23.AddEntry(&EreconE23_uk,"  UK","l");
      leg23.Draw();
  
      cEvis23.cd(2);
  
      EreconW23_uk.Draw();
      EreconW23_sim.Draw("SAME");
  
      cEvis23.Print(pdfFileBase.c_str());
  
  
      TCanvas cTypes("cTypes"," ", 1400., 600.);
      cTypes.Divide(2,1);
      cTypes.cd(1);
      gPad->SetLogy();
  
      Type_uk.Draw();
      Type_sim.Draw("SAME");
  
      Type_uk.SetMinimum(1);
  
      TLegend legTypes(0.55,0.75,0.85,0.875);
      legTypes.AddEntry(&Type_sim,"  SIM","l");
      legTypes.AddEntry(&Type_uk,"  UK","l");
      legTypes.Draw();
  
      std::cout << "UK Asym " << ukAsym << std::endl;

      cTypes.cd(2);

      TPaveText ukTxt(0.1, 0.25, 0.45, 0.75);
      Char_t temp[500];
      ukTxt.AddText("Data");
      //ukTxt.AddText(geometry.c_str());
      sprintf(temp, "Erecon = %.0f-%.0f keV", enWinLow, enWinHigh);
      ukTxt.AddText(temp);
      sprintf(temp, "Events: East %i   West %i", ukTotE, ukTotW);
      ukTxt.AddText(temp);
      sprintf(temp, "Total:  %i",ukTotE+ukTotW);
      ukTxt.AddText(temp);
      ukTxt.AddText(" ");
      sprintf(temp, "Raw Asymmetry:  %.5f",float(ukTotE-ukTotW)/float(ukTotE+ukTotW));
      ukTxt.AddText(temp);
      sprintf(temp, "Octet Asymmetry:  %.5f +/- %.5f",ukAsym,ukAsymError);
      ukTxt.AddText(temp);
      ukTxt.AddText(" ");
      ukTxt.AddText("Type:     East:        West:       "); 
      sprintf(temp, "0         %.5f   %.5f",  (float)uk0E/(float)ukTotE, (float)uk0W/(float)ukTotW);
      ukTxt.AddText(temp); 
      sprintf(temp, "I        %.5f   %.5f", (float)uk1E/(float)ukTotE, (float)uk1W/(float)ukTotW);
      ukTxt.AddText(temp); 
      sprintf(temp, "II/III       %.5f   %.5f",  (float)uk23E/(float)ukTotE, (float)uk23W/(float)ukTotW);
      ukTxt.AddText(temp); 
      ukTxt.Draw();

      TPaveText simTxt(0.55, 0.25, 0.9, 0.75);
      simTxt.AddText("Simulation");
      //simTxt.AddText(geometry.c_str());
      sprintf(temp, "Erecon = %.0f-%.0f keV", enWinLow, enWinHigh);
      simTxt.AddText(temp);
      sprintf(temp, "Events: East %i   West %i", simTotE, simTotW);
      simTxt.AddText(temp);
      sprintf(temp, "Total: %i",simTotE+simTotW);
      simTxt.AddText(temp);
      simTxt.AddText(" ");
      sprintf(temp, "Raw Asymmetry:  %.5f",float(simTotE-simTotW)/float(simTotE+simTotW));
      simTxt.AddText(temp);
      sprintf(temp, "Octet Asymmetry:  %.5f +/- %.5f",simAsym,simAsymError);
      simTxt.AddText(temp);
      simTxt.AddText(" ");
      simTxt.AddText("Type:     East:        West:       "); 
      sprintf(temp, "0         %.5f   %.5f",  (float)sim0E/(float)simTotE, (float)sim0W/(float)simTotW);
      simTxt.AddText(temp); 
      sprintf(temp, "I        %.5f   %.5f", (float)sim1E/(float)simTotE, (float)sim1W/(float)simTotW);
      simTxt.AddText(temp); 
      sprintf(temp, "II/III       %.5f   %.5f",  (float)sim23E/(float)simTotE, (float)sim23W/(float)simTotW);
      simTxt.AddText(temp); 
      simTxt.Draw();
  
      cTypes.Print(pdfFileBase.c_str());
      std::string pdfend = pdfFileBase + "]";
      cTypes.Print(pdfend.c_str());

    }

    catch(const char* ex){
      std::cerr << "Error: " << ex << std::endl;
    }
  }

  return 0;

}
      */    
	
 
