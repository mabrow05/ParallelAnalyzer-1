/*
  Script to compare UK and MPM data which has been run through the calibrations
  from each respective analyzer

  run using... root -l 'UK2MPM_comp.C(int runNumber, float enWinLow, float enWinHigh)'
 */

#include <string>
#include "../include/MButils.hh"
  
void UK2MPM_comp(int runNumberORoctet, Float_t EnWinLow, Float_t EnWinHigh) {

  gStyle->SetOptStat(0);

  TChain *mpmTree = new TChain("phys");
  TChain *ukTree = new TChain("pass4");

  std::string pdfFileBase;
  std::string mpmData;
  std::string ukData;

  if (runNumberORoctet<16000) {
    pdfFileBase = "calibrationRunComp_Octet"+itos(runNumberORoctet)+"_evis_"+itos((int)EnWinLow)+"-"+itos((int)EnWinHigh)+".pdf";
    std::string octFile = std::string(getenv("OCTET_LIST"))+"/All_Octets/octet_list_"+itos(runNumberORoctet)+".dat";
    ifstream runs(octFile.c_str());
    cout << "Opened " << octFile << "\n";
    std::string type;
    int run;
    while (runs >> type >> run) {
      if (type=="A2" || type=="A5" || type=="A7" || type=="A10" || type=="B2" || type=="B5" || type=="B7" || type=="B10") {
	mpmData = std::string(getenv("UCNAOUTPUTDIR"))+"/hists/spec_" + itos(run) + ".root";
	ukData = std::string(getenv("REPLAY_PASS4"))+"/replay_pass4_" + itos(run) + ".root";
	mpmTree->Add(mpmData.c_str());
	ukTree->Add(ukData.c_str());
	cout << "Added " << ukData << " to chain\n";
      }
    }
  }

  else {
    pdfFileBase = "calibrationRunComp_run"+itos(runNumberORoctet)+"_evis_"+itos((int)EnWinLow)+"-"+itos((int)EnWinHigh)+".pdf";
    mpmData = std::string(getenv("UCNAOUTPUTDIR"))+"/hists/spec_" + itos(runNumberORoctet) + ".root";
    ukData = std::string(getenv("REPLAY_PASS4"))+"/replay_pass4_" + itos(runNumberORoctet) + ".root";
    mpmTree->Add(mpmData.c_str());
    ukTree->Add(ukData.c_str());
  }


  //Making histograms for comparison  
  TH1F *EvisEALL_mpm = new TH1F("EvisEALL_mpm","EvisE Type ALL", 100, 0., 800.);
  TH1F *EvisWALL_mpm = new TH1F("EvisWALL_mpm","EvisW Type ALL", 100, 0., 800.);
  TH1F *EvisEALL_uk = new TH1F("EvisEALL_uk","EvisE Type ALL", 100, 0., 800.);
  EvisEALL_uk->GetXaxis()->SetTitle("Energy (keV)");
  TH1F *EvisWALL_uk = new TH1F("EvisWALL_uk","EvisW Type ALL", 100, 0., 800.);
  EvisWALL_uk->GetXaxis()->SetTitle("Energy (keV)");

  TH1F *EvisE0_mpm = new TH1F("EvisE0_mpm","EvisE Type 0", 100, 0., 800.);
  TH1F *EvisW0_mpm = new TH1F("EvisW0_mpm","EvisW Type 0", 100, 0., 800.);
  TH1F *EvisE0_uk = new TH1F("EvisE0_uk","EvisE Type 0", 100, 0., 800.);
  EvisE0_uk->GetXaxis()->SetTitle("Energy (keV)");
  TH1F *EvisW0_uk = new TH1F("EvisW0_uk","EvisW Type 0", 100, 0., 800.);
  EvisW0_uk->GetXaxis()->SetTitle("Energy (keV)");

  TH1F *EvisE1_mpm = new TH1F("EvisE1_mpm","EvisE Type 1", 200, 0., 800.);
  TH1F *EvisW1_mpm = new TH1F("EvisW1_mpm","EvisW Type 1", 200, 0., 800.);
  TH1F *EvisE1_uk = new TH1F("EvisE1_uk","EvisE Type 1", 200, 0., 800.);
  EvisE1_uk->GetXaxis()->SetTitle("Energy (keV)");
  TH1F *EvisW1_uk = new TH1F("EvisW1_uk","EvisW Type 1", 200, 0., 800.);
  EvisW1_uk->GetXaxis()->SetTitle("Energy (keV)");

  TH1F *EvisE23_mpm = new TH1F("EvisE23_mpm","EvisE Type 2/3", 200, 0., 800.);
  TH1F *EvisW23_mpm = new TH1F("EvisW23_mpm","EvisW Type 2/3", 200, 0., 800.);
  TH1F *EvisE23_uk = new TH1F("EvisE23_uk","EvisE Type 2/3", 200, 0., 800.);
  EvisE23_uk->GetXaxis()->SetTitle("Energy (keV)");
  TH1F *EvisW23_uk = new TH1F("EvisW23_uk","EvisW Type 2/3", 200, 0., 800.); 
  EvisW23_uk->GetXaxis()->SetTitle("Energy (keV)");

  TH1I *Type_mpm = new TH1I("Type_mpm","Event Types", 4, 0, 4);
  TH1I *Type_uk = new TH1I("Type_uk","Event Types", 4, 0, 4);
  

  EvisEALL_mpm->SetLineColor(kRed);
  EvisWALL_mpm->SetLineColor(kRed);
  EvisEALL_uk->SetLineColor(kBlue);
  EvisWALL_uk->SetLineColor(kBlue);
  
  EvisE0_mpm->SetLineColor(kRed);
  EvisW0_mpm->SetLineColor(kRed);
  EvisE0_uk->SetLineColor(kBlue);
  EvisW0_uk->SetLineColor(kBlue);

  EvisE1_mpm->SetLineColor(kRed);
  EvisW1_mpm->SetLineColor(kRed);
  EvisE1_uk->SetLineColor(kBlue);
  EvisW1_uk->SetLineColor(kBlue);

  EvisE23_mpm->SetLineColor(kRed);
  EvisW23_mpm->SetLineColor(kRed);
  EvisE23_uk->SetLineColor(kBlue);
  EvisW23_uk->SetLineColor(kBlue);

  Type_mpm->SetLineColor(kRed);
  Type_uk->SetLineColor(kBlue);
  

  //Plot what is of interest

  TCanvas *cEvisALL = new TCanvas("cEvisALL"," ", 1400., 600.);
  cEvisALL->Divide(2,1);
  cEvisALL->cd(1);
  
  ukTree->Draw("EvisE>>EvisEALL_uk", "Type<4 && Type>=0 && Side==0 && PID==1 && (xE.center*xE.center+yE.center*yE.center)<50.*50.");
  mpmTree->Draw("EvisE>>EvisEALL_mpm", "Type<4 && Type>=0 && Side==0 && PID==1 && (xEmpm.center*xEmpm.center+yEmpm.center*yEmpm.center)<50.*50.","SAME");

  TLegend *legALL = new TLegend(0.55,0.75,0.85,0.875);
  legALL->AddEntry(EvisEALL_mpm,"  MPM","l");
  legALL->AddEntry(EvisEALL_uk,"  UK","l");
  legALL->Draw();

  cEvisALL->cd(2);
  
  ukTree->Draw("EvisW>>EvisWALL_uk", "Type<4 && Type>=0 && Side==1 && PID==1 && (xW.center*xW.center+yW.center*yW.center)<50.*50.");
  mpmTree->Draw("EvisW>>EvisWALL_mpm", "Type<4 && Type>=0 && Side==1 && PID==1 && (xWmpm.center*xWmpm.center+yWmpm.center*yWmpm.center)<50.*50.", "SAME");

  std::string pdfstart = pdfFileBase + "[";
  cEvisALL->Print(pdfstart.c_str());
  cEvisALL->Print(pdfFileBase.c_str());

  TCanvas *cEvis0 = new TCanvas("cEvis0"," ", 1400., 600.);
  cEvis0->Divide(2,1);
  cEvis0->cd(1);
  
  ukTree->Draw("EvisE>>EvisE0_uk", "Type==0 && Side==0 && PID==1 && (xE.center*xE.center+yE.center*yE.center)<50.*50.");
  mpmTree->Draw("EvisE>>EvisE0_mpm", "Type==0 && Side==0 && PID==1 && (xEmpm.center*xEmpm.center+yEmpm.center*yEmpm.center)<50.*50.", "SAME");

  TLegend *leg0 = new TLegend(0.55,0.75,0.85,0.875);
  leg0->AddEntry(EvisE0_mpm,"  MPM","l");
  leg0->AddEntry(EvisE0_uk,"  UK","l");
  leg0->Draw();

  cEvis0->cd(2);
  
  ukTree->Draw("EvisW>>EvisW0_uk", "Type==0 && Side==1 && PID==1 && (xW.center*xW.center+yW.center*yW.center)<50.*50.");
  mpmTree->Draw("EvisW>>EvisW0_mpm", "Type==0 && Side==1 && PID==1 && (xWmpm.center*xWmpm.center+yWmpm.center*yWmpm.center)<50.*50.", "SAME");

  cEvis0->Print(pdfFileBase.c_str());

  TCanvas *cEvis1 = new TCanvas("cEvis1"," ", 1400., 600.);
  cEvis1->Divide(2,1);
  cEvis1->cd(1);
  
  ukTree->Draw("EvisE>>EvisE1_uk", "Type==1 && Side==0 && PID==1 && (xE.center*xE.center+yE.center*yE.center)<50.*50.");
  mpmTree->Draw("EvisE>>EvisE1_mpm", "Type==1 && Side==0 && PID==1 && (xEmpm.center*xEmpm.center+yEmpm.center*yEmpm.center)<50.*50.", "SAME");

  TLegend *leg1 = new TLegend(0.55,0.75,0.85,0.875);
  leg1->AddEntry(EvisE1_mpm,"  MPM","l");
  leg1->AddEntry(EvisE1_uk,"  UK","l");
  leg1->Draw();

  cEvis1->cd(2);
  
  ukTree->Draw("EvisW>>EvisW1_uk", "Type==1 && Side==1 && PID==1 && (xW.center*xW.center+yW.center*yW.center)<50.*50.");
  mpmTree->Draw("EvisW>>EvisW1_mpm", "Type==1 && Side==1 && PID==1 && (xWmpm.center*xWmpm.center+yWmpm.center*yWmpm.center)<50.*50.", "SAME");

  cEvis1->Print(pdfFileBase.c_str());

  TCanvas *cEvis23 = new TCanvas("cEvis23"," ", 1400., 600.);
  cEvis23->Divide(2,1);
  cEvis23->cd(1);
  
  ukTree->Draw("EvisE>>EvisE23_uk", "(Type==2 || Type==3) && Side==0 && PID==1 && (xE.center*xE.center+yE.center*yE.center)<50.*50.");
  mpmTree->Draw("EvisE>>EvisE23_mpm", "(Type==2 || Type==3) && Side==0 && PID==1 && (xEmpm.center*xEmpm.center+yEmpm.center*yEmpm.center)<50.*50.", "SAME");

  TLegend *leg23 = new TLegend(0.55,0.75,0.85,0.875);
  leg23->AddEntry(EvisE23_mpm,"  MPM","l");
  leg23->AddEntry(EvisE23_uk,"  UK","l");
  leg23->Draw();

  cEvis23->cd(2);
  
  ukTree->Draw("EvisW>>EvisW23_uk", "(Type==2 || Type==3) && Side==1 && PID==1 && (xW.center*xW.center+yW.center*yW.center)<50.*50.");
  mpmTree->Draw("EvisW>>EvisW23_mpm", "(Type==2 || Type==3) && Side==1 && PID==1 && (xWmpm.center*xWmpm.center+yWmpm.center*yWmpm.center)<50.*50.", "SAME");

  cEvis23->Print(pdfFileBase.c_str());

  TCanvas *cTypes = new TCanvas("cTypes"," ", 1400., 600.);
  cTypes->Divide(2,1);
  cTypes->cd(1);
  gPad->SetLogy();

  std::string TypeCut = "PID==1 && (EvisE+EvisW)>"+ftos(EnWinLow)+" && (EvisE+EvisW)<"+ftos(EnWinHigh);
  ukTree->Draw("Type>>Type_uk", TypeCut.c_str());
  mpmTree->Draw("Type>>Type_mpm", TypeCut.c_str(), "SAME");

  Type_uk->SetMinimum(1);

  TLegend *legTypes = new TLegend(0.55,0.75,0.85,0.875);
  legTypes->AddEntry(Type_mpm,"  MPM","l");
  legTypes->AddEntry(Type_uk,"  UK","l");
  legTypes->Draw();
  
  //Get the number of each event type to compare fractions
  Int_t uk0E, mpm0E, uk1E, mpm1E, uk23E, mpm23E, ukTotE, mpmTotE, uk0W, mpm0W, uk1W, mpm1W, uk23W, mpm23W, ukTotW, mpmTotW;
  std::string eastCutBase = "PID==1 && (EvisE+EvisW)>"+ftos(EnWinLow)+" && (EvisE+EvisW)<"+ftos(EnWinHigh)+" && Side==0 && (xE.center*xE.center+yE.center*yE.center)<50.*50.";
  std::string westCutBase = "PID==1 && (EvisE+EvisW)>"+ftos(EnWinLow)+" && (EvisE+EvisW)<"+ftos(EnWinHigh)+" && Side==1 && (xW.center*xW.center+yW.center*yW.center)<50.*50.";
  //std::string eastCutBase = "PID==1 && (Erecon)>"+ftos(EnWinLow)+" && (Erecon)<"+ftos(EnWinHigh)+" && Side==0 && (xE.center*xE.center+yE.center*yE.center)<50.*50.";
  //std::string westCutBase = "PID==1 && (Erecon)>"+ftos(EnWinLow)+" && (Erecon)<"+ftos(EnWinHigh)+" && Side==1 && (xW.center*xW.center+yW.center*yW.center)<50.*50.";
  //std::string eastCutBase = "PID==1 && Side==0";
  //std::string westCutBase = "PID==1 && Side==1";
    
  std::string cut = eastCutBase + " && Type==0";
  uk0E = ukTree->GetEntries(cut.c_str());
  cut = eastCutBase + " && Type==1";
  uk1E = ukTree->GetEntries(cut.c_str());
  cut = eastCutBase + " && (Type==2 || Type==3)";
  uk23E = ukTree->GetEntries(cut.c_str());
  
  ukTotE = uk0E + uk1E + uk23E;

  cut = westCutBase + " && Type==0";
  uk0W = ukTree->GetEntries(cut.c_str());
  cut = westCutBase + " && Type==1";
  uk1W = ukTree->GetEntries(cut.c_str());
  cut = westCutBase + " && (Type==2 || Type==3)";
  uk23W = ukTree->GetEntries(cut.c_str());
  
  ukTotW = uk0W + uk1W + uk23W;

  eastCutBase = "PID==1 && (EvisE+EvisW)>"+ftos(EnWinLow)+" && (EvisE+EvisW)<"+ftos(EnWinHigh)+" && Side==0 && (xEmpm.center*xEmpm.center+yEmpm.center*yEmpm.center)<50.*50.";
  westCutBase = "PID==1 && (EvisE+EvisW)>"+ftos(EnWinLow)+" && (EvisE+EvisW)<"+ftos(EnWinHigh)+" && Side==1 && (xWmpm.center*xWmpm.center+yWmpm.center*yWmpm.center)<50.*50.";

  cut = eastCutBase + " && Type==0";
  mpm0E = mpmTree->GetEntries(cut.c_str());
  cut = eastCutBase + " && Type==1";
  mpm1E = mpmTree->GetEntries(cut.c_str());
  cut = eastCutBase + " && (Type==2 || Type==3)";
  mpm23E = mpmTree->GetEntries(cut.c_str());
  
  mpmTotE = mpm0E + mpm1E + mpm23E;

  cut = westCutBase + " && Type==0";
  mpm0W = mpmTree->GetEntries(cut.c_str());
  cut = westCutBase + " && Type==1";
  mpm1W = mpmTree->GetEntries(cut.c_str());
  cut = westCutBase + " && (Type==2 || Type==3)";
  mpm23W = mpmTree->GetEntries(cut.c_str());
  
  mpmTotW = mpm0W + mpm1W + mpm23W;
  
  
  cTypes->cd(2);

  TPaveText *ukTxt = new TPaveText(0.1, 0.25, 0.45, 0.75);
  Char_t temp[500];
  ukTxt->AddText("UK Analyzer");
  //ukTxt->AddText(geometry.c_str());
  sprintf(temp, "Evis = %.0f-%.0f keV", EnWinLow, EnWinHigh);
  ukTxt->AddText(temp);
  sprintf(temp, "Events: East %i   West %i", ukTotE, ukTotW);
  ukTxt->AddText(temp);
  sprintf(temp, "Total:  %i",ukTotE+ukTotW);
  ukTxt->AddText(temp);
  ukTxt->AddText(" ");
  sprintf(temp, "Raw Asymmetry:  %.5f",float(ukTotE-ukTotW)/float(ukTotE+ukTotW));
  ukTxt->AddText(temp);
  ukTxt->AddText(" ");
  ukTxt->AddText("Type:     East:        West:       "); 
  sprintf(temp, "0         %.5f   %.5f",  (float)uk0E/(float)ukTotE, (float)uk0W/(float)ukTotW);
  ukTxt->AddText(temp); 
  sprintf(temp, "I        %.5f   %.5f", (float)uk1E/(float)ukTotE, (float)uk1W/(float)ukTotW);
  ukTxt->AddText(temp); 
  sprintf(temp, "II/III       %.5f   %.5f",  (float)uk23E/(float)ukTotE, (float)uk23W/(float)ukTotW);
  ukTxt->AddText(temp); 
  ukTxt->Draw();

  TPaveText *mpmTxt = new TPaveText(0.55, 0.25, 0.9, 0.75);
  Char_t temp[500];
  mpmTxt->AddText("MPM Analyzer");
  //mpmTxt->AddText(geometry.c_str());
  sprintf(temp, "Evis = %.0f-%.0f keV", EnWinLow, EnWinHigh);
  mpmTxt->AddText(temp);
  sprintf(temp, "Events: East %i   West %i", mpmTotE, mpmTotW);
  mpmTxt->AddText(temp);
  sprintf(temp, "Total: %i",mpmTotE+mpmTotW);
  mpmTxt->AddText(temp);
  mpmTxt->AddText(" ");
  sprintf(temp, "Raw Asymmetry:  %.5f",float(mpmTotE-mpmTotW)/float(mpmTotE+mpmTotW));
  mpmTxt->AddText(temp);
  mpmTxt->AddText(" ");
  mpmTxt->AddText("Type:     East:        West:       "); 
  sprintf(temp, "0         %.5f   %.5f",  (float)mpm0E/(float)mpmTotE, (float)mpm0W/(float)mpmTotW);
  mpmTxt->AddText(temp); 
  sprintf(temp, "I        %.5f   %.5f", (float)mpm1E/(float)mpmTotE, (float)mpm1W/(float)mpmTotW);
  mpmTxt->AddText(temp); 
  sprintf(temp, "II/III       %.5f   %.5f",  (float)mpm23E/(float)mpmTotE, (float)mpm23W/(float)mpmTotW);
  mpmTxt->AddText(temp); 
  mpmTxt->Draw();

  cTypes->Print(pdfFileBase.c_str());
  std::string pdfend = pdfFileBase + "]";
  cTypes->Print(pdfend.c_str());
}
