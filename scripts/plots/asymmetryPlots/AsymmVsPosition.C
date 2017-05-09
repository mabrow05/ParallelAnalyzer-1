

void AsymmVsPosition(TString year, TString anaCh, TString corrections, bool withPOL, Int_t ebinLow=220, Int_t ebinHigh=680) {

  
  const Int_t numQuads = 4;
  const Int_t numRings = 5;
  TString Quads[numQuads] = {"I","II","III","IV"};//
  TString Rings[numRings] = {"0-10","10-20","20-30","30-40","40-50"};
 
  std::vector <Double_t> AsymmDataQuad;
  std::vector <Double_t> AsymmDataErrQuad; //Statistical Errors
  std::vector <Double_t> AsymmSimQuad;
  std::vector <Double_t> AsymmSimErrQuad; //Statistical Errors

  std::vector <Double_t> AsymmDataRings;
  std::vector <Double_t> AsymmDataErrRings; //Statistical Errors
  std::vector <Double_t> AsymmSimRings;
  std::vector <Double_t> AsymmSimErrRings; //Statistical Errors

  
  
  //Read in Asymmetries
  ifstream infile;
  std::string strHold;
    
  for (Int_t i=0; i<numQuads; ++i) {
    infile.open(TString::Format("%s/Asymmetries/%s%s_OctetAsymmetries_AnaCh%s_%i-%i_Octets_%s_Quadrant%i.txt",
				getenv("ANALYSIS_RESULTS"),corrections.Data(),
				(withPOL?"_withPOL":""),anaCh.Data(),
				ebinLow,ebinHigh,year==TString("2011-2012")?"0-59":"60-121",i).Data());
    
    Double_t AsymHold = 0.;
    Double_t AsymErrHold = 0.;
    if (infile.is_open()) {
      for (Int_t j=0; j<3; j++) infile >> strHold >> AsymHold >> AsymErrHold;
      cout << strHold << AsymHold << AsymErrHold << endl;
    }
    AsymmDataQuad.push_back(AsymHold);
    AsymmDataErrQuad.push_back(AsymErrHold);
    
    infile.close();
    
    infile.open(TString::Format("%s/Asymmetries/%s_OctetAsymmetries_AnaCh%s_%i-%i_Octets_%s_Quadrant%i.txt",
				getenv("SIM_ANALYSIS_RESULTS"),corrections.Data(),
				anaCh.Data(),ebinLow,ebinHigh,
				year==TString("2011-2012")?"0-59":"60-121",i).Data());
    
    AsymHold = 0.;
    AsymErrHold = 0.;
    if (infile.is_open()) {
      for (Int_t j=0; j<3; j++) infile >> strHold >> AsymHold >> AsymErrHold;
    }
    AsymmSimQuad.push_back(AsymHold);
    AsymmSimErrQuad.push_back(AsymErrHold);
    
    infile.close();
  }
  
  for (Int_t i=0; i<numRings; i++) {
    infile.open(TString::Format("%s/Asymmetries/%s%s_OctetAsymmetries_AnaCh%s_%i-%i_Octets_%s_RadialRing%i.txt",
				getenv("ANALYSIS_RESULTS"),corrections.Data(),
				(withPOL?"_withPOL":""),anaCh.Data(),
				ebinLow,ebinHigh,year==TString("2011-2012")?"0-59":"60-121",i).Data());
    
    Double_t AsymHold = 0.;
    Double_t AsymErrHold = 0.;
    if (infile.is_open()) {
      for (Int_t j=0; j<3; j++) infile >> strHold >> AsymHold >> AsymErrHold;
    }
    AsymmDataRings.push_back(AsymHold);
    AsymmDataErrRings.push_back(AsymErrHold);
    
    infile.close();
    
    infile.open(TString::Format("%s/Asymmetries/%s_OctetAsymmetries_AnaCh%s_%i-%i_Octets_%s_RadialRing%i.txt",
				getenv("SIM_ANALYSIS_RESULTS"),corrections.Data(),
				anaCh.Data(),ebinLow,ebinHigh,
				year==TString("2011-2012")?"0-59":"60-121",i).Data());
    
    AsymHold = 0.;
    AsymErrHold = 0.;
    if (infile.is_open()) {
      for (Int_t j=0; j<3; j++) infile >> strHold >> AsymHold >> AsymErrHold;
    }
    AsymmSimRings.push_back(AsymHold);
    AsymmSimErrRings.push_back(AsymErrHold);
    
    infile.close();
  }
  
  

  

  Double_t xerr[10] = {0.};

  //cout << anaChoices.size() << " " << anaChoices[0] << " " << AsymmData[0] << " " << xerr[0] << " " << AsymmDataErr[0] << endl;
    

  gStyle->SetLegendBorderSize(0);
  gStyle->SetFillStyle(0);
  gStyle->SetPadLeftMargin(0.1);
  gStyle->SetTitleXOffset(1.2);

  TCanvas *c1 = new TCanvas("c1","demo bin labels",10,10,600,600);

  c1->SetLeftMargin(0.15);
  c1->SetBottomMargin(0.15);

  
  TH1F * data = new TH1F("data","Asymmetry vs. Quadrant",AsymmDataQuad.size(),0.5,AsymmDataQuad.size()+0.5);
  data->SetMarkerColor(kBlue);
  data->SetLineColor(kBlue);
  data->SetLineWidth(3);
  data->SetMarkerStyle(kFullCircle);
  //data->GetXaxis()->SetNdivisions();
  //data->SetCanExtend(TH1::kAllAxes);
  data->SetStats(0);

  for (Int_t i=0;i<numQuads;++i) {
    data->SetBinContent(i+1,AsymmDataQuad[i]);
    data->SetBinError(i+1,AsymmDataErrQuad[i]);
  }
  //data->LabelsDeflate("X");
  //data->LabelsDeflate("Y");
  //data->LabelsOption("v");
  for (int i=1; i<=numQuads; ++i) data->GetXaxis()->SetBinLabel(i,Quads[i-1]);
  
  data->GetXaxis()->SetLabelSize(0.07);
  data->GetXaxis()->SetTitle("Quadrant");
  data->GetYaxis()->SetTitle("Asymmetry");
  data->GetYaxis()->SetTitleOffset(2.);
  data->GetXaxis()->SetTitleSize(0.05);
  data->GetXaxis()->CenterTitle();
  data->GetYaxis()->CenterTitle();
  data->SetMinimum(-0.135);
  data->SetMaximum(-0.115);

  data->Draw("EX0");

  TH1F * sim = new TH1F("sim","Asymmetry vs. Quadrant",AsymmSimQuad.size(),0.5,AsymmSimQuad.size()+0.5);
  sim->SetMarkerColor(kRed);
  sim->SetLineColor(kRed);
  sim->SetLineWidth(3);
  sim->SetMarkerStyle(kFullSquare);
  //sim->GetXaxis()->SetNdivisions();
  //sim->SetCanExtend(TH1::kAllAxes);
  sim->SetStats(0);

  for (Int_t i=0;i<numQuads;++i) {
    sim->SetBinContent(i+1,AsymmSimQuad[i]);
    sim->SetBinError(i+1,AsymmSimErrQuad[i]);
  }
  
  sim->Draw("SAME EX0");
  
  
  TLegend *leg = new TLegend(0.55,0.75,0.9,0.9);
  //leg->SetHeader("The Legend Title","C"); // option "C" allows to center the header
  leg->AddEntry(data,"Data");
  leg->AddEntry(sim,"Monte Carlo");
   //leg->AddEntry("gr","Graph with error bars","lep");
  leg->Draw();


  ///////////////// Ring analysis /////////////////////////

  TCanvas *c2 = new TCanvas("c2","demo bin labels",10,10,600,600);

  c2->SetLeftMargin(0.15);
  c2->SetBottomMargin(0.15);

  
  TH1F * dataRing = new TH1F("dataRing","Asymmetry vs. Radial Position",AsymmDataRings.size(),0.5,AsymmDataRings.size()+0.5);
  dataRing->SetMarkerColor(kBlue);
  dataRing->SetLineColor(kBlue);
  dataRing->SetLineWidth(3);
  dataRing->SetMarkerStyle(kFullCircle);
  //dataRing->GetXaxis()->SetNdivisions();
  //dataRing->SetCanExtend(TH1::kAllAxes);
  dataRing->SetStats(0);

  for (Int_t i=0;i<numRings;++i) {
    dataRing->SetBinContent(i+1,AsymmDataRings[i]);
    dataRing->SetBinError(i+1,AsymmDataErrRings[i]);
  }
  //dataRing->LabelsDeflate("X");
  //dataRing->LabelsDeflate("Y");
  //dataRing->LabelsOption("v");
  for (int i=1; i<=numRings; ++i) dataRing->GetXaxis()->SetBinLabel(i,Rings[i-1]);
  
  dataRing->GetXaxis()->SetLabelSize(0.07);
  dataRing->GetXaxis()->SetTitle("Radial Cut (mm)");
  dataRing->GetYaxis()->SetTitle("Asymmetry");
  dataRing->GetYaxis()->SetTitleOffset(2.);
  dataRing->GetXaxis()->SetTitleSize(0.05);
  dataRing->GetXaxis()->CenterTitle();
  dataRing->GetYaxis()->CenterTitle();
  dataRing->SetMinimum(-0.135);
  dataRing->SetMaximum(-0.115);

  dataRing->Draw("EX0");

  TH1F * simRing = new TH1F("simRing","Asymmetry vs. Radial Position",AsymmSimRings.size(),0.5,AsymmSimRings.size()+0.5);
  simRing->SetMarkerColor(kRed);
  simRing->SetLineColor(kRed);
  simRing->SetLineWidth(3);
  simRing->SetMarkerStyle(kFullSquare);
  //simRing->GetXaxis()->SetNdivisions();
  //simRing->SetCanExtend(TH1::kAllAxes);
  simRing->SetStats(0);

  for (Int_t i=0;i<numRings;++i) {
    simRing->SetBinContent(i+1,AsymmSimRings[i]);
    simRing->SetBinError(i+1,AsymmSimErrRings[i]);
  }
  
  simRing->Draw("SAME EX0");
  
  
  TLegend *leg = new TLegend(0.55,0.75,0.9,0.9);
  //leg->SetHeader("The Legend Title","C"); // option "C" allows to center the header
  leg->AddEntry(dataRing,"Data");
  leg->AddEntry(simRing,"Monte Carlo");
   //leg->AddEntry("gr","Graph with error bars","lep");
  leg->Draw();
  

  c1->Print(TString::Format("AsymmVsPosition_%s_anaCh%s_%s_%s_%i-%i.pdf(",year.Data(),anaCh.Data(),corrections.Data(),withPOL?"withPOL":"noPOL",ebinLow,ebinHigh));
  c2->Print(TString::Format("AsymmVsPosition_%s_anaCh%s_%s_%s_%i-%i.pdf)",year.Data(),anaCh.Data(),corrections.Data(),withPOL?"withPOL":"noPOL",ebinLow,ebinHigh));

}
