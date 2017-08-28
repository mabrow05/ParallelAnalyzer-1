

void AsymmByAnach(TString corrections, bool withPOL, Int_t ebinLow=220, Int_t ebinHigh=670) {
  
  bool color = false;
  
  bool withSim = false;//true;
  bool CorrAndUnCorr = true;//false;

  double ymin0 = -0.128;
  double ymax0 = -0.108;
  
  double yminBS = -0.2;
  double ymaxBS = 0.01;

  const Int_t numAnaChT0 = 9;
  const Int_t numAnaChBacksc = 9;
  TString anChALL[numAnaChT0] = {"A","B","C","D","F","G","H","J","K"};//
  TString anChT0[numAnaChT0] = {"A","B","C","D","","","","",""};//
  TString anChBacksc[numAnaChBacksc] = {"","","","","F","G","H","J","K"};

  Int_t markerUncorr0 = 24;
  Int_t markerCorr0 = 20;
  Int_t markerUncorrBS = 25;
  Int_t markerCorrBS = 21;
  Int_t colorBS = color?2:1;
  Int_t color0 = color?4:1;

  TString year = "2011-2012";
  

  std::vector <TString> anaChoicesALL;
  std::vector <Double_t> indicesALL;
  std::vector <Double_t> AsymmDataALL;
  std::vector <Double_t> AsymmDataErrALL; //Statistical Errors
  std::vector <Double_t> AsymmCorrDataALL;
  std::vector <Double_t> AsymmCorrDataErrALL; //Statistical Errors
 
  std::vector <TString> anaChoicesT0;
  std::vector <Double_t> indicesT0;
  std::vector <Double_t> AsymmDataT0;
  std::vector <Double_t> AsymmDataErrT0; //Statistical Errors
  std::vector <Double_t> AsymmCorrDataT0;
  std::vector <Double_t> AsymmCorrDataErrT0; //Statistical Errors
  

  std::vector <TString> anaChoicesBacksc;
  std::vector <Double_t> indicesBacksc;
  std::vector <Double_t> AsymmDataBacksc;
  std::vector <Double_t> AsymmDataErrBacksc; //Statistical Errors
  std::vector <Double_t> AsymmCorrDataBacksc;
  std::vector <Double_t> AsymmCorrDataErrBacksc; //Statistical Errors
  std::vector <Double_t> AsymmSimBacksc;
  std::vector <Double_t> AsymmSimErrBacksc; //Statistical Errors

  //Fill anaChoices
  for (Int_t i=0; i< sizeof(anChT0)/sizeof(TString); i++) 
    { anaChoicesT0.push_back(anChT0[i]); indicesT0.push_back(i+1);
      anaChoicesALL.push_back(anChALL[i]); indicesALL.push_back(i+1);}

  for (Int_t i=0; i< sizeof(anChBacksc)/sizeof(TString); i++) 
    { anaChoicesBacksc.push_back(anChBacksc[i]); indicesBacksc.push_back(i+1); }
  
  
  //Read in Asymmetries
  ifstream infile;
  std::string strHold;
  
  for (Int_t i=0; i<anaChoicesT0.size(); i++) {
    infile.open(TString::Format("%s/Asymmetries/UNBLINDED_%s%s_OctetAsymmetries_AnaCh%s_%i-%i_Octets_%s.txt",
				getenv("ANALYSIS_RESULTS"),corrections.Data(),
				(withPOL?"_withPOL":""),anaChoicesT0[i].Data(),
				ebinLow,ebinHigh,year==TString("2011-2012")?"0-59":"60-121").Data());
    
    Double_t AsymHold = 0.;
    Double_t AsymErrHold = 0.;
    if (infile.is_open()) {
      for (Int_t j=0; j<3; j++) infile >> strHold >> AsymHold >> AsymErrHold;
    }
    AsymmDataT0.push_back(AsymHold);
    AsymmDataErrT0.push_back(AsymErrHold);
    
    infile.close();
    
    infile.open(TString::Format("%s/Asymmetries/UNBLINDED_%s%s_OctetAsymmetries_AnaCh%s_%i-%i_Octets_%s.txt",
				getenv("ANALYSIS_RESULTS"),"AllCorr",
				(withPOL?"_withPOL":""),anaChoicesT0[i].Data(),
				ebinLow,ebinHigh,year==TString("2011-2012")?"0-59":"60-121").Data());
    
    AsymHold = 0.;
    AsymErrHold = 0.;
    if (infile.is_open()) {
      for (Int_t j=0; j<3; j++) infile >> strHold >> AsymHold >> AsymErrHold;
    }
    AsymmCorrDataT0.push_back(AsymHold);
    AsymmCorrDataErrT0.push_back(AsymErrHold);
    
    infile.close();
    
    infile.open(TString::Format("%s/Asymmetries/UNBLINDED_%s%s_OctetAsymmetries_AnaCh%s_%i-%i_Octets_%s.txt",
				getenv("ANALYSIS_RESULTS"),corrections.Data(),(withPOL?"_withPOL":""),
				anaChoicesALL[i].Data(),ebinLow,ebinHigh,
				year==TString("2011-2012")?"0-59":"60-121").Data() );
    
    AsymHold = 0.;
    AsymErrHold = 0.;
    if (infile.is_open()) {
      for (Int_t j=0; j<3; j++) infile >> strHold >> AsymHold >> AsymErrHold;
    }
    AsymmDataALL.push_back(AsymHold);
    AsymmDataErrALL.push_back(AsymErrHold);
    
    infile.close();

    infile.open(TString::Format("%s/Asymmetries/UNBLINDED_%s%s_OctetAsymmetries_AnaCh%s_%i-%i_Octets_%s.txt",
				getenv("ANALYSIS_RESULTS"),"AllCorr",(withPOL?"_withPOL":""),
				anaChoicesALL[i].Data(),ebinLow,ebinHigh,
				year==TString("2011-2012")?"0-59":"60-121").Data() );
    
    AsymHold = 0.;
    AsymErrHold = 0.;
    if (infile.is_open()) {
      for (Int_t j=0; j<3; j++) infile >> strHold >> AsymHold >> AsymErrHold;
    }
    AsymmCorrDataALL.push_back(AsymHold);
    AsymmCorrDataErrALL.push_back(AsymErrHold);
    
    infile.close();
  }
  
  for (Int_t i=0; i<anaChoicesBacksc.size(); i++) {
    infile.open(TString::Format("%s/Asymmetries/UNBLINDED_%s%s_OctetAsymmetries_AnaCh%s_%i-%i_Octets_%s.txt",
				getenv("ANALYSIS_RESULTS"),corrections.Data(),
				(withPOL?"_withPOL":""),anaChoicesBacksc[i].Data(),
				ebinLow,ebinHigh,year==TString("2011-2012")?"0-59":"60-121").Data());
    
    Double_t AsymHold = 0.;
    Double_t AsymErrHold = 0.;
    if (infile.is_open()) {
      for (Int_t j=0; j<3; j++) infile >> strHold >> AsymHold >> AsymErrHold;
    }
    AsymmDataBacksc.push_back(AsymHold);
    AsymmDataErrBacksc.push_back(AsymErrHold);
    
    infile.close();
    
    infile.open(TString::Format("%s/Asymmetries/UNBLINDED_%s%s_OctetAsymmetries_AnaCh%s_%i-%i_Octets_%s.txt",
				getenv("ANALYSIS_RESULTS"),"AllCorr",
				(withPOL?"_withPOL":""),anaChoicesBacksc[i].Data(),
				ebinLow,ebinHigh,year==TString("2011-2012")?"0-59":"60-121").Data());
    
    AsymHold = 0.;
    AsymErrHold = 0.;
    if (infile.is_open()) {
      for (Int_t j=0; j<3; j++) infile >> strHold >> AsymHold >> AsymErrHold;
    }
    AsymmCorrDataBacksc.push_back(AsymHold);
    AsymmCorrDataErrBacksc.push_back(AsymErrHold);
    
    infile.close();
    
  }

  year = "2012-2013";
  

  std::vector <TString> anaChoicesALL_2012;
  std::vector <Double_t> indicesALL_2012;
  std::vector <Double_t> AsymmDataALL_2012;
  std::vector <Double_t> AsymmDataErrALL_2012; //Statistical Errors
  std::vector <Double_t> AsymmCorrDataALL_2012;
  std::vector <Double_t> AsymmCorrDataErrALL_2012; //Statistical Errors
 
  std::vector <TString> anaChoicesT0_2012;
  std::vector <Double_t> indicesT0_2012;
  std::vector <Double_t> AsymmDataT0_2012;
  std::vector <Double_t> AsymmDataErrT0_2012; //Statistical Errors
  std::vector <Double_t> AsymmCorrDataT0_2012;
  std::vector <Double_t> AsymmCorrDataErrT0_2012; //Statistical Errors
  

  std::vector <TString> anaChoicesBacksc_2012;
  std::vector <Double_t> indicesBacksc_2012;
  std::vector <Double_t> AsymmDataBacksc_2012;
  std::vector <Double_t> AsymmDataErrBacksc_2012; //Statistical Errors
  std::vector <Double_t> AsymmCorrDataBacksc_2012;
  std::vector <Double_t> AsymmCorrDataErrBacksc_2012; //Statistical Errors
  std::vector <Double_t> AsymmSimBacksc_2012;
  std::vector <Double_t> AsymmSimErrBacksc_2012; //Statistical Errors

  //Fill anaChoices
  for (Int_t i=0; i< sizeof(anChT0)/sizeof(TString); i++) 
    { anaChoicesT0_2012.push_back(anChT0[i]); indicesT0_2012.push_back(i+1);
      anaChoicesALL_2012.push_back(anChALL[i]); indicesALL_2012.push_back(i+1);}

  for (Int_t i=0; i< sizeof(anChBacksc)/sizeof(TString); i++) 
    { anaChoicesBacksc_2012.push_back(anChBacksc[i]); indicesBacksc_2012.push_back(i+1); }
  
  
  //Read in Asymmetries
  ifstream infile;
  std::string strHold;
  
  for (Int_t i=0; i<anaChoicesT0_2012.size(); i++) {
    infile.open(TString::Format("%s/Asymmetries/UNBLINDED_%s%s_OctetAsymmetries_AnaCh%s_%i-%i_Octets_%s.txt",
				getenv("ANALYSIS_RESULTS"),corrections.Data(),
				(withPOL?"_withPOL":""),anaChoicesT0_2012[i].Data(),
				ebinLow,ebinHigh,year==TString("2011-2012")?"0-59":"60-121").Data());
    
    Double_t AsymHold = 0.;
    Double_t AsymErrHold = 0.;
    if (infile.is_open()) {
      for (Int_t j=0; j<3; j++) infile >> strHold >> AsymHold >> AsymErrHold;
    }
    AsymmDataT0_2012.push_back(AsymHold);
    AsymmDataErrT0_2012.push_back(AsymErrHold);
    
    infile.close();
    
    infile.open(TString::Format("%s/Asymmetries/UNBLINDED_%s%s_OctetAsymmetries_AnaCh%s_%i-%i_Octets_%s.txt",
				getenv("ANALYSIS_RESULTS"),"AllCorr",
				(withPOL?"_withPOL":""),anaChoicesT0_2012[i].Data(),
				ebinLow,ebinHigh,year==TString("2011-2012")?"0-59":"60-121").Data());
    
    AsymHold = 0.;
    AsymErrHold = 0.;
    if (infile.is_open()) {
      for (Int_t j=0; j<3; j++) infile >> strHold >> AsymHold >> AsymErrHold;
    }
    AsymmCorrDataT0_2012.push_back(AsymHold);
    AsymmCorrDataErrT0_2012.push_back(AsymErrHold);
    
    infile.close();
    
    infile.open(TString::Format("%s/Asymmetries/UNBLINDED_%s%s_OctetAsymmetries_AnaCh%s_%i-%i_Octets_%s.txt",
				getenv("ANALYSIS_RESULTS"),corrections.Data(),(withPOL?"_withPOL":""),
				anaChoicesALL_2012[i].Data(),ebinLow,ebinHigh,
				year==TString("2011-2012")?"0-59":"60-121").Data() );
    
    AsymHold = 0.;
    AsymErrHold = 0.;
    if (infile.is_open()) {
      for (Int_t j=0; j<3; j++) infile >> strHold >> AsymHold >> AsymErrHold;
    }
    AsymmDataALL_2012.push_back(AsymHold);
    AsymmDataErrALL_2012.push_back(AsymErrHold);
    
    infile.close();

    infile.open(TString::Format("%s/Asymmetries/UNBLINDED_%s%s_OctetAsymmetries_AnaCh%s_%i-%i_Octets_%s.txt",
				getenv("ANALYSIS_RESULTS"),"AllCorr",(withPOL?"_withPOL":""),
				anaChoicesALL_2012[i].Data(),ebinLow,ebinHigh,
				year==TString("2011-2012")?"0-59":"60-121").Data() );
    
    AsymHold = 0.;
    AsymErrHold = 0.;
    if (infile.is_open()) {
      for (Int_t j=0; j<3; j++) infile >> strHold >> AsymHold >> AsymErrHold;
    }
    AsymmCorrDataALL_2012.push_back(AsymHold);
    AsymmCorrDataErrALL_2012.push_back(AsymErrHold);
    
    infile.close();
  }
  
  for (Int_t i=0; i<anaChoicesBacksc_2012.size(); i++) {
    infile.open(TString::Format("%s/Asymmetries/UNBLINDED_%s%s_OctetAsymmetries_AnaCh%s_%i-%i_Octets_%s.txt",
				getenv("ANALYSIS_RESULTS"),corrections.Data(),
				(withPOL?"_withPOL":""),anaChoicesBacksc_2012[i].Data(),
				ebinLow,ebinHigh,year==TString("2011-2012")?"0-59":"60-121").Data());
    
    Double_t AsymHold = 0.;
    Double_t AsymErrHold = 0.;
    if (infile.is_open()) {
      for (Int_t j=0; j<3; j++) infile >> strHold >> AsymHold >> AsymErrHold;
    }
    AsymmDataBacksc_2012.push_back(AsymHold);
    AsymmDataErrBacksc_2012.push_back(AsymErrHold);
    
    infile.close();
    
    infile.open(TString::Format("%s/Asymmetries/UNBLINDED_%s%s_OctetAsymmetries_AnaCh%s_%i-%i_Octets_%s.txt",
				getenv("ANALYSIS_RESULTS"),"AllCorr",
				(withPOL?"_withPOL":""),anaChoicesBacksc_2012[i].Data(),
				ebinLow,ebinHigh,year==TString("2011-2012")?"0-59":"60-121").Data());
    
    AsymHold = 0.;
    AsymErrHold = 0.;
    if (infile.is_open()) {
      for (Int_t j=0; j<3; j++) infile >> strHold >> AsymHold >> AsymErrHold;
    }
    AsymmCorrDataBacksc_2012.push_back(AsymHold);
    AsymmCorrDataErrBacksc_2012.push_back(AsymErrHold);
    
    infile.close();
    
  }
      

  Double_t xerr[10] = {0.};

  //cout << anaChoices.size() << " " << anaChoices[0] << " " << AsymmData[0] << " " << xerr[0] << " " << AsymmDataErr[0] << endl;
    

  gStyle->SetLegendBorderSize(0);
  gStyle->SetFillStyle(0);
  //gStyle->SetPadTopMargin(0.01);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetGridStyle(8);
  gStyle->SetGridColor(kBlack);
  gStyle->SetTitleSize(0.08,"t");
  gStyle->SetTitleAlign(23);

  
  // First make one plot with same scale for all event types
  TCanvas *c1 = new TCanvas("c1","demo bin labels",900,900);
  c1->Divide(1,2);
  c1->SetLeftMargin(0.15);
  c1->SetBottomMargin(0.15);

  c1->cd(1);

  gPad->SetGridy();
  gPad->SetTopMargin(0.0);

  TH1F * dataALL = new TH1F("dataALL","",AsymmDataALL.size(),0.5,AsymmDataALL.size()+0.5);
  dataALL->SetMarkerColor(color0);
  dataALL->SetLineColor(color0);
  dataALL->SetLineWidth(3);
  dataALL->SetMarkerStyle(markerUncorr0);
  dataALL->GetYaxis()->SetNdivisions(512);
  //dataALL->SetCanExtend(TH1::kAllAxes);
  dataALL->SetStats(0);

  for (Int_t i=0;i<AsymmDataT0.size();++i) {
    dataALL->SetBinContent(i+1,AsymmDataALL[i]);
    dataALL->SetBinError(i+1,AsymmDataErrALL[i]);
  }
  //dataALL->LabelsDeflate("X");
  //dataALL->LabelsDeflate("Y");
  //dataALL->LabelsOption("v");
  dataALL->GetXaxis()->SetBinLabel(1,"0,1,2*,3*");//A
  dataALL->GetXaxis()->SetBinLabel(2,"0,1");//B
  dataALL->GetXaxis()->SetBinLabel(3,"0,1,2,3");//C
  dataALL->GetXaxis()->SetBinLabel(4,"0");//D
  dataALL->GetXaxis()->SetBinLabel(5,"1");
  dataALL->GetXaxis()->SetBinLabel(6,"2*,3*");
  dataALL->GetXaxis()->SetBinLabel(7,"2,3");
  dataALL->GetXaxis()->SetBinLabel(8,"2");
  dataALL->GetXaxis()->SetBinLabel(9,"3");
  dataALL->GetXaxis()->SetLabelSize(0.06);
  dataALL->GetXaxis()->SetTitle("Event Types Included");
  dataALL->GetYaxis()->SetTitle("Asymmetry");
  dataALL->GetXaxis()->SetTitleSize(0.05);
  dataALL->GetYaxis()->SetTitleSize(0.05);
  dataALL->GetXaxis()->CenterTitle();
  dataALL->GetYaxis()->CenterTitle();
  dataALL->SetMinimum(yminBS);
  dataALL->SetMaximum(ymaxBS);

  dataALL->Draw("EX0");

  TLine *line1 = new TLine(dataALL->GetXaxis()->GetXmin(),-0.1184,dataALL->GetXaxis()->GetXmax(),-0.1184);
  line1->SetLineStyle(7);
  line1->SetLineWidth(2);
  if (corrections==TString("AllCorr")) line1->Draw("SAME");

  TH1F * corrDataALL = new TH1F("corrDataALL","",AsymmCorrDataALL.size(),0.5,AsymmCorrDataALL.size()+0.5);
  corrDataALL->SetMarkerColor(color0);
  corrDataALL->SetLineColor(color0);
  corrDataALL->SetLineWidth(3);
  corrDataALL->SetMarkerStyle(markerCorr0);
  //corrDataALL->GetXaxis()->SetNdivisions();
  //corrDataALL->SetCanExtend(TH1::kAllAxes);
  corrDataALL->SetStats(0);

  for (Int_t i=0;i<AsymmCorrDataT0.size();++i) {
    corrDataALL->SetBinContent(i+1,AsymmCorrDataALL[i]);
    corrDataALL->SetBinError(i+1,AsymmCorrDataErrALL[i]);
  }

  if (CorrAndUnCorr) corrDataALL->Draw("SAME EX0");
  
  
  TLegend *leg = new TLegend(0.2,0.65,0.5,0.9);
  leg->SetHeader("2011-2012"); // option "C" allows to center the header
  TLegendEntry *header = (TLegendEntry*)leg->GetListOfPrimitives()->First();
  header->SetTextSize(0.05);
  leg->AddEntry(dataALL,"Uncorrected Data","p");
  leg->AddEntry(corrDataALL,"Corrected Data","p");
  leg->Draw();

  c1->Print(TString::Format("CorrVsUncorr_singleAxis_%s%s.pdf",year.Data(),(color?"_color":"")));

  /// 2012-2013

  // First make one plot with same scale for all event types
  c1->cd(2);
  gPad->SetGridy();
  gPad->SetTopMargin(0.0);


  TH1F * dataALL_2012 = new TH1F("dataALL_2012","",AsymmDataALL_2012.size(),0.5,AsymmDataALL_2012.size()+0.5);
  dataALL_2012->SetMarkerColor(color0);
  dataALL_2012->SetLineColor(color0);
  dataALL_2012->SetLineWidth(3);
  dataALL_2012->SetMarkerStyle(markerUncorr0);
  dataALL_2012->GetYaxis()->SetNdivisions(512);
  //dataALL_2012->SetCanExtend(TH1::kAllAxes);
  dataALL_2012->SetStats(0);

  for (Int_t i=0;i<AsymmDataT0.size();++i) {
    dataALL_2012->SetBinContent(i+1,AsymmDataALL_2012[i]);
    dataALL_2012->SetBinError(i+1,AsymmDataErrALL_2012[i]);
  }
  //dataALL_2012->LabelsDeflate("X");
  //dataALL_2012->LabelsDeflate("Y");
  //dataALL_2012->LabelsOption("v");
  dataALL_2012->GetXaxis()->SetBinLabel(1,"0,1,2*,3*");//A
  dataALL_2012->GetXaxis()->SetBinLabel(2,"0,1");//B
  dataALL_2012->GetXaxis()->SetBinLabel(3,"0,1,2,3");//C
  dataALL_2012->GetXaxis()->SetBinLabel(4,"0");//D
  dataALL_2012->GetXaxis()->SetBinLabel(5,"1");
  dataALL_2012->GetXaxis()->SetBinLabel(6,"2*,3*");
  dataALL_2012->GetXaxis()->SetBinLabel(7,"2,3");
  dataALL_2012->GetXaxis()->SetBinLabel(8,"2");
  dataALL_2012->GetXaxis()->SetBinLabel(9,"3");
  dataALL_2012->GetXaxis()->SetLabelSize(0.06);
  dataALL_2012->GetXaxis()->SetTitle("Event Types Included");
  dataALL_2012->GetYaxis()->SetTitle("Asymmetry");
  dataALL_2012->GetXaxis()->SetTitleSize(0.05);
  dataALL_2012->GetYaxis()->SetTitleSize(0.05);
  dataALL_2012->GetXaxis()->CenterTitle();
  dataALL_2012->GetYaxis()->CenterTitle();
  dataALL_2012->SetMinimum(yminBS);
  dataALL_2012->SetMaximum(ymaxBS);

  dataALL_2012->Draw("EX0");

  TLine *line1_2012 = new TLine(dataALL_2012->GetXaxis()->GetXmin(),-0.1184,dataALL_2012->GetXaxis()->GetXmax(),-0.1184);
  line1_2012->SetLineStyle(7);
  line1_2012->SetLineWidth(2);
  if (corrections==TString("AllCorr")) line1->Draw("SAME");

  TH1F * corrDataALL_2012 = new TH1F("corrDataALL_2012","",AsymmCorrDataALL_2012.size(),0.5,AsymmCorrDataALL_2012.size()+0.5);
  corrDataALL_2012->SetMarkerColor(color0);
  corrDataALL_2012->SetLineColor(color0);
  corrDataALL_2012->SetLineWidth(3);
  corrDataALL_2012->SetMarkerStyle(markerCorr0);
  //corrDataALL_2012->GetXaxis()->SetNdivisions();
  //corrDataALL_2012->SetCanExtend(TH1::kAllAxes);
  corrDataALL_2012->SetStats(0);

  for (Int_t i=0;i<AsymmCorrDataT0.size();++i) {
    corrDataALL_2012->SetBinContent(i+1,AsymmCorrDataALL_2012[i]);
    corrDataALL_2012->SetBinError(i+1,AsymmCorrDataErrALL_2012[i]);
  }

  if (CorrAndUnCorr) corrDataALL_2012->Draw("SAME EX0");
  
  
  TLegend *leg_2012 = new TLegend(0.2,0.65,0.5,0.9);
  leg_2012->SetHeader("2012-2013"); // option "C" allows to center the header
  header = (TLegendEntry*)leg->GetListOfPrimitives()->First();
  header->SetTextSize(0.05);
  leg_2012->AddEntry(dataALL_2012,"Uncorrected Data","p");
  leg_2012->AddEntry(corrDataALL_2012,"Corrected Data","p");
  leg_2012->Draw();

  c1->Print(TString::Format("CorrVsUncorr_singleAxis%s.pdf",(color?"_color":"")));

  ///////////////// Multiple Axes /////////////////////////

  TCanvas *c2 = new TCanvas("c2","demo bin labels",900,600);
  c2->SetLeftMargin(0.15);
  c2->SetBottomMargin(0.15);

  TPad *pad1 = new TPad("pad1","",0,0,1,1);
  TPad *pad2 = new TPad("pad2","",0,0,1,1);
  pad2->SetFillStyle(4000);
  pad2->SetFrameFillStyle(0);

  TH1F * data = new TH1F("data","",AsymmDataT0.size(),0.5,AsymmDataT0.size()+0.5);
  data->SetMarkerColor(color0);
  data->SetLineColor(color0);
  data->SetLineWidth(3);
  data->SetMarkerStyle(markerUncorr0);
  data->GetYaxis()->SetNdivisions(512);
  //data->SetCanExtend(TH1::kAllAxes);
  data->SetStats(0);

  for (Int_t i=0;i<AsymmDataT0.size();++i) {
    data->SetBinContent(i+1,AsymmDataT0[i]);
    data->SetBinError(i+1,AsymmDataErrT0[i]);
  }
  //data->LabelsDeflate("X");
  //data->LabelsDeflate("Y");
  //data->LabelsOption("v");
  data->GetXaxis()->SetBinLabel(1,"0,1,2*,3*");//A
  data->GetXaxis()->SetBinLabel(2,"0,1");//B
  data->GetXaxis()->SetBinLabel(3,"0,1,2,3");//C
  data->GetXaxis()->SetBinLabel(4,"0");//D
  data->GetXaxis()->SetBinLabel(5,"1");
  data->GetXaxis()->SetBinLabel(6,"2*,3*");
  data->GetXaxis()->SetBinLabel(7,"2,3");
  data->GetXaxis()->SetBinLabel(8,"2");
  data->GetXaxis()->SetBinLabel(9,"3");
  data->GetXaxis()->SetLabelSize(0.05);
  data->GetXaxis()->SetTitle("Event Types Included");
  data->GetYaxis()->SetTitle("Asymmetry");
  data->GetYaxis()->SetTitleOffset(1.2);
  data->GetXaxis()->SetTitleSize(0.05);
  data->GetYaxis()->SetTitleSize(0.05);
  data->GetXaxis()->CenterTitle();
  data->GetYaxis()->CenterTitle();
  data->SetMinimum(ymin0);
  data->SetMaximum(ymax0);


  TLine *line1 = new TLine(data->GetXaxis()->GetXmin(),-0.1184,data->GetXaxis()->GetXmax(),-0.1184);
  line1->SetLineStyle(7);
  line1->SetLineWidth(2);
  if (corrections==TString("AllCorr")) line1->Draw("SAME");

  TH1F * corrData = new TH1F("corrData","",AsymmCorrDataT0.size(),0.5,AsymmCorrDataT0.size()+0.5);
  corrData->SetMarkerColor(color0);
  corrData->SetLineColor(color0);
  corrData->SetLineWidth(3);
  corrData->SetMarkerStyle(markerCorr0);
  //corrData->GetXaxis()->SetNdivisions();
  //corrData->SetCanExtend(TH1::kAllAxes);
  corrData->SetStats(0);

  for (Int_t i=0;i<AsymmCorrDataT0.size();++i) {
    corrData->SetBinContent(i+1,AsymmCorrDataT0[i]);
    corrData->SetBinError(i+1,AsymmCorrDataErrT0[i]);
  }


  TH1F * dataBacksc = new TH1F("dataBacksc","",AsymmDataBacksc.size(),0.5,AsymmDataBacksc.size()+0.5);
  dataBacksc->SetMarkerColor(colorBS);
  dataBacksc->SetLineColor(colorBS);
  dataBacksc->SetLineWidth(3);
  dataBacksc->SetMarkerStyle(markerUncorrBS);
  dataBacksc->GetYaxis()->SetNdivisions(512);
  //dataBacksc->SetCanExtend(TH1::kAllAxes);
  dataBacksc->SetStats(0);

  for (Int_t i=0;i<AsymmDataBacksc.size();++i) {
    dataBacksc->SetBinContent(i+1,AsymmDataBacksc[i]);
    dataBacksc->SetBinError(i+1,AsymmDataErrBacksc[i]);
  }
  //dataBacksc->LabelsDeflate("X");
  //dataBacksc->LabelsDeflate("Y");
  //dataBacksc->LabelsOption("v");
  dataBacksc->GetXaxis()->SetBinLabel(1,"0,1,2*,3*");//A
  dataBacksc->GetXaxis()->SetBinLabel(2,"0,1");//B
  dataBacksc->GetXaxis()->SetBinLabel(3,"0,1,2,3");//C
  dataBacksc->GetXaxis()->SetBinLabel(4,"0");//D
  dataBacksc->GetXaxis()->SetBinLabel(5,"1");
  dataBacksc->GetXaxis()->SetBinLabel(6,"2*,3*");
  dataBacksc->GetXaxis()->SetBinLabel(7,"2,3");
  dataBacksc->GetXaxis()->SetBinLabel(8,"2");
  dataBacksc->GetXaxis()->SetBinLabel(9,"3");
  dataBacksc->GetXaxis()->SetLabelSize(0.05);
  dataBacksc->GetXaxis()->SetTitle("Event Types Included");
  //dataBacksc->GetYaxis()->SetTitle("Asymmetry");
  //dataBacksc->GetYaxis()->SetTitleOffset(2.);
  dataBacksc->GetXaxis()->SetTitleSize(0.05);
  dataBacksc->GetYaxis()->SetTitleSize(0.05);
  dataBacksc->GetXaxis()->CenterTitle();
  dataBacksc->GetYaxis()->CenterTitle();
  dataBacksc->SetMinimum(yminBS);
  dataBacksc->SetMaximum(ymaxBS);



  TH1F * corrDataBacksc = new TH1F("","CorrDataBacksc vs MC Blinded Asymmetries",AsymmCorrDataBacksc.size(),0.5,AsymmCorrDataBacksc.size()+0.5);
  corrDataBacksc->SetMarkerColor(colorBS);
  corrDataBacksc->SetLineColor(colorBS);
  corrDataBacksc->SetLineWidth(3);
  corrDataBacksc->SetMarkerStyle(markerCorrBS);
  //corrDataBacksc->GetXaxis()->SetNdivisions();
  //corrDataBacksc->SetCanExtend(TH1::kAllAxes);
  corrDataBacksc->SetStats(0);

  for (Int_t i=0;i<AsymmCorrDataBacksc.size();++i) {
    corrDataBacksc->SetBinContent(i+1,AsymmCorrDataBacksc[i]);
    corrDataBacksc->SetBinError(i+1,AsymmCorrDataErrBacksc[i]);
  }

  if (CorrAndUnCorr) corrDataBacksc->Draw("SAME EX0");

  TH1F * simBacksc = new TH1F("","Sim vs MC Blinded Asymmetries",AsymmSimBacksc.size(),0.5,AsymmSimBacksc.size()+0.5);
  simBacksc->SetMarkerColor(kRed);
  simBacksc->SetLineColor(kRed);
  simBacksc->SetLineWidth(3);
  simBacksc->SetMarkerStyle(kFullSquare);
  //simBacksc->GetXaxis()->SetNdivisions();
  //simBacksc->SetCanExtend(TH1::kAllAxes);
  simBacksc->SetStats(0);

  for (Int_t i=0;i<AsymmSimBacksc.size();++i) {
    simBacksc->SetBinContent(i+1,AsymmSimBacksc[i]);
    simBacksc->SetBinError(i+1,AsymmSimErrBacksc[i]);
  }
  
  if (withSim) simBacksc->Draw("SAME EX0");

  pad1->Draw();
  pad1->SetGridy();
  pad1->cd();

  data->GetYaxis()->SetAxisColor(color0);
  data->GetYaxis()->SetLabelColor(color0);
  data->Draw("EX0");
  corrData->Draw("EX0 SAME");

  pad2->Draw();
  pad2->cd();
  dataBacksc->GetYaxis()->SetAxisColor(colorBS);
  dataBacksc->GetYaxis()->SetLabelColor(colorBS);
  dataBacksc->Draw("EX0 Y+");
  corrDataBacksc->Draw("EX0 SAME");

  TLegend *leg1 = new TLegend(0.15,0.65,0.5,0.85);
  //leg->SetHeader("The Legend Title","C"); // option "C" allows to center the header
  leg1->AddEntry(data,"Uncorrected Data Left Scale","p");
  leg1->AddEntry(dataBacksc,"Uncorrected Data Right Scale","p");
  leg1->AddEntry(corrData,"Corrected Data Left Scale","p");
  leg1->AddEntry(corrDataBacksc,"Corrected Data Right Scale","p");
  leg1->Draw("SAME");

 
  // TLegend *leg = new TLegend(0.65,0.65,0.95,0.8);
  //leg->SetHeader("The Legend Title","C"); // option "C" allows to center the header
  //leg->AddEntry(dataBacksc,"Data");
  //leg->AddEntry(simBacksc,"Monte Carlo");
   //leg->AddEntry("gr","Graph with error bars","lep");
  //if (withSim) leg->Draw();

  c2->Print(TString::Format("CorrVsUncorr_doubleAxis_%s%s.pdf",year.Data(),(color?"_color":"")));

  
}
