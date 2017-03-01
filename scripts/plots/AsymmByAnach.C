

void AsymmByAnach(TString year, TString corrections, bool withPOL, Int_t ebinLow=220, Int_t ebinHigh=680) {

  const Int_t numAnaCh = 4;
  TString anCh[numAnaCh] = {"A","D","F","G"};//{"A","B","C","D","E","F","G","H","J","K"};

  std::vector <TString> anaChoices;
  std::vector <Double_t> indices;
  std::vector <Double_t> AsymmData;
  std::vector <Double_t> AsymmDataErr; //Statistical Errors
  std::vector <Double_t> AsymmSim;
  std::vector <Double_t> AsymmSimErr; //Statistical Errors

  //Fill anaChoices
  for (Int_t i=0; i< sizeof(anCh)/sizeof(TString); i++) 
    { anaChoices.push_back(anCh[i]); indices.push_back(i+1); }


  //Read in Asymmetries
  ifstream infile;
  std::string strHold; 

  for (Int_t i=0; i<anaChoices.size(); i++) {
    infile.open(TString::Format("%s/Asymmetries/%s%s_OctetAsymmetries_AnaCh%s_%i-%i_Octets_%s.txt",
				getenv("ANALYSIS_RESULTS"),corrections.Data(),
				(withPOL?"_withPOL":""),anaChoices[i].Data(),
				ebinLow,ebinHigh,year==TString("2011-2012")?"0-59":"60-121").Data());

    Double_t AsymHold = 0.;
    Double_t AsymErrHold = 0.;
    if (infile.is_open()) {
      for (Int_t j=0; j<3; j++) infile >> strHold >> AsymHold >> AsymErrHold;
    }
    AsymmData.push_back(AsymHold);
    AsymmDataErr.push_back(AsymErrHold);
    
    infile.close();

    infile.open(TString::Format("%s/Asymmetries/%s_OctetAsymmetries_AnaCh%s_%i-%i_Octets_%s.txt",
				getenv("SIM_ANALYSIS_RESULTS"),corrections.Data(),
				anaChoices[i].Data(),ebinLow,ebinHigh,
				year==TString("2011-2012")?"0-59":"60-121").Data() );

    AsymHold = 0.;
    AsymErrHold = 0.;
    if (infile.is_open()) {
      for (Int_t j=0; j<3; j++) infile >> strHold >> AsymHold >> AsymErrHold;
    }
    AsymmSim.push_back(AsymHold);
    AsymmSimErr.push_back(AsymErrHold);
    
    infile.close();
  }

  

  

  Double_t xerr[10] = {0.};

  cout << anaChoices.size() << " " << anaChoices[0] << " " << AsymmData[0] << " " << xerr[0] << " " << AsymmDataErr[0] << endl;
    
  TCanvas *c1 = new TCanvas("c1","demo bin labels",10,10,600,600);

  c1->SetLeftMargin(0.15);
  c1->SetBottomMargin(0.15);

  
  TH1F * data = new TH1F("data","Data vs MC Blinded Asymmetries",4,0.5,4.5);
  data->SetMarkerColor(kBlue);
  data->SetLineColor(kBlue);
  data->SetLineWidth(3);
  data->SetMarkerStyle(kFullCircle);
  //data->GetXaxis()->SetNdivisions();
  data->SetCanExtend(TH1::kAllAxes);
  data->SetStats(0);

  for (Int_t i=1;i<5;++i) {
    data->SetBinContent(i,AsymmData[i-1]);
    data->SetBinError(i,AsymmDataErr[i-1]);
  }
  //data->LabelsDeflate("X");
  //data->LabelsDeflate("Y");
  //data->LabelsOption("v");
  data->GetXaxis()->SetBinLabel(1,"A");
  data->GetXaxis()->SetBinLabel(2,"D");
  data->GetXaxis()->SetBinLabel(3,"F");
  data->GetXaxis()->SetBinLabel(4,"G");
  data->GetXaxis()->SetLabelSize(0.07);
  data->GetXaxis()->SetTitle("Analysis Choice");
  data->GetYaxis()->SetTitle("Asymmetry");
  data->GetYaxis()->SetTitleOffset(2.);
  data->GetXaxis()->SetTitleSize(0.05);
  data->GetXaxis()->CenterTitle();
  data->GetYaxis()->CenterTitle();
  data->Draw("EX0");

  TH1F * sim = new TH1F("","Sim vs MC Blinded Asymmetries",4,0.5,4.5);
  sim->SetMarkerColor(kRed);
  sim->SetLineColor(kRed);
  sim->SetLineWidth(3);
  sim->SetMarkerStyle(kFullSquare);
  //sim->GetXaxis()->SetNdivisions();
  sim->SetCanExtend(TH1::kAllAxes);
  sim->SetStats(0);

  for (Int_t i=1;i<5;++i) {
    sim->SetBinContent(i,AsymmSim[i-1]);
    sim->SetBinError(i,AsymmSimErr[i-1]);
  }
  
  sim->Draw("SAME EX0");
  gStyle->SetLegendBorderSize(0);
  
  TLegend *leg = new TLegend(0.25,0.65,0.6,0.8);
  //leg->SetHeader("The Legend Title","C"); // option "C" allows to center the header
  leg->AddEntry(data,"Data");
  leg->AddEntry(sim,"Monte Carlo");
   //leg->AddEntry("gr","Graph with error bars","lep");
  leg->Draw();
  
  /*TGraphErrors *data0 = new TGraphErrors(numAnaCh,&indices[0],&AsymmData[0],xerr,&AsymmDataErr[0]);
  data0->SetMarkerColor(kBlue);
  data0->SetLineColor(kBlue);
  data0->SetMarkerStyle(kFullTriangleUp);
  data0->GetXaxis()->SetNdivisions(-400);
  data0->GetXaxis()->ChangeLabel(0,-1,-1,-1,-1,-1,"A");
  data0->GetXaxis()->ChangeLabel(1,-1,-1,-1,-1,-1,"D");
  data0->GetXaxis()->ChangeLabel(2,-1,-1,-1,-1,-1,"F");
  data0->GetXaxis()->ChangeLabel(3,-1,-1,-1,-1,-1,"G");
  data0->Draw("AP");

  TGraphErrors *sim0 = new TGraphErrors(numAnaCh,&indices[0],&AsymmSim[0],xerr,&AsymmSimErr[0]);
  sim0->SetMarkerColor(kRed);
  sim0->SetLineColor(kRed);
  sim0->SetMarkerStyle(kFullSquare);
  sim0->GetXaxis()->SetNdivisions(400);
  sim0->Draw("PSAME");


  TMultiGraph *mg0 = new TMultiGraph();
  mg0->Add(data0,"P");
  mg0->Add(sim0,"P");
  TAxis *a = mg0->GetXaxis();
  //a->SetNdivisions(-400);
  //mg0->GetXaxis()->ChangeLabel(0,-1,-1,-1,-1,-1,"A");
  //mg0->GetXaxis()->ChangeLabel(1,-1,-1,-1,-1,-1,"D");
  //mg0->GetXaxis()->ChangeLabel(2,-1,-1,-1,-1,-1,"F");
  //mg0->GetXaxis()->ChangeLabel(3,-1,-1,-1,-1,-1,"G");
  //mg0->SetMaximum(-0.116);
  //mg0->SetMinimum(-0.13);
  mg0->Draw("A");*/


  //Backscatters only
  /*TCanvas *c2 = new TCanvas("c2");

  TGraphErrors *dataBS = new TGraphErrors(2,&anaChoices[5],&AsymmData[5],xerr,&AsymmDataErr[5]);
  dataBS->SetMarkerColor(kBlue);
  dataBS->SetLineColor(kBlue);
  dataBS->SetMarkerStyle(kFullTriangleUp);


  TGraphErrors *simBS = new TGraphErrors(2,&anaChoices[5],&AsymmSim[5],xerr,&AsymmSimErr[5]);
  simBS->SetMarkerColor(kRed);
  simBS->SetLineColor(kRed);
  simBS->SetMarkerStyle(kFullSquare);


  TMultiGraph *mgBS = new TMultiGraph();
  mgBS->Add(dataBS,"P");
  mgBS->Add(simBS,"P");
  mgBS->SetMaximum(0.01);
  mgBS->SetMinimum(-0.15);

  mgBS->Draw("A");*/
}
