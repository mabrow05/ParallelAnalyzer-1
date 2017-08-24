

void betaSpectrum(int octetStart, int octetEnd) {
  bool color = false;
  
  gStyle->SetOptStat(0);
  //  gStyle->SetLabelSize(0.6);
  gStyle->SetTitleSize(0.08,"t");
  gStyle->SetTitleSize(0.08,"x");
  gStyle->SetTitleSize(0.08,"y");
  //gStyle->SetPadTopMargin(0.0);
  //gStyle->SetPadBottomMargin(0.3);
  //gStyle->SetPadLeftMargin(0.12);
  gStyle->SetLabelSize(0.07,"xyz");

  gStyle->SetLegendBorderSize(0);
  gStyle->SetFillStyle(0);

  TString normType = "ALL";

  Double_t normLow = 190.;
  Double_t normHigh = 740.;
  
  Double_t xAxisMax = 800.;

  Double_t minEnBin = 19; //Energy bins to integrate event types over
  Double_t maxEnBin = 73;

  //Storing event fractions for data, E/W and sfON/OFF
  Double_t t0E_sfON, t0E_sfOFF, t1E_sfON, t1E_sfOFF, t2E_sfON, t2E_sfOFF, t3E_sfON, t3E_sfOFF;
  Double_t t0W_sfON, t0W_sfOFF, t1W_sfON, t1W_sfOFF, t2W_sfON, t2W_sfOFF, t3W_sfON, t3W_sfOFF;
  t0E_sfON = t0E_sfOFF = t1E_sfON = t1E_sfOFF = t23E_sfON = t2E_sfOFF = t3E_sfOFF = 0;
  t0W_sfON = t0W_sfOFF = t1W_sfON = t1W_sfOFF = t23W_sfON = t2W_sfOFF = t3E_sfOFF = 0.;

  //Storing event fractions for sim, E/W and sfON/OFF
  Double_t sim_t0E_sfON, sim_t0E_sfOFF, sim_t1E_sfON, sim_t1E_sfOFF, sim_t2E_sfON, sim_t2E_sfOFF, sim_t3E_sfON, sim_t3E_sfOFF;
  Double_t sim_t0W_sfON, sim_t0W_sfOFF, sim_t1W_sfON, sim_t1W_sfOFF, sim_t2W_sfON, sim_t2W_sfOFF, sim_t3W_sfON, sim_t3W_sfOFF;
  sim_t0E_sfON = sim_t0E_sfOFF = sim_t1E_sfON = sim_t1E_sfOFF = sim_t23E_sfON = sim_t2E_sfOFF = sim_t3E_sfOFF = 0;
  sim_t0W_sfON = sim_t0W_sfOFF = sim_t1W_sfON = sim_t1W_sfOFF = sim_t23W_sfON = sim_t2W_sfOFF = sim_t3E_sfOFF = 0.;

  
  TFile *data_file = new TFile(TString::Format("../../Asymmetry/SpectraComparisons/Octets_%i-%i_DATA.root",octetStart,octetEnd),"READ");
  TFile *sim_file = new TFile(TString::Format("../../Asymmetry/SpectraComparisons/Octets_%i-%i_SIM.root",octetStart,octetEnd),"READ");


  //Load the histograms

  TH1D *dataALL = (TH1D*)data_file->Get("EreconALL");
  TH1D *data0 = (TH1D*)data_file->Get("Erecon0");
  TH1D *data1 = (TH1D*)data_file->Get("Erecon1");
  TH1D *data2 = (TH1D*)data_file->Get("Erecon2");
  TH1D *data3 = (TH1D*)data_file->Get("Erecon3");

  TH1D *BG_dataALL = (TH1D*)data_file->Get("BG_EreconALL");
  TH1D *BG_data0 = (TH1D*)data_file->Get("BG_Erecon0");
  TH1D *BG_data1 = (TH1D*)data_file->Get("BG_Erecon1");
  TH1D *BG_data2 = (TH1D*)data_file->Get("BG_Erecon2");
  TH1D *BG_data3 = (TH1D*)data_file->Get("BG_Erecon3");
  
  TH1D *simALL = (TH1D*)sim_file->Get("EreconALL");
  TH1D *sim0 = (TH1D*)sim_file->Get("Erecon0");
  TH1D *sim1 = (TH1D*)sim_file->Get("Erecon1");
  TH1D *sim2 = (TH1D*)sim_file->Get("Erecon2");
  TH1D *sim3 = (TH1D*)sim_file->Get("Erecon3");
  
  TH1D *dataALL_sfOFF_E = (TH1D*)data_file->Get("EreconALL_sfOFF_E");
  TH1D *data0_sfOFF_E = (TH1D*)data_file->Get("Erecon0_sfOFF_E");
  TH1D *data1_sfOFF_E = (TH1D*)data_file->Get("Erecon1_sfOFF_E");
  TH1D *data2_sfOFF_E = (TH1D*)data_file->Get("Erecon2_sfOFF_E");
  TH1D *data3_sfOFF_E = (TH1D*)data_file->Get("Erecon3_sfOFF_E");

  TH1D *BG_dataALL_sfOFF_E = (TH1D*)data_file->Get("BG_EreconALL_sfOFF_E");
  TH1D *BG_data0_sfOFF_E = (TH1D*)data_file->Get("BG_Erecon0_sfOFF_E");
  TH1D *BG_data1_sfOFF_E = (TH1D*)data_file->Get("BG_Erecon1_sfOFF_E");
  TH1D *BG_data2_sfOFF_E = (TH1D*)data_file->Get("BG_Erecon2_sfOFF_E");
  TH1D *BG_data3_sfOFF_E = (TH1D*)data_file->Get("BG_Erecon3_sfOFF_E");
  
  TH1D *simALL_sfOFF_E = (TH1D*)sim_file->Get("EreconALL_sfOFF_E");
  TH1D *sim0_sfOFF_E = (TH1D*)sim_file->Get("Erecon0_sfOFF_E");
  TH1D *sim1_sfOFF_E = (TH1D*)sim_file->Get("Erecon1_sfOFF_E");
  TH1D *sim2_sfOFF_E = (TH1D*)sim_file->Get("Erecon2_sfOFF_E");
  TH1D *sim3_sfOFF_E = (TH1D*)sim_file->Get("Erecon3_sfOFF_E");

  TH1D *dataALL_sfOFF_W = (TH1D*)data_file->Get("EreconALL_sfOFF_W");
  TH1D *data0_sfOFF_W = (TH1D*)data_file->Get("Erecon0_sfOFF_W");
  TH1D *data1_sfOFF_W = (TH1D*)data_file->Get("Erecon1_sfOFF_W");
  TH1D *data2_sfOFF_W = (TH1D*)data_file->Get("Erecon2_sfOFF_W");
  TH1D *data3_sfOFF_W = (TH1D*)data_file->Get("Erecon3_sfOFF_W");

  TH1D *BG_dataALL_sfOFF_W = (TH1D*)data_file->Get("BG_EreconALL_sfOFF_W");
  TH1D *BG_data0_sfOFF_W = (TH1D*)data_file->Get("BG_Erecon0_sfOFF_W");
  TH1D *BG_data1_sfOFF_W = (TH1D*)data_file->Get("BG_Erecon1_sfOFF_W");
  TH1D *BG_data2_sfOFF_W = (TH1D*)data_file->Get("BG_Erecon2_sfOFF_W");
  TH1D *BG_data3_sfOFF_W = (TH1D*)data_file->Get("BG_Erecon3_sfOFF_W");
  
  TH1D *simALL_sfOFF_W = (TH1D*)sim_file->Get("EreconALL_sfOFF_W");
  TH1D *sim0_sfOFF_W = (TH1D*)sim_file->Get("Erecon0_sfOFF_W");
  TH1D *sim1_sfOFF_W = (TH1D*)sim_file->Get("Erecon1_sfOFF_W");
  TH1D *sim2_sfOFF_W = (TH1D*)sim_file->Get("Erecon2_sfOFF_W");
  TH1D *sim3_sfOFF_W = (TH1D*)sim_file->Get("Erecon3_sfOFF_W");
  
  TH1D *dataALL_sfON_E = (TH1D*)data_file->Get("EreconALL_sfON_E");
  TH1D *data0_sfON_E = (TH1D*)data_file->Get("Erecon0_sfON_E");
  TH1D *data1_sfON_E = (TH1D*)data_file->Get("Erecon1_sfON_E");
  TH1D *data2_sfON_E = (TH1D*)data_file->Get("Erecon2_sfON_E");
  TH1D *data3_sfON_E = (TH1D*)data_file->Get("Erecon3_sfON_E");

  TH1D *BG_dataALL_sfON_E = (TH1D*)data_file->Get("BG_EreconALL_sfON_E");
  TH1D *BG_data0_sfON_E = (TH1D*)data_file->Get("BG_Erecon0_sfON_E");
  TH1D *BG_data1_sfON_E = (TH1D*)data_file->Get("BG_Erecon1_sfON_E");
  TH1D *BG_data2_sfON_E = (TH1D*)data_file->Get("BG_Erecon2_sfON_E");
  TH1D *BG_data3_sfON_E = (TH1D*)data_file->Get("BG_Erecon3_sfON_E");
  
  TH1D *simALL_sfON_E = (TH1D*)sim_file->Get("EreconALL_sfON_E");
  TH1D *sim0_sfON_E = (TH1D*)sim_file->Get("Erecon0_sfON_E");
  TH1D *sim1_sfON_E = (TH1D*)sim_file->Get("Erecon1_sfON_E");
  TH1D *sim2_sfON_E = (TH1D*)sim_file->Get("Erecon2_sfON_E");
  TH1D *sim3_sfON_E = (TH1D*)sim_file->Get("Erecon3_sfON_E");

  TH1D *dataALL_sfON_W = (TH1D*)data_file->Get("EreconALL_sfON_W");
  TH1D *data0_sfON_W = (TH1D*)data_file->Get("Erecon0_sfON_W");
  TH1D *data1_sfON_W = (TH1D*)data_file->Get("Erecon1_sfON_W");
  TH1D *data2_sfON_W = (TH1D*)data_file->Get("Erecon2_sfON_W");
  TH1D *data3_sfON_W = (TH1D*)data_file->Get("Erecon3_sfON_W");

  TH1D *BG_dataALL_sfON_W = (TH1D*)data_file->Get("BG_EreconALL_sfON_W");
  TH1D *BG_data0_sfON_W = (TH1D*)data_file->Get("BG_Erecon0_sfON_W");
  TH1D *BG_data1_sfON_W = (TH1D*)data_file->Get("BG_Erecon1_sfON_W");
  TH1D *BG_data2_sfON_W = (TH1D*)data_file->Get("BG_Erecon2_sfON_W");
  TH1D *BG_data3_sfON_W = (TH1D*)data_file->Get("BG_Erecon3_sfON_W");
  
  TH1D *simALL_sfON_W = (TH1D*)sim_file->Get("EreconALL_sfON_W");
  TH1D *sim0_sfON_W = (TH1D*)sim_file->Get("Erecon0_sfON_W");
  TH1D *sim1_sfON_W = (TH1D*)sim_file->Get("Erecon1_sfON_W");
  TH1D *sim2_sfON_W = (TH1D*)sim_file->Get("Erecon2_sfON_W");
  TH1D *sim3_sfON_W = (TH1D*)sim_file->Get("Erecon3_sfON_W");
  

  t0E_sfOFF = data0_sfOFF_E->Integral(minEnBin,maxEnBin);
  t0E_sfON = data0_sfON_E->Integral(minEnBin,maxEnBin);
  t1E_sfOFF = data1_sfOFF_E->Integral(minEnBin,maxEnBin);
  t1E_sfON = data1_sfON_E->Integral(minEnBin,maxEnBin);
  t2E_sfOFF = data2_sfOFF_E->Integral(minEnBin,maxEnBin);
  t2E_sfON = data2_sfON_E->Integral(minEnBin,maxEnBin);
  t3E_sfOFF = data3_sfOFF_E->Integral(minEnBin,maxEnBin);
  t3E_sfON = data3_sfON_E->Integral(minEnBin,maxEnBin);
  t0W_sfOFF = data0_sfOFF_W->Integral(minEnBin,maxEnBin);
  t0W_sfON = data0_sfON_W->Integral(minEnBin,maxEnBin);
  t1W_sfOFF = data1_sfOFF_W->Integral(minEnBin,maxEnBin);
  t1W_sfON = data1_sfON_W->Integral(minEnBin,maxEnBin);
  t2W_sfOFF = data2_sfOFF_W->Integral(minEnBin,maxEnBin);
  t2W_sfON = data2_sfON_W->Integral(minEnBin,maxEnBin);
  t3W_sfOFF = data3_sfOFF_W->Integral(minEnBin,maxEnBin);
  t3W_sfON = data3_sfON_W->Integral(minEnBin,maxEnBin);
  

  double totalE_sfOFF = t0E_sfOFF;// + t1E_sfOFF + t2E_sfOFF + t3E_sfOFF;
  double totalE_sfON = t0E_sfON;// + t1E_sfON + t2E_sfON + t3E_sfON;
  double totalW_sfOFF = t0W_sfOFF;// + t1W_sfOFF + t2W_sfOFF + t3W_sfOFF;
  double totalW_sfON = t0W_sfON;// + t1W_sfON + t2W_sfON + t3W_sfON;

  double fracData1 = 100.*data1->Integral(minEnBin,maxEnBin)/data0->Integral(minEnBin,maxEnBin);
  double fracData2 = 100.*data2->Integral(minEnBin,maxEnBin)/data0->Integral(minEnBin,maxEnBin);
  double fracData3 = 100.*data3->Integral(minEnBin,maxEnBin)/data0->Integral(minEnBin,maxEnBin);  
  double fracSim1 = 100.*sim1->Integral(minEnBin,maxEnBin)/sim0->Integral(minEnBin,maxEnBin);
  double fracSim2 = 100.*sim2->Integral(minEnBin,maxEnBin)/sim0->Integral(minEnBin,maxEnBin);
  double fracSim3 = 100.*sim3->Integral(minEnBin,maxEnBin)/sim0->Integral(minEnBin,maxEnBin);  
  
  std::cout << std::fixed;
  std::cout << std::setprecision(2);
  std::cout << "********************************************\n"
	    << "          Event ratios (220-670 keV)\n\n"
	    << "\t\tsfON\t\t\tsfOFF\n"
	    << "Type 0E:\t" << t0E_sfON/totalE_sfON << "\t\t" << t0E_sfOFF/totalE_sfOFF << "\n"
	    << "Type 0W:\t" << t0W_sfON/totalW_sfON << "\t\t" << t0W_sfOFF/totalW_sfOFF << "\n\n"
	    << "Type 1E:\t" << t1E_sfON/totalE_sfON << "\t\t" << t1E_sfOFF/totalE_sfOFF << "\n"
	    << "Type 1W:\t" << t1W_sfON/totalW_sfON << "\t\t" << t1W_sfOFF/totalW_sfOFF << "\n\n"
	    << "Type 2E:\t" << t2E_sfON/totalE_sfON << "\t\t" << t2E_sfOFF/totalE_sfOFF << "\n"
	    << "Type 2W:\t" << t2W_sfON/totalW_sfON << "\t\t" << t2W_sfOFF/totalW_sfOFF << "\n\n"
	    << "Type 3E:\t" << t3E_sfON/totalE_sfON << "\t\t" << t3E_sfOFF/totalE_sfOFF << "\n"
	    << "Type 3W:\t" << t3W_sfON/totalW_sfON << "\t\t" << t3W_sfOFF/totalW_sfOFF << "\n\n"
	    << "********************************************\n\n"
	    << "          Total Fractions (as % of type0)\n"
	    << "Type\t\t1\t2\t3\n"
	    << "Data\t\t"<<fracData1<<"\t"<<fracData2<<"\t"<<fracData3<<"\n"
	    << "Sim\t\t"<<fracSim1<<"\t"<<fracSim2<<"\t"<<fracSim3<<"\n"
	    << "% Diff\t\t"
	    << (fracSim1-fracData1)/fracData1*100.<<"\t"
	    << (fracSim2-fracData2)/fracData2*100.<<"\t"
	    << (fracSim3-fracData3)/fracData3*100.<<"\n\n***************************************\n\n";
 
    
  //Normalize
  Double_t normFactor = 1.;
  Double_t normFactor_sfOFF_E = 1.;
  Double_t normFactor_sfON_E = 1.;
  Double_t normFactor_sfOFF_W = 1.;
  Double_t normFactor_sfON_W = 1.;

  Double_t normFactor_ALL = 1.;
  Double_t normFactor_sfOFF_E_ALL = 1.;
  Double_t normFactor_sfON_E_ALL = 1.;
  Double_t normFactor_sfOFF_W_ALL = 1.;
  Double_t normFactor_sfON_W_ALL = 1.;

  Double_t normFactor_0 = 1.;
  Double_t normFactor_sfOFF_E_0 = 1.;
  Double_t normFactor_sfON_E_0 = 1.;
  Double_t normFactor_sfOFF_W_0 = 1.;
  Double_t normFactor_sfON_W_0 = 1.;

  Double_t normFactor_1 = 1.;
  Double_t normFactor_sfOFF_E_1 = 1.;
  Double_t normFactor_sfON_E_1 = 1.;
  Double_t normFactor_sfOFF_W_1 = 1.;
  Double_t normFactor_sfON_W_1 = 1.;
  
  Double_t normFactor_2 = 1.;
  Double_t normFactor_sfOFF_E_2 = 1.;
  Double_t normFactor_sfON_E_2 = 1.;
  Double_t normFactor_sfOFF_W_2 = 1.;
  Double_t normFactor_sfON_W_2 = 1.;
  
  Double_t normFactor_3 = 1.;
  Double_t normFactor_sfOFF_E_3 = 1.;
  Double_t normFactor_sfON_E_3 = 1.;
  Double_t normFactor_sfOFF_W_3 = 1.;
  Double_t normFactor_sfON_W_3 = 1.;



  
  normFactor_ALL = dataALL->Integral(dataALL->GetXaxis()->FindFixBin(normLow),dataALL->GetXaxis()->FindFixBin(normHigh))/simALL->Integral(simALL->GetXaxis()->FindFixBin(normLow),simALL->GetXaxis()->FindFixBin(normHigh));
  normFactor_sfOFF_E_ALL = dataALL_sfOFF_E->Integral(dataALL_sfOFF_E->GetXaxis()->FindFixBin(normLow),dataALL_sfOFF_E->GetXaxis()->FindFixBin(normHigh))/simALL_sfOFF_E->Integral(simALL_sfOFF_E->GetXaxis()->FindFixBin(normLow),simALL_sfOFF_E->GetXaxis()->FindFixBin(normHigh));
  normFactor_sfON_E_ALL = dataALL_sfON_E->Integral(dataALL_sfON_E->GetXaxis()->FindFixBin(normLow),dataALL_sfON_E->GetXaxis()->FindFixBin(normHigh))/simALL_sfON_E->Integral(simALL_sfON_E->GetXaxis()->FindFixBin(normLow),simALL_sfON_E->GetXaxis()->FindFixBin(normHigh));
  normFactor_sfOFF_W_ALL = dataALL_sfOFF_W->Integral(dataALL_sfOFF_W->GetXaxis()->FindFixBin(normLow),dataALL_sfOFF_W->GetXaxis()->FindFixBin(normHigh))/simALL_sfOFF_W->Integral(simALL_sfOFF_W->GetXaxis()->FindFixBin(normLow),simALL_sfOFF_W->GetXaxis()->FindFixBin(normHigh));
  normFactor_sfON_W_ALL = dataALL_sfON_W->Integral(dataALL_sfON_W->GetXaxis()->FindFixBin(normLow),dataALL_sfON_W->GetXaxis()->FindFixBin(normHigh))/simALL_sfON_W->Integral(simALL_sfON_W->GetXaxis()->FindFixBin(normLow),simALL_sfON_W->GetXaxis()->FindFixBin(normHigh));
  
  normFactor_0 = data0->Integral(data0->GetXaxis()->FindFixBin(normLow),data0->GetXaxis()->FindFixBin(normHigh))/sim0->Integral(sim0->GetXaxis()->FindFixBin(normLow),sim0->GetXaxis()->FindFixBin(normHigh));
  normFactor_sfOFF_E_0 = data0_sfOFF_E->Integral(data0_sfOFF_E->GetXaxis()->FindFixBin(normLow),data0_sfOFF_E->GetXaxis()->FindFixBin(normHigh))/sim0_sfOFF_E->Integral(sim0_sfOFF_E->GetXaxis()->FindFixBin(normLow),sim0_sfOFF_E->GetXaxis()->FindFixBin(normHigh));
  normFactor_sfON_E_0 = data0_sfON_E->Integral(data0_sfON_E->GetXaxis()->FindFixBin(normLow),data0_sfON_E->GetXaxis()->FindFixBin(normHigh))/sim0_sfON_E->Integral(sim0_sfON_E->GetXaxis()->FindFixBin(normLow),sim0_sfON_E->GetXaxis()->FindFixBin(normHigh));
  normFactor_sfOFF_W_0 = data0_sfOFF_W->Integral(data0_sfOFF_W->GetXaxis()->FindFixBin(normLow),data0_sfOFF_W->GetXaxis()->FindFixBin(normHigh))/sim0_sfOFF_W->Integral(sim0_sfOFF_W->GetXaxis()->FindFixBin(normLow),sim0_sfOFF_W->GetXaxis()->FindFixBin(normHigh));
  normFactor_sfON_W_0 = data0_sfON_W->Integral(data0_sfON_W->GetXaxis()->FindFixBin(normLow),data0_sfON_W->GetXaxis()->FindFixBin(normHigh))/sim0_sfON_W->Integral(sim0_sfON_W->GetXaxis()->FindFixBin(normLow),sim0_sfON_W->GetXaxis()->FindFixBin(normHigh));
  
  normFactor_1 = data1->Integral(data1->GetXaxis()->FindFixBin(normLow),data1->GetXaxis()->FindFixBin(normHigh))/sim1->Integral(sim1->GetXaxis()->FindFixBin(normLow),sim1->GetXaxis()->FindFixBin(normHigh));
  normFactor_sfOFF_E_1 = data1_sfOFF_E->Integral(data1_sfOFF_E->GetXaxis()->FindFixBin(normLow),data1_sfOFF_E->GetXaxis()->FindFixBin(normHigh))/sim1_sfOFF_E->Integral(sim1_sfOFF_E->GetXaxis()->FindFixBin(normLow),sim1_sfOFF_E->GetXaxis()->FindFixBin(normHigh));
  normFactor_sfON_E_1 = data1_sfON_E->Integral(data1_sfON_E->GetXaxis()->FindFixBin(normLow),data1_sfON_E->GetXaxis()->FindFixBin(normHigh))/sim1_sfON_E->Integral(sim1_sfON_E->GetXaxis()->FindFixBin(normLow),sim1_sfON_E->GetXaxis()->FindFixBin(normHigh));
  normFactor_sfOFF_W_1 = data1_sfOFF_W->Integral(data1_sfOFF_W->GetXaxis()->FindFixBin(normLow),data1_sfOFF_W->GetXaxis()->FindFixBin(normHigh))/sim1_sfOFF_W->Integral(sim1_sfOFF_W->GetXaxis()->FindFixBin(normLow),sim1_sfOFF_W->GetXaxis()->FindFixBin(normHigh));
  normFactor_sfON_W_1 = data1_sfON_W->Integral(data1_sfON_W->GetXaxis()->FindFixBin(normLow),data1_sfON_W->GetXaxis()->FindFixBin(normHigh))/sim1_sfON_W->Integral(sim1_sfON_W->GetXaxis()->FindFixBin(normLow),sim1_sfON_W->GetXaxis()->FindFixBin(normHigh));
  
  normFactor_2 = data2->Integral(data2->GetXaxis()->FindFixBin(normLow),data2->GetXaxis()->FindFixBin(normHigh))/sim2->Integral(sim2->GetXaxis()->FindFixBin(normLow),sim2->GetXaxis()->FindFixBin(normHigh));
  normFactor_sfOFF_E_2 = data2_sfOFF_E->Integral(data2_sfOFF_E->GetXaxis()->FindFixBin(normLow),data2_sfOFF_E->GetXaxis()->FindFixBin(normHigh))/sim2_sfOFF_E->Integral(sim2_sfOFF_E->GetXaxis()->FindFixBin(normLow),sim2_sfOFF_E->GetXaxis()->FindFixBin(normHigh));
  normFactor_sfON_E_2 = data2_sfON_E->Integral(data2_sfON_E->GetXaxis()->FindFixBin(normLow),data2_sfON_E->GetXaxis()->FindFixBin(normHigh))/sim2_sfON_E->Integral(sim2_sfON_E->GetXaxis()->FindFixBin(normLow),sim2_sfON_E->GetXaxis()->FindFixBin(normHigh));
  normFactor_sfOFF_W_2 = data2_sfOFF_W->Integral(data2_sfOFF_W->GetXaxis()->FindFixBin(normLow),data2_sfOFF_W->GetXaxis()->FindFixBin(normHigh))/sim2_sfOFF_W->Integral(sim2_sfOFF_W->GetXaxis()->FindFixBin(normLow),sim2_sfOFF_W->GetXaxis()->FindFixBin(normHigh));
  normFactor_sfON_W_2 = data2_sfON_W->Integral(data2_sfON_W->GetXaxis()->FindFixBin(normLow),data2_sfON_W->GetXaxis()->FindFixBin(normHigh))/sim2_sfON_W->Integral(sim2_sfON_W->GetXaxis()->FindFixBin(normLow),sim2_sfON_W->GetXaxis()->FindFixBin(normHigh));
  
  normFactor_3 = data3->Integral(data3->GetXaxis()->FindFixBin(normLow),data3->GetXaxis()->FindFixBin(normHigh))/sim3->Integral(sim3->GetXaxis()->FindFixBin(normLow),sim3->GetXaxis()->FindFixBin(normHigh));
  normFactor_sfOFF_E_3 = data3_sfOFF_E->Integral(data3_sfOFF_E->GetXaxis()->FindFixBin(normLow),data3_sfOFF_E->GetXaxis()->FindFixBin(normHigh))/sim3_sfOFF_E->Integral(sim3_sfOFF_E->GetXaxis()->FindFixBin(normLow),sim3_sfOFF_E->GetXaxis()->FindFixBin(normHigh));
  normFactor_sfON_E_3 = data3_sfON_E->Integral(data3_sfON_E->GetXaxis()->FindFixBin(normLow),data3_sfON_E->GetXaxis()->FindFixBin(normHigh))/sim3_sfON_E->Integral(sim3_sfON_E->GetXaxis()->FindFixBin(normLow),sim3_sfON_E->GetXaxis()->FindFixBin(normHigh));
  normFactor_sfOFF_W_3 = data3_sfOFF_W->Integral(data3_sfOFF_W->GetXaxis()->FindFixBin(normLow),data3_sfOFF_W->GetXaxis()->FindFixBin(normHigh))/sim3_sfOFF_W->Integral(sim3_sfOFF_W->GetXaxis()->FindFixBin(normLow),sim3_sfOFF_W->GetXaxis()->FindFixBin(normHigh));
  normFactor_sfON_W_3 = data3_sfON_W->Integral(data3_sfON_W->GetXaxis()->FindFixBin(normLow),data3_sfON_W->GetXaxis()->FindFixBin(normHigh))/sim3_sfON_W->Integral(sim3_sfON_W->GetXaxis()->FindFixBin(normLow),sim3_sfON_W->GetXaxis()->FindFixBin(normHigh));
  
  
  
  
  if (normType==TString("ALL")) {     
    normFactor = normFactor_ALL;
    normFactor_sfOFF_E = normFactor_sfOFF_E_ALL;
    normFactor_sfOFF_W = normFactor_sfOFF_W_ALL;
    normFactor_sfON_E = normFactor_sfON_E_ALL;
    normFactor_sfON_W = normFactor_sfON_W_ALL;
  }
  
  if (normType==TString("0")) {     
    normFactor = normFactor_0;
    normFactor_sfOFF_E = normFactor_sfOFF_E_0;
    normFactor_sfOFF_W = normFactor_sfOFF_W_0;
    normFactor_sfON_E = normFactor_sfON_E_0;
    normFactor_sfON_W = normFactor_sfON_W_0;
  }
  
  if (normType==TString("1")) {     
    normFactor = normFactor_1;
    normFactor_sfOFF_E = normFactor_sfOFF_E_1;
    normFactor_sfOFF_W = normFactor_sfOFF_W_1;
    normFactor_sfON_E = normFactor_sfON_E_1;
    normFactor_sfON_W = normFactor_sfON_W_1;
  }

  if (normType==TString("2")) {     
    normFactor = normFactor_2;
    normFactor_sfOFF_E = normFactor_sfOFF_E_2;
    normFactor_sfOFF_W = normFactor_sfOFF_W_2;
    normFactor_sfON_E = normFactor_sfON_E_2;
    normFactor_sfON_W = normFactor_sfON_W_2;
  }

  if (normType==TString("3")) {     
    normFactor = normFactor_3;
    normFactor_sfOFF_E = normFactor_sfOFF_E_3;
    normFactor_sfOFF_W = normFactor_sfOFF_W_3;
    normFactor_sfON_E = normFactor_sfON_E_3;
    normFactor_sfON_W = normFactor_sfON_W_3;
  }
  
  
  simALL->Scale(normFactor);
  sim0->Scale(normFactor);
  sim1->Scale(normFactor);
  sim2->Scale(normFactor);
  sim3->Scale(normFactor);

  simALL_sfOFF_E->Scale(normFactor_sfOFF_E);
  sim0_sfOFF_E->Scale(normFactor_sfOFF_E);
  sim1_sfOFF_E->Scale(normFactor_sfOFF_E);
  sim2_sfOFF_E->Scale(normFactor_sfOFF_E);
  sim3_sfOFF_E->Scale(normFactor_sfOFF_E);

  simALL_sfON_E->Scale(normFactor_sfON_E);
  sim0_sfON_E->Scale(normFactor_sfON_E);
  sim1_sfON_E->Scale(normFactor_sfON_E);
  sim2_sfON_E->Scale(normFactor_sfON_E);
  sim3_sfON_E->Scale(normFactor_sfON_E);

  simALL_sfOFF_W->Scale(normFactor_sfOFF_W);
  sim0_sfOFF_W->Scale(normFactor_sfOFF_W);
  sim1_sfOFF_W->Scale(normFactor_sfOFF_W);
  sim2_sfOFF_W->Scale(normFactor_sfOFF_W);
  sim3_sfOFF_W->Scale(normFactor_sfOFF_W);

  simALL_sfON_W->Scale(normFactor_sfON_W);
  sim0_sfON_W->Scale(normFactor_sfON_W);
  sim1_sfON_W->Scale(normFactor_sfON_W);
  sim2_sfON_W->Scale(normFactor_sfON_W);
  sim3_sfON_W->Scale(normFactor_sfON_W);
  

  ////////////////////// SUPER_SUMS ////////////////////////////////////////////

  //////// Super-Sum of all data ///////

  // All Types

  //Residuals...
  TCanvas *c1 = new TCanvas("c1","c1", 500, 600);
  //c1->Divide(1,2);
  //c1->cd(2);
  TPad *p1 = new TPad("p1","p1",0.,0.3,1.,1.);
  TPad *p2 = new TPad("p2","p2",0.,0.,1.,0.3);
  p2->SetBottomMargin(0.3);
  p1->Draw();
  p2->Draw();

  p2->cd();
  Int_t nBins = dataALL->GetNbinsX();
  Double_t Min = dataALL->GetXaxis()->GetBinLowEdge(dataALL->GetXaxis()->GetFirst());
  Double_t Max = dataALL->GetXaxis()->GetBinUpEdge(dataALL->GetXaxis()->GetLast());

  //residual
  TH1D *residALL = new TH1D("residALL","Residuals: MC-Data", nBins, Min, Max);
  residALL->Add(simALL,dataALL,1,-1);
  residALL->GetXaxis()->SetRangeUser(0., xAxisMax);
  residALL->GetXaxis()->SetTitleOffset(1.);
  residALL->GetYaxis()->SetTitleOffset(0.5);
  residALL->GetYaxis()->SetTitle("event rate (mHz/keV)");
  residALL->GetXaxis()->SetTitle("Energy (keV)");
  //residALL->SetMaximum(2.);
  //residALL->SetMinimum(-2.);
  residALL->SetLineWidth(2);
  residALL->SetMarkerColor(color?kBlue:kBlack);
  residALL->SetMarkerStyle(34);
  residALL->SetMarkerSize(1.);

  TH1D *perc_residALL = new TH1D("perc_residALL","Fractional Residuals: (MC-Data)/MC", nBins, Min, Max);
  perc_residALL->Add(simALL,dataALL,1,-1);
  perc_residALL->Divide(perc_residALL,dataALL,1.,1.);
  perc_residALL->GetXaxis()->SetRangeUser(0., xAxisMax);
  //perc_residALL->GetYaxis()->SetTitle("event rate (mHz/keV)");
  perc_residALL->SetMaximum(.1);
  perc_residALL->SetMinimum(-.1);
  perc_residALL->SetLineWidth(2);
  


  //perc_residALL->Draw();
  residALL->Draw("E0");
  //c1->Update();
  TLine *l = new TLine(Min, 0., xAxisMax, 0.);
  l->SetLineStyle(8);
  l->Draw();

  //c1->cd(1);
  p1->cd();
  
  dataALL->SetTitle("");
  dataALL->SetMarkerColor(color?kBlue:kBlack);
  dataALL->SetMarkerStyle(24);
  dataALL->SetMarkerSize(0.75);
  dataALL->SetLineColor(color?kBlue:kBlack); 
  //dataALL->SetFillColor(kBlue);
  dataALL->SetLineWidth(3);
  dataALL->GetXaxis()->SetRangeUser(0., xAxisMax);
  simALL->SetTitle("");
  simALL->SetLineWidth(2);
  simALL->SetMarkerColor(kBlack);
  simALL->SetLineColor(kBlack);
  simALL->SetMarkerSize(0);
  simALL->SetMarkerStyle(24);
  simALL->SetMaximum(dataALL->GetMaximum()*1.2);
  simALL->GetYaxis()->SetTitle("event rate (mHz/keV)");
  simALL->GetXaxis()->SetRangeUser(0., xAxisMax);
  simALL->Draw("HIST C");
  dataALL->Draw("SAMEE0");

  BG_dataALL->SetMarkerColor(1);
  BG_dataALL->SetMarkerStyle(20);
  BG_dataALL->SetMarkerSize(0.6);
  BG_dataALL->SetLineColor(1);
  BG_dataALL->SetLineWidth(3);
  BG_dataALL->Draw("SAMEE0");

  TLegend *leg = new TLegend(0.58,0.65,0.85,0.85);
  leg->AddEntry(dataALL,"Data","p");
  leg->AddEntry(simALL,"Monte Carlo","l");
  leg->AddEntry(BG_dataALL,"Background","p");
  leg->Draw("SAME");

  c1->Print(TString::Format("Spectra_Octets_%i-%i%s.pdf",octetStart,octetEnd,(color?"_color":"")));
  
   
}
  
