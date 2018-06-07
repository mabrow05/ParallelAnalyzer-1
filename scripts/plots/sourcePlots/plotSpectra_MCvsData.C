






void plotSpectra_MCvsData(int run) {

  TH1::AddDirectory(kFALSE);

  gStyle->SetPadTopMargin(0.1);
  gStyle->SetOptStat(0000);
  gStyle->SetPadRightMargin(0.05);
  gStyle->SetPadBottomMargin(0.12);
  gStyle->SetLabelSize(0.05,"xy");
  gStyle->SetTitleSize(0.06,"t");

  //All histograms to be drawn
  TH1D *hisEvisTot_data[3][2];
  TH1D *hisEvis_data[3][8];
  TH1D *hisEreconTot_data[3][2];
  
  TH1D *hisEvisTot_sim[3][2];
  TH1D *hisEvis_sim[3][8];
  TH1D *hisEreconTot_sim[3][2];

  //open data file and sim file
  TFile *fData = new TFile(TString::Format("%s/source_peaks_%i_Erecon.root",
					   getenv("SOURCE_PEAKS"),run),
			   "READ");


  TFile *fSim = new TFile(TString::Format("%s/source_peaks/source_peaks_%i_Erecon.root",
					   getenv("REVCALSIM"),run),
			   "READ");

  if ( !(fData->IsOpen() && fSim->IsOpen()) ) exit(0);

  for ( int i=0; i<3; ++i ) {

    for ( int p=0; p<4; ++p ) {
      
      hisEvis_data[i][p] = (TH1D*)fData->Get(TString::Format("hisEvis%i_E%i",
							     i+1,p));
      hisEvis_data[i][p+4] = (TH1D*)fData->Get(TString::Format("hisEvis%i_W%i",
							     i+1,p));

      hisEvis_sim[i][p] = (TH1D*)fSim->Get(TString::Format("hisEvis%i_E%i",
							   i+1,p));
      hisEvis_sim[i][p+4] = (TH1D*)fSim->Get(TString::Format("hisEvis%i_W%i",
							     i+1,p));
    }
    
    hisEreconTot_data[i][0] = (TH1D*)fData->Get(TString::Format("hisErecon%iE",
								i+1));
    hisEreconTot_data[i][1] = (TH1D*)fData->Get(TString::Format("hisErecon%iW",
								i+1));

    hisEreconTot_sim[i][0] = (TH1D*)fSim->Get(TString::Format("hisErecon%iE",
								i+1));
    hisEreconTot_sim[i][1] = (TH1D*)fSim->Get(TString::Format("hisErecon%iW",
								i+1));
    
    hisEvisTot_data[i][0] = (TH1D*)fData->Get(TString::Format("hisEvis%iE",
								i+1));
    hisEvisTot_data[i][1] = (TH1D*)fData->Get(TString::Format("hisEvis%iW",
								i+1));

    hisEvisTot_sim[i][0] = (TH1D*)fSim->Get(TString::Format("hisEvis%iE",
								i+1));
    hisEvisTot_sim[i][1] = (TH1D*)fSim->Get(TString::Format("hisEvis%iW",
								i+1));
  }

  delete fData;
  delete fSim;

  TFile *ofile = new TFile(TString::Format("run_%i.root",run),"RECREATE");
  //hisEreconTot_data[0][0]->Write();
  for (int jj=0;jj<2;++jj) {
    for (int j=0;j<3;++j) {
      hisEreconTot_data[j][jj]->Write();
      hisEreconTot_sim[j][jj]->SetName(TString::Format("SIM_hisErecon%i%s",j+1,jj==0?"E":"W"));
      hisEreconTot_sim[j][jj]->SetLineColor(kRed);
      hisEreconTot_sim[j][jj]->Write();
    }
  }
  delete ofile;


  int nBins = hisEreconTot_data[0][0]->GetNbinsX();
  int bin0 = hisEvis_data[0][0]->GetXaxis()->FindBin(10.);
  int binLast = hisEvis_data[0][0]->GetNbinsX();

  TCanvas *cPMT[3];
  cPMT[0] = new TCanvas("cPMT1","cPMT1",1200, 800);
  cPMT[0]->Divide(4,2);
  cPMT[1] = new TCanvas("cPMT2","cPMT2",1200, 800);
  cPMT[1]->Divide(4,2);
  cPMT[2] = new TCanvas("cPMT3","cPMT3",1200, 800);
  cPMT[2]->Divide(4,2);

  TCanvas *cErecon = new TCanvas("cErecon","cErecon",1200, 800);
  cErecon->Divide(3,2);

  TCanvas *cEvis = new TCanvas("cEvis","cEvis",1200, 800);
  cEvis->Divide(3,2);

  double normFactor = 0.;

  for ( int i=0; i<3; ++i ) {

    for ( int p=0; p<4; ++p ) {
      
      /////////////////////////////////////////////////
      cPMT[i]->cd(p+1);
      

      if ( hisEvis_data[i][p]->Integral(bin0,binLast) > 0 && hisEvis_sim[i][p]->Integral(bin0,binLast) > 0 ) {
	hisEvis_data[i][p]->Draw("hist"); 
	hisEvis_data[i][p]->SetTitle(TString::Format("PMT East %i",p)); 
	
	hisEvis_data[i][p]->GetXaxis()->SetTitle("E_{vis} [keV]"); 
	hisEvis_data[i][p]->GetXaxis()->SetLabelSize(0.05);
	hisEvis_data[i][p]->GetYaxis()->SetLabelSize(0.05);
	hisEvis_data[i][p]->GetXaxis()->SetTitleSize(0.06); 
 
	hisEvis_data[i][p]->GetXaxis()->SetRangeUser(10.,1500.);
	hisEvis_sim[i][p]->GetXaxis()->SetRangeUser(10.,1500.);
	hisEvis_sim[i][p]->SetLineColor(kRed);
	hisEvis_sim[i][p]->SetLineStyle(2);
	
	normFactor = ( hisEvis_data[i][p]->GetMaximum() / 
		       hisEvis_sim[i][p]->GetMaximum() );
	
	hisEvis_data[i][p]->GetXaxis()->SetRangeUser(-20.,1300.);
	hisEvis_sim[i][p]->GetXaxis()->SetRangeUser(-20.,1300.);
	
	hisEvis_sim[i][p]->Scale(normFactor);
	hisEvis_sim[i][p]->Draw("SAME HIST");
      }
      /////////////////////////////////////////////////
      cPMT[i]->cd(p+5);
      
      if ( hisEvis_data[i][p+4]->Integral(bin0,binLast) > 0  && hisEvis_sim[i][p+4]->Integral(bin0,binLast) > 0 ) {
	hisEvis_data[i][p+4]->Draw("hist");
	hisEvis_data[i][p+4]->GetXaxis()->SetLabelSize(0.05);
	hisEvis_data[i][p+4]->GetYaxis()->SetLabelSize(0.05);
	hisEvis_data[i][p+4]->GetXaxis()->SetTitleSize(0.06); 

	hisEvis_data[i][p+4]->SetTitle(TString::Format("PMT West %i",p)); 
	hisEvis_data[i][p+4]->GetXaxis()->SetTitle("E_{vis} [keV]"); 

	if ( p==1 ) cout << hisEvis_sim[i][p+4]->Integral(bin0,binLast) << endl;
	hisEvis_data[i][p+4]->GetXaxis()->SetRangeUser(10.,1500.);
	hisEvis_sim[i][p+4]->GetXaxis()->SetRangeUser(10.,1500.);
	hisEvis_sim[i][p+4]->SetLineColor(kRed);
	hisEvis_sim[i][p+4]->SetLineStyle(2);
	
	normFactor = ( hisEvis_data[i][p+4]->GetMaximum() / 
		       hisEvis_sim[i][p+4]->GetMaximum() );
	
	hisEvis_data[i][p+4]->GetXaxis()->SetRangeUser(-20.,1300.);
	hisEvis_sim[i][p+4]->GetXaxis()->SetRangeUser(-20.,1300.);
	
	hisEvis_sim[i][p+4]->Scale(normFactor);
	hisEvis_sim[i][p+4]->Draw("SAME HIST");
      }
    }
   
    //////////////////////////////////////////////////////////
    cErecon->cd(i+1);

    if ( hisEreconTot_data[i][0]->Integral(bin0,binLast) > 0 && hisEreconTot_sim[i][0]->Integral(bin0,binLast) > 0 ) {
      hisEreconTot_data[i][0]->Draw("hist");
      hisEreconTot_data[i][0]->SetTitle(i==0?"^{137}Ce East":(i==1?"^{113}Sn East":("^{207}Bi East"))); 
      hisEreconTot_data[i][0]->GetXaxis()->SetTitle("E_{recon} [keV]"); 
      hisEreconTot_data[i][0]->GetXaxis()->SetLabelSize(0.04);
      hisEreconTot_data[i][0]->GetYaxis()->SetLabelSize(0.04);
      hisEreconTot_data[i][0]->GetXaxis()->SetTitleSize(0.05); 

      hisEreconTot_data[i][0]->GetXaxis()->SetRange(2,nBins);
      hisEreconTot_sim[i][0]->GetXaxis()->SetRange(2,nBins);
      hisEreconTot_sim[i][0]->SetLineColor(kRed);
      hisEreconTot_sim[i][0]->SetLineStyle(2);
	
      normFactor = ( hisEreconTot_data[i][0]->GetMaximum() / 
		     hisEreconTot_sim[i][0]->GetMaximum() );
      
      hisEreconTot_data[i][0]->GetXaxis()->SetRangeUser(-20.,1300.);
      hisEreconTot_sim[i][0]->GetXaxis()->SetRangeUser(-20.,1300.);
      
      hisEreconTot_sim[i][0]->Scale(normFactor);
      hisEreconTot_sim[i][0]->Draw("SAME HIST");
    }

    //////////////////////////////////////////////////////////
    cErecon->cd(i+4);
    if ( hisEreconTot_data[i][1]->Integral(bin0,binLast) > 0 && hisEreconTot_sim[i][1]->Integral(bin0,binLast) > 0 ) {
      hisEreconTot_data[i][1]->Draw("hist");
      hisEreconTot_data[i][1]->SetTitle(i==0?"^{137}Ce West":(i==1?"^{113}Sn West":("^{207}Bi West"))); 
      hisEreconTot_data[i][1]->GetXaxis()->SetTitle("E_{recon} [keV]"); 
      hisEreconTot_data[i][1]->GetXaxis()->SetLabelSize(0.04);
      hisEreconTot_data[i][1]->GetYaxis()->SetLabelSize(0.04);
      hisEreconTot_data[i][1]->GetXaxis()->SetTitleSize(0.05); 
      
      hisEreconTot_data[i][1]->GetXaxis()->SetRange(2,nBins);
      hisEreconTot_sim[i][1]->GetXaxis()->SetRange(2,nBins);
      hisEreconTot_sim[i][1]->SetLineColor(kRed);
      hisEreconTot_sim[i][1]->SetLineStyle(2);
      
      normFactor = ( hisEreconTot_data[i][1]->GetMaximum() / 
		     hisEreconTot_sim[i][1]->GetMaximum() );
      
      hisEreconTot_data[i][1]->GetXaxis()->SetRangeUser(-20.,1300.);;
      hisEreconTot_sim[i][1]->GetXaxis()->SetRangeUser(-20.,1300.);
      
      hisEreconTot_sim[i][1]->Scale(normFactor);
      hisEreconTot_sim[i][1]->Draw("SAME HIST");
    }
    //////////////////////////////////////////////////////////
    cEvis->cd(i+1);
    if ( hisEvisTot_data[i][0]->Integral(bin0,binLast) > 0 && hisEvisTot_sim[i][0]->Integral(bin0,binLast) > 0 ) {
      hisEvisTot_data[i][0]->Draw("hist");
      
      hisEvisTot_data[i][0]->GetXaxis()->SetRange(2,nBins);
      hisEvisTot_sim[i][0]->GetXaxis()->SetRange(2,nBins);
      hisEvisTot_sim[i][0]->SetLineColor(kRed);
      hisEvisTot_sim[i][0]->SetLineStyle(2);
      
      normFactor = ( hisEvisTot_data[i][0]->GetMaximum() / 
		     hisEvisTot_sim[i][0]->GetMaximum() );
      
      hisEvisTot_data[i][0]->GetXaxis()->SetRangeUser(-20.,1300.);
      hisEvisTot_sim[i][0]->GetXaxis()->SetRangeUser(-20.,1300.);
      
      hisEvisTot_sim[i][0]->Scale(normFactor);
      hisEvisTot_sim[i][0]->Draw("SAME HIST");
    }
    //////////////////////////////////////////////////////////
    cEvis->cd(i+4);
    if ( hisEvisTot_data[i][1]->Integral(bin0,binLast) > 0 && hisEvisTot_sim[i][1]->Integral(bin0,binLast) > 0 ) {
      hisEvisTot_data[i][1]->Draw("hist");
      
      hisEvisTot_data[i][1]->GetXaxis()->SetRange(2,nBins);
      hisEvisTot_sim[i][1]->GetXaxis()->SetRange(2,nBins);
      hisEvisTot_sim[i][1]->SetLineColor(kRed);
      hisEvisTot_sim[i][1]->SetLineStyle(2);

      normFactor = ( hisEvisTot_data[i][1]->GetMaximum() / 
		     hisEvisTot_sim[i][1]->GetMaximum() );
      
      hisEvisTot_data[i][1]->GetXaxis()->SetRangeUser(-20.,1300.);
      hisEvisTot_sim[i][1]->GetXaxis()->SetRangeUser(-20.,1300.);
      
      hisEvisTot_sim[i][1]->Scale(normFactor);
      hisEvisTot_sim[i][1]->Draw("SAME HIST");
    }
  }


  
  TString file = TString::Format("spectraPlots/run_%i.pdf",run);

  cErecon->Print(file+TString("("));
  cPMT[0]->Print(file);
  cPMT[1]->Print(file);
  cPMT[2]->Print(file);
  cEvis->Print(file+TString(")"));

  std::ofstream ofileSIME(TString::Format("SIM_run_%i_binned_east.txt"),run);
  std::ofstream ofileSIMW(TString::Format("SIM_run_%i_binned_west.txt"),run);

  std::ofstream ofileDATAE(TString::Format("DATA_run_%i_binned_east.txt"),run);
  std::ofstream ofileDATAW(TString::Format("DATA_run_%i_binned_west.txt"),run);
  std::ofstream ofileSource1(TString::Format("Source1_run_%i.txt"),run);
  std::ofstream ofileSource2(TString::Format("Source2_run_%i.txt"),run);
  std::ofstream ofileSource3(TString::Format("Source3_run_%i.txt"),run);

  /*for (int j=0; j<3 ; ++j) {
    for (int jj=1; jj<=hisEvisTot_data[j][0]; ++jj) {
    if ( hisEvisTot_data[j][1]->Integral(bin0,binLast) > 0 && hisEvisTot_sim[j][1]->Integral(bin0,binLast) > 0 ) {      
	ofileSource1 << hisEvisTot_data[j][1]->GetBinContent(jj) << "\t";
	ofileSIMW << hisEvisTot_sim[j][1]->GetBinContent(jj) << "\t";
      }
      else { 
	ofileDATAW << 0. << "\t";
	ofileSIMW << 0. << "\t";
      }
      if ( hisEvisTot_data[j][0]->Integral(bin0,binLast) > 0 && hisEvisTot_sim[j][0]->Integral(bin0,binLast) > 0 ) {      
	ofileDATAE << hisEvisTot_data[j][0]->GetBinContent(jj) << "\t";
	ofileSIME << hisEvisTot_sim[j][0]->GetBinContent(jj) << "\t";
      }
      else { 
	ofileDATAE << 0. << "\t";
	ofileSIME << 0. << "\t";
      }
      
    }
  }
  ofileSIME.close();
  ofileSIMW.close();
  ofileDATAE.close();
  ofileDATAW.close();*/
};
