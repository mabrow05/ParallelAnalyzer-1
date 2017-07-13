{


  int run = 18222;//21769;//
  
  TFile *f = new TFile(TString::Format("%s/beta_highStatistics/revCalSim_%i_Beta.root",getenv("REVCALSIM"),run));

  TTree *t = (TTree*)f->Get("revCalSim");

  TH1D *triggE = new TH1D("triggE","East Trigger",80,0.,800.);  
  TH1D *AllE = new TH1D("AllE","East Trigger",80,0.,800.);
  TH1D *triggW = new TH1D("triggW","West Trigger",80,0.,800.);  
  TH1D *AllW = new TH1D("AllW","West Trigger",80,0.,800.);

  t->Draw("EdepQE>>triggE","(MWPCPosE[0]*MWPCPosE[0]+MWPCPosE[1]*MWPCPosE[1])<25./0.6 && MWPCEnergyE>0. && EdepQE>0. && (side==0 || type==1)","GOFF"); 
  t->Draw("EdepQE>>AllE","(MWPCPosE[0]*MWPCPosE[0]+MWPCPosE[1]*MWPCPosE[1])<25./0.6 && MWPCEnergyE>0. && EdepQE>0.","GOFF"); 

  t->Draw("EdepQW>>triggW","(MWPCPosW[0]*MWPCPosW[0]+MWPCPosW[1]*MWPCPosW[1])<25./0.6 && MWPCEnergyW>0. && EdepQW>0. && (side==1 || type==1 )","GOFF"); 
  t->Draw("EdepQW>>AllW","(MWPCPosW[0]*MWPCPosW[0]+MWPCPosW[1]*MWPCPosW[1])<25./0.6 && MWPCEnergyW>0. && EdepQW>0.","GOFF"); 

  TCanvas c1;
  TGraphAsymmErrors efficE;
  efficE.Divide(triggE,AllE);
  efficE.Draw("AP");

  TF1 *fE = new TF1("fE","0.5+0.5*TMath::Erf((x-[0])/[1])",0.,600.);
  fE->SetParameters(-20.,50.);
  efficE->Fit(fE,"R");

  TCanvas c2;
  TGraphAsymmErrors efficW;
  efficW.Divide(triggW,AllW);
  efficW.Draw("AP");

  TF1 *fW = new TF1("fW","0.5+0.5*TMath::Erf((x-[0])/[1])",0.,600.);
  fW->SetParameters(-20.,50.);
  efficW->Fit(fW,"R");
}
