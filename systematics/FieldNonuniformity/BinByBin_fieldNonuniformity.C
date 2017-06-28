
// Function to return beta when given the kinetic energy of an electron
Double_t returnBeta(Double_t En) { 
  Double_t me = 510.998928; //rest mass energy of electron in keV
  return sqrt(En*En+2.*En*me)/(En+me);
};


void BinByBin_fieldNonuniformity(Double_t emin, Double_t emax) {

  //std::vector <Double_t> A;

  TString field = "good"; // "good", "bad", "flat"

  //TFile *fout = new TFile(TString::Format("%sField_BinByBin.root",
  //field.Data()),"RECREATE");
  //TH1D *h = new TH1D("h",TString::Format("Field Non-uniformity A_{raw} for %0.0f-%0.0f keV",emin,emax),250, -0.05,0.05);
  //TH1D *h2 = new TH1D("h2",TString::Format("Initial A_{raw} for %0.0f-%0.0f keV",emin,emax),250, -0.05,0.05);
  //TCanvas *c1 = new TCanvas("c1","c1");

  //TH1D* h[100];

  /*Double_t rE[100] = {0.};
  Double_t rW[100] = {0.};
  Double_t rEtrue[100] = {0.};
  Double_t rWtrue[100] = {0.};
  */

  Int_t binMin = emin/10 + 1;
  Int_t binMax = emax/10 ;

  std::cout << binMin << " " << binMax << endl;
  Double_t rE_polE = 0.;
  Double_t rW_polE = 0.;
  Double_t rE_polW = 0.;
  Double_t rW_polW = 0.;
 
  /*for (int i=0; i<100; ++i) {
    binMin[i] = i*10.;
    binMax[i] = binMin[i]+10.;
    h = new TH1D(TString::Format("h%0.0f-%0.0f",binMin,binMax),
		 TString::Format("A %0.0f-%0.0f",binMin,binMax),100,-0.05,0.05);
		 }*/
		 

  

  TString basepath = TString::Format("/extern/mabrow05/ucna/geant4work/output/bigSim/%sField_2011-2012/",field.Data());

  TFile *f;
  

  //Start with East polarization
  for (int i=0; i<5; ++i) {
    std::cout << "In Folder " << i << std::endl;
    for (int j=0; j<1000; ++j) {
      std::cout << j << " rW_polE = " << rW_polE <<  std::endl;
      f = new TFile(TString::Format("%s/%s%i/analyzed_%i.root",
					   basepath.Data(),
					   (field==TString("good")?"G":
					    (field==TString("bad")?"B":"F")),i,j),"READ");
      TTree *t = (TTree*)f->Get("anaTree");
      TH1D *his = new TH1D("his","his",100,0.,1000.);
      //t->Draw("primKE>>his",TString::Format("(EdepSD[5]>0. && hitCountSD[15]==0 && EdepQE>0.)*(1+0.12*sqrt(primKE*primKE+2.*primKE*510.998928)/(primKE+510.998928)*TMath::Cos(primTheta))"),"goff");
      t->Draw("primKE>>his",TString::Format("((EdepQE>0. && EdepQW==0. && MWPCEnergyE>0. && MWPCEnergyW==0.))*(1+0.12*sqrt(primKE*primKE+2.*primKE*510.998928)/(primKE+510.998928)*TMath::Cos(primTheta))"),"goff");
      
      //cout << his->Integral(binMin,binMax) << endl;
      rE_polE += his->Integral(binMin,binMax);

      //t->Draw("primKE>>his",TString::Format("(EdepSD[15]>0. && hitCountSD[5]==0 && EdepQW>0.)*(1+0.12*sqrt(primKE*primKE+2.*primKE*510.998928)/(primKE+510.998928)*TMath::Cos(primTheta))"),"goff");
      t->Draw("primKE>>his",TString::Format("((EdepQW>0. && EdepQE==0. && MWPCEnergyW>0. && MWPCEnergyE==0.))*(1+0.12*sqrt(primKE*primKE+2.*primKE*510.998928)/(primKE+510.998928)*TMath::Cos(primTheta))"),"goff");

      rW_polE += his->Integral(binMin,binMax);
      
      delete his;
      delete f;
      
    }
  }

  //Now West polarization
  for (int i=5; i<10; ++i) {
    std::cout << "In Folder " << i << std::endl;
    for (int j=0; j<1000; ++j) {
      std::cout << j << " rW_polW = " << rW_polW <<  std::endl;
      f = new TFile(TString::Format("%s/%s%i/analyzed_%i.root",
					   basepath.Data(),
					   (field==TString("good")?"G":
					    (field==TString("bad")?"B":"F")),i,j),"READ");
      TTree *t = (TTree*)f->Get("anaTree");
      TH1D *his = new TH1D("his","his",100,0.,1000.);

      //t->Draw("primKE>>his",TString::Format("(EdepSD[5]>0. && hitCountSD[15]==0 && EdepQE>0.)*(1-0.12*sqrt(primKE*primKE+2.*primKE*510.998928)/(primKE+510.998928)*TMath::Cos(primTheta))"),"goff");
      t->Draw("primKE>>his",TString::Format("((EdepQE>0. && EdepQW==0. && MWPCEnergyE>0. && MWPCEnergyW==0.))*(1-0.12*sqrt(primKE*primKE+2.*primKE*510.998928)/(primKE+510.998928)*TMath::Cos(primTheta))"),"goff");
      //t->Draw("primKE>>his",TString::Format("((EdepQE>0. && EdepQW==0.) || (EdepQE>0. && EdepQW>0. && timeE<timeW))*(1-0.12*sqrt(primKE*primKE+2.*primKE*510.998928)/(primKE+510.998928)*TMath::Cos(primTheta))"),"goff");
      rE_polW += his->Integral(binMin,binMax);

      //t->Draw("primKE>>his",TString::Format("(EdepSD[15]>0. && hitCountSD[5]==0 && EdepQW>0.)*(1-0.12*sqrt(primKE*primKE+2.*primKE*510.998928)/(primKE+510.998928)*TMath::Cos(primTheta))"),"goff");
      //t->Draw("primKE>>his",TString::Format("((EdepQW>0. && EdepQE==0.) || (EdepQE>0. && EdepQW>0. && timeW<timeE))*(1-0.12*sqrt(primKE*primKE+2.*primKE*510.998928)/(primKE+510.998928)*TMath::Cos(primTheta))"),"goff");
      t->Draw("primKE>>his",TString::Format("((EdepQW>0. && EdepQE==0. && MWPCEnergyW>0. && MWPCEnergyE==0.))*(1-0.12*sqrt(primKE*primKE+2.*primKE*510.998928)/(primKE+510.998928)*TMath::Cos(primTheta))"),"goff");
      rW_polW += his->Integral(binMin,binMax);
      
      delete his;
      delete f;
      
    }
  }
  
  Double_t rE_polE_err = TMath::Sqrt(rE_polE);
  Double_t rW_polE_err = TMath::Sqrt(rW_polE);
  Double_t rE_polW_err = TMath::Sqrt(rE_polW);
  Double_t rW_polW_err = TMath::Sqrt(rW_polW);

  Double_t R = (rE_polE*rW_polW) / (rE_polW*rW_polE);
  Double_t deltaR = TMath::Sqrt( TMath::Power( rE_polE_err*rW_polW/(rE_polW*rW_polE), 2 ) + 
				 TMath::Power( rW_polW_err*rE_polE/(rE_polW*rW_polE), 2 ) +
				 TMath::Power( rE_polW_err*(rE_polE*rW_polW) / (rE_polW*rE_polW*rW_polE), 2 ) +
				 TMath::Power( rW_polE_err*(rE_polE*rW_polW) / (rE_polW*rW_polE*rW_polE), 2 ) );
  

  Double_t asymm = (1.-TMath::Sqrt(R))/(1.+TMath::Sqrt(R));
  Double_t asymmErr = (deltaR)/(TMath::Sqrt(R)*TMath::Power((TMath::Sqrt(R)+1.),2));
 
  std::cout << " A = " << asymm << " +/- " << asymmErr << std::endl;

  ofstream ofile(TString::Format("asymm_%sField_%0.0f-%0.0f.txt",field.Data(),emin,emax).Data());
  ofile << " A = " << asymm << " +/- " << asymmErr << std::endl;
  ofile.close();
  


}
