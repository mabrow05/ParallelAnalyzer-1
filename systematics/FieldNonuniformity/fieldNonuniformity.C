

void fieldNonuniformity(Double_t emin, Double_t emax) {

  //std::vector <Double_t> A;

  TString field = "good"; // "good", "bad", "flat"

  TFile *fout = new TFile(TString::Format("%sField_%0.0f-%0.0f.root",
					  field.Data(),emin,emax),"RECREATE");

  //TCanvas *c1 = new TCanvas("c1","c1");

  TH1D *h = new TH1D("h",TString::Format("Field Non-uniformity A_{raw} for %0.0f-%0.0f keV",emin,emax),250, -0.05,0.05);
  TH1D *h2 = new TH1D("h2",TString::Format("Initial A_{raw} for %0.0f-%0.0f keV",emin,emax),250, -0.05,0.05);

  TString basepath = TString::Format("/extern/mabrow05/ucna/geant4work/output/bigSim/%sField_2011-2012/",field.Data());

  for (int i=0; i<10; ++i) {
    std::cout << "In Folder " << i << std::endl;
    for (int j=0; j<1000; ++j) {
      
      TFile *f = new TFile(TString::Format("%s/%s%i/analyzed_%i.root",
					   basepath.Data(),
					   (field==TString("good")?"G":
					    (field==TString("bad")?"B":"F")),i,j),"READ");
      TTree *t = (TTree*)f->Get("anaTree");

      Long_t rE = t->GetEntries(TString::Format("EdepSD[5]>0. && hitCountSD[15]==0 && EdepE>0. && primKE>%f && primKE<%f",emin, emax));
      Long_t rW = t->GetEntries(TString::Format("EdepSD[15]>0. && hitCountSD[5]==0 && EdepW>0. && primKE>%f && primKE<%f",emin, emax));
      Double_t asymm = (Double_t)(rE-rW)/(Double_t)(rE+rW);

      Long_t rEtrue = t->GetEntries(TString::Format("((EdepSD[5]>0. && hitCountSD[15]==0 && EdepE>0.) || (EdepSD[15]>0. && hitCountSD[5]==0 && EdepW>0.)) && primKE>%f && primKE<%f && primTheta>TMath::Pi()/2.",emin, emax));
      Long_t rWtrue = t->GetEntries(TString::Format("((EdepSD[5]>0. && hitCountSD[15]==0 && EdepE>0.) || (EdepSD[15]>0. && hitCountSD[5]==0 && EdepW>0.)) && primKE>%f && primKE<%f && primTheta<TMath::Pi()/2.",emin, emax));
      Double_t asymm_true = (Double_t)(rEtrue-rWtrue)/(Double_t)(rEtrue+rWtrue);

      h->Fill(asymm);
      h2->Fill(asymm_true);
      //if (j%1==0) std::cout << j << " deltaA = " << asymm-asymm_true << std::endl;
      if (j%1==0) std::cout << j << " deltaA = " << asymm << std::endl;
      
      delete f;

    }

  }

  
  //h->Draw();
  h->Fit("gaus");
  h2->Fit("gaus");

  fout->Write();
  delete fout;

}
