/*//////////////////////////////////////////////////////////////////////////////

This script reads in the events from the processed beta files and 
creates energy spectra from which b can be extracted

//////////////////////////////////////////////////////////////////////////////*/

//#include <string>

void spectraMaker(std::string file)
{
  
  std::string filePath = "analyzed_files/"+file;

  TFile *fout = new TFile("outputHistsTest.root","RECREATE");
  TH1D *hEreconALL = new TH1D("hEreconALL","Reconstructed Energy Spectra for All events",120,0., 1200.);
  TH1D *hErecon0 = new TH1D("hErecon0","Reconstructed Energy Spectra for All events",120,0., 1200.);
  TH1D *hErecon1 = new TH1D("hErecon1","Reconstructed Energy Spectra for All events",120,0., 1200.);
  TH1D *hErecon23 = new TH1D("hErecon23","Reconstructed Energy Spectra for All events",120,0., 1200.);

  hEreconALL->GetXaxis()->SetTitle("E_{recon} (keV)");
  hErecon0->GetXaxis()->SetTitle("E_{recon} (keV)");
  hErecon1->GetXaxis()->SetTitle("E_{recon} (keV)");
  hErecon23->GetXaxis()->SetTitle("E_{recon} (keV)");
  
  TFile *fin = new TFile(filePath.c_str(),"READ");
  TTree *tin = (TTree*)(fin->Get("SimAnalyzed"));

  //variables to be read in
  Int_t PID, side, type;
  Double_t Erecon;
  
  tin->SetBranchAddress("PID",&PID);
  tin->SetBranchAddress("side",&side);
  tin->SetBranchAddress("type",&type);
  tin->SetBranchAddress("Erecon",&Erecon);

  UInt_t nevents = tin->GetEntriesFast();

  for (UInt_t evt=0; evt<nevents; evt++) {
    tin->GetEvent(evt);

    if (PID!=1) continue;

    hEreconALL->Fill(Erecon);
    if (type==0) hErecon0->Fill(Erecon);
    else if (type==1) hErecon1->Fill(Erecon);
    else if (type==2 || type==3) hErecon23->Fill(Erecon);  
  }

  fin->Close();

  fout->Write();
  fout->Close();
}
