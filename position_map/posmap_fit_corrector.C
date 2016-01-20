

{
  int iRunPeriod = 10;
  int pmt = 7; 
  float xpos =7.5;
  float ypos = 0.;

  //open file
  Char_t temp[500];
  sprintf(temp,"position_map_%i_RC_123.root",iRunPeriod);
  TFile *f = new TFile(temp,"READ");
  
  //Load histogram
  TCanvas *c1 = new TCanvas("c1");
  sprintf(temp,"%s%i_%i_%i",pmt<4?"e":"w",pmt<4?pmt:(pmt-4),(int)xpos,(int)ypos);
  cout << temp << endl;
  TH1F *h = new TH1F("h","h",800, 0., 4000.);
  h = (TH1F*)(f->Get(temp));

  h->Draw();

  TF1 *f1 = new TF1("f1","[0]*exp(-(x-[1])*(x-[1])/(2.*[2]*[2]))",100., 400.);
  f1->SetParameter(0,100.);
  f1->SetParameter(1,300.);
  f1->SetParameter(2,100.);
  f1->SetParLimits(0,0.,5000.);
  f1->SetParLimits(1,0.,500.);	
  f1->SetParLimits(2,0.,0.);
  h->Fit("f1", "LR");
}
