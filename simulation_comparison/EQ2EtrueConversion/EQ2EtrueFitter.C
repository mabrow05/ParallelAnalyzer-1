#include <vector>
#include <iomanip>

void EQ2EtrueFitter(TString geom) {

  ofstream params(TString::Format("%s_EQ2EtrueFitParams.dat",geom.Data()).Data()); //Output file for fit parameters
  TString fileName  = TString::Format("HistMeans_%s.dat",geom.Data());
  ifstream infile(fileName.Data());

  int type0lowOffset = 10;
  int type1lowOffset = 13;
  int type23lowOffset = 13;
  
  int type0highOffset = 2;
  int type1highOffset = 4;
  int type23highOffset = 4;

  vector <string> nameRow(23);

  vector <double> EtrueMin;
  vector <double> EtrueMax;
  vector <double> EtrueMid;
  vector <double> EQtype0E;
  vector <double> EQtype1E;
  vector <double> EQtype2E;
  vector <double> EQtype3E;
  vector <double> EQtype23E;
  vector <double> EQtype0W;
  vector <double> EQtype1W;
  vector <double> EQtype2W;
  vector <double> EQtype3W;
  vector <double> EQtype23W;
  // Errors
  vector <double> err_EQtype0E;
  vector <double> err_EQtype1E;
  vector <double> err_EQtype2E;
  vector <double> err_EQtype3E;
  vector <double> err_EQtype23E;
  vector <double> err_EQtype0W;
  vector <double> err_EQtype1W;
  vector <double> err_EQtype2W;
  vector <double> err_EQtype3W;
  vector <double> err_EQtype23W;

  double nEvents;
  int nHist,Emin,Emax;
  double type0E,type1E,type2E,type3E,type23E;
  double type0W,type1W,type2W,type3W,type23W;
  
  double err_type0E,err_type1E,err_type2E,err_type3E,err_type23E;
  double err_type0W,err_type1W,err_type2W,err_type3W,err_type23W;
  
  for (int i=0; i<23; i++) {
    infile >> nameRow[i];
    //cout << nameRow[i] << " " ;
  }
  //cout << endl;
  
  while (infile >> nHist >> Emin >> Emax >> type0E >> err_type0E >> type0W >> err_type0W >>
	 type1E >> err_type1E >> type1W >> err_type1W >> type2E >> err_type2E >> type2W >> err_type2W >>
	 type3E >> err_type3E >> type3W >> err_type3W >> type23E >> err_type23E >> type23W >> err_type23W) {
    //if (Emin<=50.) continue;
    cout << Emin << " " << type0E << endl;
    EtrueMin.push_back(Emin);
    EtrueMax.push_back(Emax);
    EtrueMid.push_back((Emin+Emax)/2.);
    EQtype0E.push_back(type0E);
    EQtype1E.push_back(type1E);
    EQtype2E.push_back(type2E);
    EQtype3E.push_back(type3E);
    EQtype23E.push_back(type23E);
    EQtype0W.push_back(type0W);
    EQtype1W.push_back(type1W);
    EQtype2W.push_back(type2W);
    EQtype3W.push_back(type3W);
    EQtype23W.push_back(type23W);

    err_EQtype0E.push_back(err_type0E);
    err_EQtype1E.push_back(err_type1E);
    err_EQtype2E.push_back(err_type2E);
    err_EQtype3E.push_back(err_type3E);
    err_EQtype23E.push_back(err_type23E);
    err_EQtype0W.push_back(err_type0W);
    err_EQtype1W.push_back(err_type1W);
    err_EQtype2W.push_back(err_type2W);
    err_EQtype3W.push_back(err_type3W);
    err_EQtype23W.push_back(err_type23W);

    if (infile.eof()) break;
  }
  infile.close();

  
  
  int numDataPoints = EtrueMin.size();
  //cout << numDataPoints << " " << EQtype3E.size() << " " <<EtrueMid.size() <<  endl;

  TF1 *func0 = new TF1("func0","[0]+[1]*x+[2]/x+[3]/(x*x)",40.,1000.);
  //TF1 *func0 = new TF1("func0","[0]+[1]*x+[2]/(x+[3])+[4]/((x+[5])*(x+[5]))",10.,1000.);
  func0->SetParameters(100.,1.,0.,0.);
  //func0->FixParameter(0,0.);
  func0->SetLineColor(kBlue);

  TF1 *func1 = new TF1("func1","[0]+[1]*x+[2]/x+[3]/(x*x)",40.,1000.);
  //TF1 *func1 = new TF1("func0","[0]+[1]*x+[2]/(x+[3])+[4]/((x+[5])*(x+[5]))",10.,1000.);
  func1->SetParameters(100.,1.,0.,0.);
  //func1->SetParameters(100.,1.,1.,0.,1.,0.);
  //func1->FixParameter(0,0.);
  func1->SetLineColor(kRed);

  TF1 *func23 = new TF1("func23","[0]+[1]*x+[2]/x+[3]/(x*x)",40.,1000.);
  //TF1 *func23 = new TF1("func0","[0]+[1]*x+[2]/(x+[3])+[4]/((x+[5])*(x+[5]))",10.,1000.);
  func23->SetParameters(100.,1.,0.,0.);
  //func23->SetParameters(100.,1.,1.,10.,1.,10.);
  //func23->SetParLimits(0,10., 200.);
  //func23->SetParLimits(3,0., 10000.);
  //func23->SetParLimits(5,0., 10000.);
  //func23->FixParameter(0,0.);
  func23->SetLineColor(kGreen);


  TGraphErrors *t0E = new TGraphErrors(numDataPoints-(type0lowOffset+type0highOffset),&EQtype0E[type0lowOffset-1],&EtrueMid[type0lowOffset-1],&err_EQtype0E[type0lowOffset-1],0); 
  t0E->SetMarkerColor(kBlue);
  t0E->SetMarkerStyle(21);
  //t0E->SetMarkerSize(.8);
  TGraphErrors *t1E = new TGraphErrors(numDataPoints-(type1lowOffset+type1highOffset),&EQtype1E[type1lowOffset-1],&EtrueMid[type1lowOffset-1],&err_EQtype1E[type1lowOffset-1],0);
  t1E->SetMarkerColor(kRed);
  t1E->SetMarkerStyle(20);
  TGraphErrors *t2E = new TGraphErrors(numDataPoints-11,&EQtype2E[8],&EtrueMid[8],&err_EQtype2E[8],0);
  TGraphErrors *t3E = new TGraphErrors(numDataPoints-11,&EQtype3E[8],&EtrueMid[8],&err_EQtype3E[8],0);
  TGraphErrors *t23E = new TGraphErrors(numDataPoints-(type23lowOffset+type23highOffset),&EQtype23E[type23lowOffset-1],&EtrueMid[type23lowOffset-1],&err_EQtype23E[type23lowOffset-1],0);
  t23E->SetMarkerColor(kGreen);
  t23E->SetMarkerStyle(22);

  t0E->Fit(func0,"R");
  t1E->Fit(func1,"R");
  t23E->Fit(func23,"R");

  params << std::fixed << std::setprecision(17);
  
  params << "Type0E: " << func0->GetParameter(0) << " " << func0->GetParameter(1) << " " << func0->GetParameter(2) << " " << func0->GetParameter(3) 
	 << " " << func0->GetParameter(4) << " " << func0->GetParameter(5) << endl;
  params << "Type1E: " << func1->GetParameter(0) << " " << func1->GetParameter(1) << " " << func1->GetParameter(2) << " " << func1->GetParameter(3) 
	 << " " << func1->GetParameter(4) << " " << func1->GetParameter(5) << endl;
  params << "Type23E: " << func23->GetParameter(0) << " " << func23->GetParameter(1) << " " << func23->GetParameter(2) << " " << func23->GetParameter(3) 
	 << " " << func23->GetParameter(4) << " " << func23->GetParameter(5) << endl;

  TGraphErrors *t0W = new TGraphErrors(numDataPoints-(type0lowOffset+type0highOffset),&EQtype0W[type0lowOffset-1],&EtrueMid[type0lowOffset-1],&err_EQtype0W[type0lowOffset-1],0); 
  t0W->SetMarkerColor(kBlue);
  t0W->SetMarkerStyle(21);
  TGraphErrors *t1W = new TGraphErrors(numDataPoints-(type1lowOffset+type1highOffset),&EQtype1W[type1lowOffset-1],&EtrueMid[type1lowOffset-1],&err_EQtype1W[type1lowOffset-1],0);
  t1W->SetMarkerColor(kRed);
  t1W->SetMarkerStyle(20);
  TGraphErrors *t2W = new TGraphErrors(numDataPoints-11,&EQtype2W[8],&EtrueMid[8],&err_EQtype2W[8],0);
  TGraphErrors *t3W = new TGraphErrors(numDataPoints-11,&EQtype3W[8],&EtrueMid[8],&err_EQtype3W[8],0);
  TGraphErrors *t23W = new TGraphErrors(numDataPoints-(type23lowOffset+type23highOffset),&EQtype23W[type23lowOffset-1],&EtrueMid[type23lowOffset-1],&err_EQtype23W[type23lowOffset-1],0);
  t23W->SetMarkerColor(kGreen);
  t23W->SetMarkerStyle(22);

  t0W->Fit(func0,"R");
  t1W->Fit(func1,"R");
  t23W->Fit(func23,"R");

   params << "Type0W: " << func0->GetParameter(0) << " " << func0->GetParameter(1) << " " << func0->GetParameter(2) << " " << func0->GetParameter(3) 
	 << " " << func0->GetParameter(4) << " " << func0->GetParameter(5) << endl;
  params << "Type1W: " << func1->GetParameter(0) << " " << func1->GetParameter(1) << " " << func1->GetParameter(2) << " " << func1->GetParameter(3) 
	 << " " << func1->GetParameter(4) << " " << func1->GetParameter(5) << endl;
  params << "Type23W: " << func23->GetParameter(0) << " " << func23->GetParameter(1) << " " << func23->GetParameter(2) << " " << func23->GetParameter(3) 
	 << " " << func23->GetParameter(4) << " " << func23->GetParameter(5) << endl;
 
  TMultiGraph *east = new TMultiGraph();
  east->Add(t0E);
  east->Add(t1E);
  east->Add(t23E);
  east->SetTitle("East Detector");

  TMultiGraph *west = new TMultiGraph();
  west->Add(t0W);
  west->Add(t1W);
  west->Add(t23W);
  west->SetTitle("West Detector");
  
  TCanvas *c1 = new TCanvas("c1", "c1", 1200, 500);
  c1->Divide(2,1);
  c1->cd(1);
  east->Draw("AP");
  east->SetMinimum(0.);
  east->SetMaximum(800.);
  east->GetXaxis()->SetLimits(0.,800.);
  east->GetXaxis()->SetTitle("E_{Q} (keV)");
  east->GetYaxis()->SetTitle("E_{true} (keV)");
  east->GetXaxis()->CenterTitle();
  east->GetYaxis()->CenterTitle();
  east->GetXaxis()->SetTitleOffset(1.4);
  east->GetYaxis()->SetTitleOffset(1.4);
  east->Draw("AP");

  TLegend *leg0 = new TLegend(0.55,0.25,0.85,0.375);
  leg0->AddEntry(t0E,"  Type 0","p");
  leg0->AddEntry(t1E,"  Type 1","p");
  leg0->AddEntry(t23E,"  Type 2/3","p");
  leg0->Draw();

  c1->cd(2);
  west->Draw("AP");
  west->SetMinimum(0.);
  west->SetMaximum(800.);
  west->GetXaxis()->SetLimits(0.,800.);
  west->GetXaxis()->SetTitle("E_{Q} (keV)");
  west->GetYaxis()->SetTitle("E_{true} (keV)");
  west->GetXaxis()->CenterTitle();
  west->GetYaxis()->CenterTitle();
  west->GetXaxis()->SetTitleOffset(1.4);
  west->GetYaxis()->SetTitleOffset(1.4);
  west->Draw("AP");

  params.close();
  
  

}
