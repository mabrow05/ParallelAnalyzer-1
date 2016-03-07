#include <vector>

void EQ2EtrueFitter() {
  
  ofstream params("2011-2012_EQ2EtrueFitParams.dat"); //Output file for fit parameters
  string fileName  = "HistMeans.dat";
  ifstream infile(fileName.c_str());

  int type0lowOffset = 5;
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
  double nEvents;
  int nHist,Emin,Emax;
  double type0E,type1E,type2E,type3E,type23E;
  double type0W,type1W,type2W,type3W,type23W;
  
  for (int i=0; i<23; i++) {
    infile >> nameRow[i];
    //cout << nameRow[i] << " " ;
  }
  //cout << endl;
  
  while (infile >> nHist >> Emin >> Emax >> type0E >> nEvents >> type0W >> nEvents >>
	 type1E >> nEvents >> type1W >> nEvents >> type2E >> nEvents >> type2W >> nEvents >>
	 type3E >> nEvents >> type3W >> nEvents >> type23E >> nEvents >> type23W >> nEvents) {
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

    if (infile.eof()) break;
  }
  infile.close();

  
  
  int numDataPoints = EtrueMin.size();
  //cout << numDataPoints << " " << EQtype3E.size() << " " <<EtrueMid.size() <<  endl;

  //TF1 *func0 = new TF1("func0","[0]+[1]*x+[2]/x+[3]/(x*x)",20.,1000.);
  TF1 *func0 = new TF1("func0","[0]+[1]*x+[2]/(x+[3])+[4]/((x+[5])*(x+[5]))",0.,1000.);
  func0->SetParameters(100.,1.,1.,0.,1.,0.);
  //func0->FixParameter(0,0.);
  func0->SetLineColor(kBlue);

  //TF1 *func1 = new TF1("func1","[0]+[1]*x+[2]/x+[3]/(x*x)",40.,1000.);
  TF1 *func1 = new TF1("func0","[0]+[1]*x+[2]/(x+[3])+[4]/((x+[5])*(x+[5]))",0.,1000.);
  func1->SetParameters(100.,1.,1.,0.,1.,0.);
  //func1->FixParameter(0,0.);
  func1->SetLineColor(kRed);

  //TF1 *func23 = new TF1("func23","[0]+[1]*x+[2]/x+[3]/(x*x)",40.,1000.);
  TF1 *func23 = new TF1("func0","[0]+[1]*x+[2]/(x+[3])+[4]/((x+[5])*(x+[5]))",0.,1000.);
  func23->SetParameters(100.,1.,1.,0.,1.,0.);
  //func23->FixParameter(0,0.);
  func23->SetLineColor(kGreen);


  TGraph *t0E = new TGraph(numDataPoints-(type0lowOffset+type0highOffset),&EQtype0E[type0lowOffset-1],&EtrueMid[type0lowOffset-1]); 
  t0E->SetMarkerColor(kBlue);
  t0E->SetMarkerStyle(21);
  //t0E->SetMarkerSize(.8);
  TGraph *t1E = new TGraph(numDataPoints-(type1lowOffset+type1highOffset),&EQtype1E[type1lowOffset-1],&EtrueMid[type1lowOffset-1]);
  t1E->SetMarkerColor(kRed);
  t1E->SetMarkerStyle(20);
  TGraph *t2E = new TGraph(numDataPoints-11,&EQtype2E[8],&EtrueMid[8]);
  TGraph *t3E = new TGraph(numDataPoints-11,&EQtype3E[8],&EtrueMid[8]);
  TGraph *t23E = new TGraph(numDataPoints-(type23lowOffset+type23highOffset),&EQtype23E[type23lowOffset-1],&EtrueMid[type23lowOffset-1]);
  t23E->SetMarkerColor(kGreen);
  t23E->SetMarkerStyle(22);

  t0E->Fit(func0,"R");
  t1E->Fit(func1,"R");
  t23E->Fit(func23,"R");
  
  params << "Type0E: " << func0->GetParameter(0) << " " << func0->GetParameter(1) << " " << func0->GetParameter(2) << " " << func0->GetParameter(3) 
	 << " " << func0->GetParameter(4) << " " << func0->GetParameter(5) << endl;
  params << "Type1E: " << func1->GetParameter(0) << " " << func1->GetParameter(1) << " " << func1->GetParameter(2) << " " << func1->GetParameter(3) 
	 << " " << func1->GetParameter(4) << " " << func1->GetParameter(5) << endl;
  params << "Type23E: " << func23->GetParameter(0) << " " << func23->GetParameter(1) << " " << func23->GetParameter(2) << " " << func23->GetParameter(3) 
	 << " " << func23->GetParameter(4) << " " << func23->GetParameter(5) << endl;

  TGraph *t0W = new TGraph(numDataPoints-(type0lowOffset+type0highOffset),&EQtype0W[type0lowOffset-1],&EtrueMid[type0lowOffset-1]); 
  t0W->SetMarkerColor(kBlue);
  t0W->SetMarkerStyle(21);
  TGraph *t1W = new TGraph(numDataPoints-(type1lowOffset+type1highOffset),&EQtype1W[type1lowOffset-1],&EtrueMid[type1lowOffset-1]);
  t1W->SetMarkerColor(kRed);
  t1W->SetMarkerStyle(20);
  TGraph *t2W = new TGraph(numDataPoints-11,&EQtype2W[8],&EtrueMid[8]);
  TGraph *t3W = new TGraph(numDataPoints-11,&EQtype3W[8],&EtrueMid[8]);
  TGraph *t23W = new TGraph(numDataPoints-(type23lowOffset+type23highOffset),&EQtype23W[type23lowOffset-1],&EtrueMid[type23lowOffset-1]);
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
