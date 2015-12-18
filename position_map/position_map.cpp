//SWANKS EDIT
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <string>

// ROOT libraries
#include "TRandom3.h"
#include <TTree.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TF1.h>
#include <TH1D.h>
#include <TSpectrum.h>

#include "fullTreeVariables.h"
#include "MWPCGeometry.h"
#include "pedestals.h"
#include "cuts.h"
#include "basic_reconstruction.h"
#include "DataTree.hh"

#include "replay_pass2.h"

using namespace std;

string itos(int val) {
  char temp[32];
  sprintf(temp,"%i",val);
  return string(temp);
};

int main(int argc, char *argv[])
{
  cout.setf(ios::fixed, ios::floatfield);
  cout.precision(12);

  // Prompt for filename of run numbers
  int iXeRunPeriod;
  bool allResponseClasses = true;
  int numResponseClasses = 0;
  vector <int> responseClasses;
  int nRuns;
  int runList[500];
  cout << "Enter Xenon run period: " << endl;
  cin  >> iXeRunPeriod;
  cout << endl;

  cout << "All Response Classes? (true=1/false=0): " << endl;
  cin  >> allResponseClasses;
  cout << endl;

  if (!allResponseClasses) {
    
    cout << "Enter number of response classes: " << endl;
    cin >> numResponseClasses;
    cout << endl;

    if (numResponseClasses<1 || numResponseClasses>9) {
      cout << "Bad number of response classes to include!!\n";
      exit(1);
    }
    responseClasses.resize(numResponseClasses,0);
    char quest[30];
    for (int i=0; i<numResponseClasses; i++) {
      sprintf(quest,"Enter class #%i: ",i+1);
      cout << quest;
      cin >> responseClasses[i];
      cout << endl;
      
      if (responseClasses[i]<0 || responseClasses[i]>8) {
	cout << "You entered a non-existent response class!\n";
	exit(1);
      }
    }
  }


  char temp[500];
  sprintf(temp, "xenon_runs_%i.dat", iXeRunPeriod);
  ifstream fileRuns(temp);

  fileRuns >> nRuns;
  for (int i=0; i<nRuns; i++) {
    fileRuns >> runList[i];
  }

  cout << "... Number of runs: " << nRuns << endl;
  for (int i=0; i<nRuns; i++) {
    cout << runList[i] << endl;
  }

  // Position bins
  double xBinWidth = 2.5;
  double yBinWidth = 2.5;
  int nPosBinsX = 43; //Be sure these two are odd (10->11 and 5->21 and 4->27 and 2->53)
  int nPosBinsY = 43;
  double xBinLower[nPosBinsX];
  double xBinUpper[nPosBinsX];
  double xBinCenter[nPosBinsX];
  double yBinLower[nPosBinsY];
  double yBinUpper[nPosBinsY];
  double yBinCenter[nPosBinsY];
  int intXBinCenter[nPosBinsX];
  int intYBinCenter[nPosBinsY];

  for (int k=0; k<nPosBinsX; k++) {
    xBinLower[k]     = -(double)nPosBinsX*xBinWidth/2. + ((double) k)*xBinWidth;
    xBinUpper[k]     = -(double)nPosBinsX*xBinWidth/2. + ((double) k)*xBinWidth + xBinWidth;
    xBinCenter[k]    = (xBinLower[k] + xBinUpper[k])/2.;
    intXBinCenter[k] = (int) xBinCenter[k];
    //cout << xBinLower[k] << " " << intXBinCenter[k] << " " << xBinUpper[k] << endl;
  }

  for (int k=0; k<nPosBinsY; k++) {
    yBinLower[k]     = -(double)nPosBinsY*yBinWidth/2. + ((double) k)*yBinWidth;
    yBinUpper[k]     = -(double)nPosBinsY*yBinWidth/2. + ((double) k)*yBinWidth + yBinWidth;
    yBinCenter[k]    = (yBinLower[k] + yBinUpper[k])/2.;
    intYBinCenter[k] = (int) yBinCenter[k];
    //cout << yBinLower[k] << " " << intYBinCenter[k] << " " << yBinUpper[k] << endl;
  }

  // Open output ntuple
  string tempOutBase;
  string tempOut;
  //sprintf(tempOut, "position_map_%s.root", argv[1]);
  tempOutBase = "position_map_" + itos(iXeRunPeriod);
  if (!allResponseClasses) {
    tempOutBase+="_RC_";
    for (int i=0; i< numResponseClasses; i++) {
      tempOutBase+=itos(responseClasses[i]);
    }
  }
  tempOut = tempOutBase+".root";
  TFile *fileOut = new TFile(tempOut.c_str(),"RECREATE");

  // Output histograms
  int nPMT = 8;
  int nBin = 800;

  TH1F *hisxy[nPMT][nPosBinsX][nPosBinsY];
  char *hisxyName = new char[10];

  for (int p=0; p<nPMT; p++) {
    for (int i=0; i<nPosBinsX; i++) {
      for (int j=0; j<nPosBinsY; j++) {
        if (p == 0)
          sprintf(hisxyName, "e0_%i_%i", intXBinCenter[i], intYBinCenter[j]);
        if (p == 1)
          sprintf(hisxyName, "e1_%i_%i", intXBinCenter[i], intYBinCenter[j]);
        if (p == 2)
          sprintf(hisxyName, "e2_%i_%i", intXBinCenter[i], intYBinCenter[j]);
        if (p == 3)
          sprintf(hisxyName, "e3_%i_%i", intXBinCenter[i], intYBinCenter[j]);
        if (p == 4)
          sprintf(hisxyName, "w0_%i_%i", intXBinCenter[i], intYBinCenter[j]);
        if (p == 5)
          sprintf(hisxyName, "w1_%i_%i", intXBinCenter[i], intYBinCenter[j]);
        if (p == 6)
          sprintf(hisxyName, "w2_%i_%i", intXBinCenter[i], intYBinCenter[j]);
        if (p == 7)
          sprintf(hisxyName, "w3_%i_%i", intXBinCenter[i], intYBinCenter[j]);

        hisxy[p][i][j] = new TH1F(hisxyName, "", nBin,0.0,4000.0);

      }
    }
  }

  // Loop through input ntuples
  char tempIn[500];
  for (int i=0; i<nRuns; i++) {

    // Open input ntuple
    sprintf(tempIn, "%s/replay_pass2_%i.root",getenv("REPLAY_PASS2"), runList[i]);
    DataTree *t = new DataTree();
    t->setupInputTree(std::string(tempIn),"pass2");
    //TFile *fileIn = new TFile(tempIn, "READ");
    //TTree *Tin = (TTree*)(fileIn->Get("pass2"));
    int nEvents = t->getEntries();
    cout << "Processing " << runList[i] << " ... " << endl;
    cout << "... nEvents = " << nEvents << endl;


    // Loop over events
    for (int i=0; i<nEvents; i++) {
      t->getEvent(i);  

      // Select Type 0 events
      if (t->PID != 1) continue;
      if (t->Type > 0) continue;

	  
		

      // Type 0 East Trigger
      int intBinX, intBinY; 
      bool moveOnX = true, moveOnY=true; // Determining if the event is of the correct response class in x and y
      if (t->Side == 0) {
	//Swank addition: Wire Chamber Response class. 
	for (int j=0; j<numResponseClasses; j++) {
	  if (t->xeRC == responseClasses[j]) {moveOnX=false;}
	  if (t->yeRC == responseClasses[j]) {moveOnY=false;}
	}
      
	if (moveOnX || moveOnY) continue;

	intBinX = -1;
        intBinY = -1;

        // Determine position bin
        for (int m=0; m<nPosBinsX; m++) {
          if ( (t->xE.center >= xBinLower[m]) && (t->xE.center < xBinUpper[m]) ) intBinX = m;
        }
        for (int m=0; m<nPosBinsY; m++) {
          if ( (t->yE.center >= yBinLower[m]) && (t->yE.center < yBinUpper[m]) ) intBinY = m;
        }

        // Fill PMT histograms
        if (intBinX>-1 && intBinY>-1) hisxy[0][intBinX][intBinY]->Fill(t->ScintE.q1);
        if (intBinX>-1 && intBinY>-1) hisxy[1][intBinX][intBinY]->Fill(t->ScintE.q2);
        if (intBinX>-1 && intBinY>-1) hisxy[2][intBinX][intBinY]->Fill(t->ScintE.q3);
        if (intBinX>-1 && intBinY>-1) hisxy[3][intBinX][intBinY]->Fill(t->ScintE.q4);
      }

      // Type 0 West Trigger
      moveOnX=moveOnY=true;
      if (t->Side == 1) {
	//Swank Only Allow triangles!!!	  	
	for (int j=0; j<numResponseClasses; j++) {
	  if (t->xwRC == responseClasses[j]) {moveOnX=false;}
	  if (t->ywRC == responseClasses[j]) {moveOnY=false;}
	}
      
	if (moveOnX || moveOnY) continue;
	
        intBinX = -1;
        intBinY = -1;
        // Determine position bin
        for (int m=0; m<nPosBinsX; m++) {
          if ( (t->xW.center >= xBinLower[m]) && (t->xW.center < xBinUpper[m]) ) intBinX = m;
        }
        for (int m=0; m<nPosBinsY; m++) {
          if ( (t->yW.center >= yBinLower[m]) && (t->yW.center < yBinUpper[m]) ) intBinY = m;
        }

	// Fill PMT histograms 
        if (intBinX>-1 && intBinY>-1) hisxy[4][intBinX][intBinY]->Fill(t->ScintW.q1);
        if (intBinX>-1 && intBinY>-1) hisxy[5][intBinX][intBinY]->Fill(t->ScintW.q2);
        if (intBinX>-1 && intBinY>-1) hisxy[6][intBinX][intBinY]->Fill(t->ScintW.q3);
        if (intBinX>-1 && intBinY>-1) hisxy[7][intBinX][intBinY]->Fill(t->ScintW.q4);
      }


    }

    // Close input ntuple
    delete t;

  }

  // Gaussian fits to "200 keV" peak

  // Define fit ranges
  double xLow[nPMT][nPosBinsX][nPosBinsY];
  double xHigh[nPMT][nPosBinsX][nPosBinsY];
  int maxBin[nPMT][nPosBinsX][nPosBinsY];
  double maxCounts[nPMT][nPosBinsX][nPosBinsY];
  double binCenterMax[nPMT][nPosBinsX][nPosBinsY];

  TSpectrum *spec;

  for (int p=0; p<nPMT; p++) {
    for (int i=0; i<nPosBinsX; i++) {
      for (int j=0; j<nPosBinsY; j++) {	
	
	double r = sqrt(yBinCenter[j]*yBinCenter[j]+xBinCenter[i]*xBinCenter[i]);
        // Find bin with maximum content
        maxBin[p][i][j] = hisxy[p][i][j]->GetMaximumBin();
        maxCounts[p][i][j] = hisxy[p][i][j]->GetBinContent(maxBin[p][i][j]);
        binCenterMax[p][i][j] = hisxy[p][i][j]->GetBinCenter(maxBin[p][i][j]);
	
	  
	if (r<=52.5)
	      {
		spec = new TSpectrum(20);
		Int_t npeaks = spec->Search(hisxy[p][i][j],1.5,"",0.5);
		
		if (npeaks==0)
		  {
		    cout << "No peaks identified at PMT" << p << " position " << xBinCenter[i] << ", " << yBinCenter[j] << endl;
		  }
		else
		  {
		    Float_t *xpeaks = spec->GetPositionX();
		    TAxis *xaxis = (TAxis*)(hisxy[p][i][j]->GetXaxis());
		    Int_t peakBin=0;
		    Double_t BinSum=0.;
		    Double_t BinSumHold = 0.;
		    Int_t maxPeak=0.;
		    for (int pk=0;pk<npeaks;pk++) {
		      peakBin = xaxis->FindBin(xpeaks[pk]);
		      //Sum over 3 center bins of peak and compare to previos BinSum to see which peak is higher
		      BinSum = hisxy[p][i][j]->GetBinContent(peakBin) + hisxy[p][i][j]->GetBinContent(peakBin-1) + hisxy[p][i][j]->GetBinContent(peakBin+1);
		      if (BinSum>BinSumHold) {
			BinSumHold=BinSum;
			maxPeak=pk;
		      }
		    }
		    binCenterMax[p][i][j] = xpeaks[maxPeak];
		  }
		delete spec;
	      }
	xLow[p][i][j] = binCenterMax[p][i][j]/3.;
	xHigh[p][i][j] = 2.*binCenterMax[p][i][j];

	/*bool flag = false;
        // Define fit range
        for (int m=maxBin[p][i][j]; m<nBin; m++) {
          if (hisxy[p][i][j]->GetBinContent(m+1) < 0.5*maxCounts[p][i][j]) {
            xHigh[p][i][j] = hisxy[p][i][j]->GetBinCenter(m+1);
	    
	    if (((m+1)-maxBin[p][i][j])>4) break;
	    else {
	      if (r<50.) 
		{
		  flag=true;
		  cout << "Max Bin was " << m+1 << " at pmt " << p << " position " << xBinCenter[i] << ", " << yBinCenter[j] << endl;
		}
	      //xHigh[p][i][j] = 2.*binCenterMax[p][i][j];
	      break;
	    }
	  }
        }

        for (int m=maxBin[p][i][j]; m>0; m--) {
          if (hisxy[p][i][j]->GetBinContent(m-1) < 0.5*maxCounts[p][i][j]) {
            xLow[p][i][j] = hisxy[p][i][j]->GetBinCenter(m-1);
	    
	    if (m!=1 && (maxBin[p][i][j]-(m-1))>4) break;
	    else {
	      if (r<50.) 
		{
		  flag=true;
		  cout << "Min Bin was " << m-1 << " at pmt " << p << " position " << xBinCenter[i] << ", " << yBinCenter[j] << endl;
		}
	      //xLow[p][i][j] = binCenterMax[p][i][j]/3.;
	      break;
	    }
	  }
	}
	
	//xLow[p][i][j] = binCenterMax[p][i][j]/2.;
	//xHigh[p][i][j] = 3.*binCenterMax[p][i][j]/2.;
	if (flag)
	  {
	    spec = new TSpectrum(5);
	    Int_t npeaks = spec->Search(hisxy[p][i][j],1,"",0.99);
	    
	    if (npeaks==1)
	      {
		Float_t *xpeaks = spec->GetPositionX();
		Int_t max = 0;
		binCenterMax[p][i][j] = xpeaks[max];
		xLow[p][i][j] = binCenterMax[p][i][j]/3.;
		xHigh[p][i][j] = 2.*binCenterMax[p][i][j];
		cout << "No peaks identified at PMT" << p << " position " << xBinCenter[i] << ", " << yBinCenter[j] << endl;		    
	      }
	    //else if (npeaks>1) cout << "More than 1 peak identified at PMT" << p << " position " << xBinCenter[i] << ", " << yBinCenter[j] << endl;
	    else
	      {
		xLow[p][i][j] = binCenterMax[p][i][j]/3.;
		xHigh[p][i][j] = 2.*binCenterMax[p][i][j];
	      }
	    delete spec;
	    }*/
      }
    }
  }
  
  // Fit histograms
  TF1 *gaussian_fit[nPMT][nPosBinsX][nPosBinsY];
  double fitMean[nPMT][nPosBinsX][nPosBinsY];

  double sigmaMax[8] = {400.,250.,250.,250.,250.,300.,300.,500};

  for (int p=0; p<nPMT; p++) {
    for (int i=0; i<nPosBinsX; i++) {
      for (int j=0; j<nPosBinsY; j++) {

        char fitName[500];
        sprintf(fitName, "gaussian_fit_%i_%i_%i.root", p, i, j);

        gaussian_fit[p][i][j] = new TF1(fitName,
       	     	                        "[0]*exp(-(x-[1])*(x-[1])/(2.*[2]*[2]))",
                                        xLow[p][i][j], xHigh[p][i][j]);
        gaussian_fit[p][i][j]->SetParameter(0,maxCounts[p][i][j]);
        gaussian_fit[p][i][j]->SetParameter(1,binCenterMax[p][i][j]);
        gaussian_fit[p][i][j]->SetParameter(2,100.);
	gaussian_fit[p][i][j]->SetParLimits(0,0.,5000.);
	gaussian_fit[p][i][j]->SetParLimits(1,0.,4000.);	
	gaussian_fit[p][i][j]->SetParLimits(2,0.,500.);

	double r = sqrt(xBinCenter[i]*xBinCenter[i]+yBinCenter[j]*yBinCenter[j]);

        if (maxCounts[p][i][j] > 0. && r<=52.5) {
          hisxy[p][i][j]->Fit(fitName, "LRQ");
          fitMean[p][i][j] = gaussian_fit[p][i][j]->GetParameter(1);
	  //attempt at manual fix of fitting problem. Constrain p1 and p2
	  if (fitMean[p][i][j]<xLow[p][i][j] || fitMean[p][i][j]>xHigh[p][i][j])
	    { 
	      string pmt;
	      if (p<4) pmt = "East";
	      else pmt = "West";
	      int pmtNum = p<4?p:(p-4);
	      cout << "Bad Fit in PMT " << pmt << " " << pmtNum << " at " << xBinCenter[i] << ", " << yBinCenter[j] <<  " -> " << fitMean[p][i][j] << endl;
	      Float_t lowerBoundPar1 = (float)xLow[p][i][j];
	      Float_t upperBoundPar1 = (float)xHigh[p][i][j];
	      Float_t lowerBoundPar2 = 0.;
	      Float_t upperBoundPar2 = sigmaMax[p];
	      gaussian_fit[p][i][j]->SetParameter(0,maxCounts[p][i][j]);
	      gaussian_fit[p][i][j]->SetParameter(1,binCenterMax[p][i][j]);
	      gaussian_fit[p][i][j]->SetParameter(2,100.);
	      gaussian_fit[p][i][j]->SetParLimits(1,lowerBoundPar1,upperBoundPar1);
	      gaussian_fit[p][i][j]->SetParLimits(2,lowerBoundPar2,upperBoundPar2);
	      hisxy[p][i][j]->Fit(fitName, "LRQ");
	      fitMean[p][i][j] = gaussian_fit[p][i][j]->GetParameter(1);
	      cout << "New Fit Mean -> " << fitMean[p][i][j] << endl;

	      //Check no. 2. Attempt to constrain p0 
	      if (fitMean[p][i][j]<(1.2*xLow[p][i][j]) || fitMean[p][i][j]>(0.8*xHigh[p][i][j]))
		{
		  Float_t lowerBoundPar0 = (float)(0.8*maxCounts[p][i][j]);
		  Float_t upperBoundPar0 = (float)(2.*maxCounts[p][i][j]);
		  gaussian_fit[p][i][j]->SetParameter(0,maxCounts[p][i][j]);
		  gaussian_fit[p][i][j]->SetParameter(1,binCenterMax[p][i][j]);
		  gaussian_fit[p][i][j]->SetParameter(2,100.);
		  gaussian_fit[p][i][j]->SetParLimits(0,lowerBoundPar0,upperBoundPar0);
		  gaussian_fit[p][i][j]->SetParLimits(1,lowerBoundPar1,upperBoundPar1);	
		  gaussian_fit[p][i][j]->SetParLimits(2,lowerBoundPar2,upperBoundPar2);
		  hisxy[p][i][j]->Fit(fitName, "LRQ");
		  fitMean[p][i][j] = gaussian_fit[p][i][j]->GetParameter(1);
		  cout << "New Fit Mean -> " << fitMean[p][i][j] << endl;
		
		  //Checking if fit mean is still out of range. If so, set value to max bin
		  if (fitMean[p][i][j]<xLow[p][i][j] || fitMean[p][i][j]>xHigh[p][i][j])
		    {
		      fitMean[p][i][j] = binCenterMax[p][i][j];
		      cout << "**** replaced fit mean with max bin " << fitMean[p][i][j] << endl;
		    }
		}
	    }	      
        }
        else {
          fitMean[p][i][j] = 0.;
        }
	//cout << "PMT " << p << " at " << xBinCenter[i] << ", " << yBinCenter[j] << " fitMean = " << fitMean[p][i][j] << endl;

      }
    }
  }

  // Extract position maps
  double norm[nPMT];
  for (int p=0; p<nPMT; p++) {
    norm[p] = fitMean[p][nPosBinsX/2][nPosBinsY/2];
    cout << norm[p] << endl;
  }

  double positionMap[nPMT][nPosBinsX][nPosBinsY];
  for (int p=0; p<nPMT; p++) {
    for (int i=0; i<nPosBinsX; i++) {
      for (int j=0; j<nPosBinsY; j++) {
        positionMap[p][i][j] = fitMean[p][i][j] / norm[p];
      }
    }
  }

  // Write position maps to file
  string tempMap;
  tempMap = tempOutBase + ".dat";
  ofstream outMap(tempMap.c_str());

  for (int i=0; i<nPosBinsX; i++) {
    for (int j=0; j<nPosBinsY; j++) {
      outMap << xBinCenter[i]        << "  "
             << yBinCenter[j]        << "  "
             << positionMap[0][i][j] << "  "
             << positionMap[1][i][j] << "  "
             << positionMap[2][i][j] << "  "
             << positionMap[3][i][j] << "  "
             << positionMap[4][i][j] << "  "
             << positionMap[5][i][j] << "  "
             << positionMap[6][i][j] << "  "
             << positionMap[7][i][j] << endl;
    }
  }
  outMap.close();

  // Write norms to file
  string tempNorm;
  tempNorm = "norm_"+tempMap;
  ofstream outNorm(tempNorm.c_str());

  for (int p=0; p<nPMT; p++) {
    outNorm << norm[p] << endl;
  }
  outNorm.close();

  // Write output ntuple
  fileOut->Write();
  fileOut->Close();

  return 0;
}
