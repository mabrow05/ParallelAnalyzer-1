#include <iostream>
#include <fstream>
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

#include "replay_pass2.h"

using namespace std;

int main(int argc, char *argv[])
{
  cout.setf(ios::fixed, ios::floatfield);
  cout.precision(12);

  // Prompt for filename of run numbers
  int iXeRunPeriod;
  int nRuns;
  int runList[500];
  cout << "Enter Xenon run period: " << endl;
  cin  >> iXeRunPeriod;
  cout << endl;

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
  char tempOut[500];
  //sprintf(tempOut, "position_map_%s.root", argv[1]);
  sprintf(tempOut, "position_map_%i.root", iXeRunPeriod);
  TFile *fileOut = new TFile(tempOut,"RECREATE");

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
    sprintf(tempIn, "/extern/UCNA/replay_pass2/replay_pass2_%i.root", runList[i]);
    TFile *fileIn = new TFile(tempIn, "READ");
    TTree *Tin = (TTree*)(fileIn->Get("pass2"));
    int nEvents = Tin->GetEntriesFast();
    cout << "Processing " << runList[i] << " ... " << endl;
    cout << "... nEvents = " << nEvents << endl;

    // Input variables (moved to here)
    Tin->SetBranchAddress("pmt0_pass2", &pmt_pass2[0]);
    Tin->SetBranchAddress("pmt1_pass2", &pmt_pass2[1]);
    Tin->SetBranchAddress("pmt2_pass2", &pmt_pass2[2]);
    Tin->SetBranchAddress("pmt3_pass2", &pmt_pass2[3]);
    Tin->SetBranchAddress("pmt4_pass2", &pmt_pass2[4]);
    Tin->SetBranchAddress("pmt5_pass2", &pmt_pass2[5]);
    Tin->SetBranchAddress("pmt6_pass2", &pmt_pass2[6]);
    Tin->SetBranchAddress("pmt7_pass2", &pmt_pass2[7]);
    
    Tin->SetBranchAddress("xE_pass2", &xE_pass2);
    Tin->SetBranchAddress("yE_pass2", &yE_pass2);
    Tin->SetBranchAddress("xW_pass2", &xW_pass2);
    Tin->SetBranchAddress("yW_pass2", &yW_pass2);
    
    Tin->SetBranchAddress("PID_pass2", &PID_pass2);
    Tin->SetBranchAddress("type_pass2", &type_pass2);
    Tin->SetBranchAddress("side_pass2", &side_pass2);
    Tin->SetBranchAddress("posError_pass2", &posError_pass2);

    // Loop over events
    for (int i=0; i<nEvents; i++) {
      Tin->GetEvent(i);  

      // Select Type 0 events
      if (PID_pass2 != 1) continue;
      if (type_pass2 > 0) continue;

      // Type 0 East Trigger
      int intBinX, intBinY;

      if (side_pass2 == 0) {
        intBinX = -1;
        intBinY = -1;

        // Determine position bin
        for (int m=0; m<nPosBinsX; m++) {
          if ( (xE_pass2 >= xBinLower[m]) && (xE_pass2 < xBinUpper[m]) ) intBinX = m;
        }
        for (int m=0; m<nPosBinsY; m++) {
          if ( (yE_pass2 >= yBinLower[m]) && (yE_pass2 < yBinUpper[m]) ) intBinY = m;
        }

        // Fill PMT histograms
        if (intBinX>-1 && intBinY>-1) hisxy[0][intBinX][intBinY]->Fill(pmt_pass2[0]);
        if (intBinX>-1 && intBinY>-1) hisxy[1][intBinX][intBinY]->Fill(pmt_pass2[1]);
        if (intBinX>-1 && intBinY>-1) hisxy[2][intBinX][intBinY]->Fill(pmt_pass2[2]);
        if (intBinX>-1 && intBinY>-1) hisxy[3][intBinX][intBinY]->Fill(pmt_pass2[3]);
      }

      // Type 0 West Trigger
      if (side_pass2 == 1) {
        intBinX = -1;
        intBinY = -1;
        // Determine position bin
        for (int m=0; m<nPosBinsX; m++) {
          if ( (xW_pass2 >= xBinLower[m]) && (xW_pass2 < xBinUpper[m]) ) intBinX = m;
        }
        for (int m=0; m<nPosBinsY; m++) {
          if ( (yW_pass2 >= yBinLower[m]) && (yW_pass2 < yBinUpper[m]) ) intBinY = m;
        }

	// Fill PMT histograms 
        if (intBinX>-1 && intBinY>-1) hisxy[4][intBinX][intBinY]->Fill(pmt_pass2[4]);
        if (intBinX>-1 && intBinY>-1) hisxy[5][intBinX][intBinY]->Fill(pmt_pass2[5]);
        if (intBinX>-1 && intBinY>-1) hisxy[6][intBinX][intBinY]->Fill(pmt_pass2[6]);
        if (intBinX>-1 && intBinY>-1) hisxy[7][intBinX][intBinY]->Fill(pmt_pass2[7]);
      }


    }

    // Close input ntuple
    fileIn->Close();

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
	
	  
	if (r<=50.)
	      {
		spec = new TSpectrum(20);
		Int_t npeaks = spec->Search(hisxy[p][i][j],2,"",0.5);
		
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

	double r = sqrt(xBinCenter[i]*xBinCenter[i]+yBinCenter[j]*yBinCenter[j]);

        if (maxCounts[p][i][j] > 0. && r<=50.) {
          hisxy[p][i][j]->Fit(fitName, "LRQ");
          fitMean[p][i][j] = gaussian_fit[p][i][j]->GetParameter(1);
	  //attempt at manual fix of fitting problem
	  if (fitMean[p][i][j]<xLow[p][i][j] || fitMean[p][i][j]>xHigh[p][i][j])
	    { 
	      string pmt;
	      if (p<4) pmt = "East";
	      else pmt = "West";
	      int pmtNum = p<4?p:(p-4);
	      cout << "Bad Fit in PMT " << pmt << " " << pmtNum << " at " << xBinCenter[i] << ", " << yBinCenter[j] <<  " -> " << fitMean[p][i][j] << endl;
	      Float_t lowerBoundPar1 = (float)xLow[p][i][j];
	      Float_t upperBoundPar1 = (float)xHigh[p][i][j];
	      gaussian_fit[p][i][j]->SetParameter(0,maxCounts[p][i][j]);
	      gaussian_fit[p][i][j]->SetParameter(1,binCenterMax[p][i][j]);
	      gaussian_fit[p][i][j]->SetParameter(2,100.);
	      gaussian_fit[p][i][j]->SetParLimits(1,lowerBoundPar1,upperBoundPar1);	      
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
  char tempMap[500];
  sprintf(tempMap, "position_map_%i.dat", iXeRunPeriod);
  ofstream outMap(tempMap);

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

  // Write norms to file
  char tempNorm[500];
  sprintf(tempNorm, "norm_%i.dat", iXeRunPeriod);
  ofstream outNorm(tempNorm);

  for (int p=0; p<nPMT; p++) {
    outNorm << norm[p] << endl;
  }

  // Write output ntuple
  fileOut->Write();
  fileOut->Close();

  return 0;
}
