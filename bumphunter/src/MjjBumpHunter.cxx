// ---------------------------------------------------------

#include "inc/MjjBumpHunter.h"

// ---------------------------------------------------------
MjjBumpHunter::MjjBumpHunter() :
  fSpectrumTomography(), fErrHist() {

  // Set defaults as specified in header

  fAllowDeficit = false;
  fUseSidebands = false;
  fMinBinsInBump = 1;
  fMaxBinsInBump = 1e5;
  fNBinsInSideband = 1;
  fNBinsInShift = 1;
  fDoErr = false;
  fLowEdge = fHighEdge = 0;
  fLowEdgesAllBumps.clear();
  fHighEdgesAllBumps.clear();
  fProbAllBumps.clear();
  fNoisy = false;
}

// ---------------------------------------------------------
MjjBumpHunter::~MjjBumpHunter() 
   { }

// ---------------------------------------------------------
double MjjBumpHunter::DoTest(MjjHistogram & dataHist, MjjHistogram & bkgHist, int firstBinToUse, int lastBinToUse){

  TH1D h_data(dataHist.GetEffectiveHistogram());
  TString dataname(Form("%s_bh_data",h_data.GetName()));
  h_data.SetName(dataname); 

  TH1D h_bkg(bkgHist.GetEffectiveHistogram());
  TString bkgname(Form("%s_bh_bkg",h_bkg.GetName()));
  h_bkg.SetName(bkgname);

//  std::cout << "h_data:" << std::endl;
//  h_data.Print("all");
//  std::cout << "h_bkg:" << std::endl;
//  h_bkg.Print("all");

  int nBinsData = h_data.GetNbinsX();
  int nBinsBkg = h_bkg.GetNbinsX();
  assert(nBinsData==nBinsBkg);

  // Find first and last bins with data
  // If reasonable, overwrite with user's choice
  int firstBin = dataHist.GetFirstBinWithData();
  int lastBin = dataHist.GetLastBinWithData();
  if (firstBinToUse>0 && firstBinToUse > firstBin && firstBinToUse < lastBin) firstBin = firstBinToUse;
  if (lastBinToUse > firstBinToUse && lastBinToUse>0 && lastBinToUse > firstBin && lastBinToUse < lastBin) lastBin = lastBinToUse;

  vector<std::pair<int,int> > regionsdef;
  if (fExcludeWindow) {
    regionsdef.push_back(std::make_pair(firstBin,fFirstBinToExclude-1));
    regionsdef.push_back(std::make_pair(fLastBinToExclude+1,lastBin));
  } else {
    regionsdef.push_back(std::make_pair(firstBin,lastBin));
  }

  // To hold final p-value, low edge, and high edge of most interesting bump  
  fMostInterestingP = 1;
  fLowEdge = fHighEdge = 0;

  // Store low edges, high edges, p-values of all tested bumps
  fLowEdgesAllBumps.clear();
  fHighEdgesAllBumps.clear();
  fProbAllBumps.clear();

  for (unsigned int index = 0; index < regionsdef.size(); index++) {

    int thisfirstBin = regionsdef.at(index).first;
    int thislastBin = regionsdef.at(index).second;

    //std::cout << "Testing range " << thisfirstBin << ", " << thislastBin << std::endl;

    // Set number of bins in spectrum, max number of bins in bump
    int nBins = thislastBin-thisfirstBin+1;
    int minWidth = std::max(fMinBinsInBump,1);
    int maxWidth = (int) std::min(fMaxBinsInBump,nBins/2);
 
    DoTestCore(h_data, h_bkg, minWidth,maxWidth,thisfirstBin,thislastBin);

    //std::cout << "Current biggest bump has p-val "<< fMostInterestingP
    //          << " and is between " << fLowEdge << ", " << fHighEdge << std::endl;

  }

  // Store fProbAllBumps contents into tomography plot
  TGraphErrors tomography(fProbAllBumps.size());
  for (unsigned int i=0; i<fProbAllBumps.size(); i++) {
    tomography.SetPoint(i,(fLowEdgesAllBumps.at(i)+fHighEdgesAllBumps.at(i))/2., fProbAllBumps.at(i));
    tomography.SetPointError(i,(fHighEdgesAllBumps.at(i)-fLowEdgesAllBumps.at(i))/2., 0);
  }
  fSpectrumTomography = tomography;

  if (fMostInterestingP==0) FindBumpInCaseOfIncalculable(h_data, h_bkg, firstBin,lastBin);

//  std::cout << "Finished full run of BH and chose most interesting p=" << fMostInterestingP << std::endl;
  return -log(fMostInterestingP);

}

void MjjBumpHunter::DoTestCore(TH1D h_data, TH1D h_bkg, int minWidth, int maxWidth, int firstBin, int lastBin) {

  // Declare variables to hold everything during bump hunt:
  double sidebandWidth = 0;
  double smallestPforWidth, lowEdgeForWidth, highEdgeForWidth;
  // Left hand side of leftmost bin in bump at minimum of its range
  int minBinL = 0;
  // Left hand side of leftmost bin in bump at maximum of its range
  int maxBinL = 0;
  // Number of bins by which we shift window between tests
  int nbinsinstep = 0;
  // Left of bump window, right of bump window, left of left sideband, right of right sideband
  int binL , binR , binLL, binRR;
  // Background, background error in center, left sideband, right sideband
  double bC, bL, bR, deltaBC, deltaBL, deltaBR;
  // Data in center, left sideband, right sideband
  double dC, dL, dR;
  // Newly added: data error
  double deltaDC, deltaDL, deltaDR;
  // p-value in left sideband, center, right sideband
  double probL, probC, probR;

  // Initialize all.
  dC = dL = dR = 0.;
  deltaBC = deltaBL = deltaBR = 0.;
  bC = bL = bR = deltaBC = deltaBL = deltaBR = 0.;
  binL = binR = binLL = binRR = 0;
  smallestPforWidth = lowEdgeForWidth = highEdgeForWidth = 0.;
  probL = probC = probR = 1.;

  // Test each bin width in range
  int counter = 0;
  for (int width = minWidth; width < maxWidth+1; width+=1) {

    // Use user specification for sideband width if reasonable
    // Else use default of width/2 to a minimum of 1 bin
    if (fNBinsInSideband>=1) sidebandWidth = fNBinsInSideband;
    else sidebandWidth = std::max(1,(int)width/2);

    // Set number of bins by which to shift window in checks
    if (!(fNBinsInShift<1)) nbinsinstep = fNBinsInShift;
    else nbinsinstep = std::max(1,(int)width/2);

    // Reset values before looping over bins
    smallestPforWidth = 1;
    lowEdgeForWidth = 0;
    highEdgeForWidth = 0;
    if (fUseSidebands) {
       minBinL = firstBin + sidebandWidth;
       maxBinL = lastBin - width - sidebandWidth + 1;
    } else {
       minBinL = firstBin;
       maxBinL = lastBin - width + 1;
    }

    // Loop while leftmost edge of bump is below or equal to leftmost allowed limit
    for (binL = minBinL; binL <= maxBinL; binL += nbinsinstep) {
      // Find upper edge of bump, sideband edges
      binR = binL+width-1;
      binLL = binL-sidebandWidth;
      binRR = binR+sidebandWidth;

      counter++;

      //std::cout << "Doing number " << counter << ": bins " << binL << " to " << binR << std::endl;


      // Get effective data and background content and uncertainty in bump window and sidebands
      // from histogram and errors
      GetEffectiveBandContentsWithError(h_data, h_bkg, binL, binR, dC, deltaDC, bC, deltaBC);
      if (fUseSidebands) {
        GetEffectiveBandContentsWithError(h_data, h_bkg, binLL, binL-1, dL, deltaDL, bL, deltaBL);
        GetEffectiveBandContentsWithError(h_data, h_bkg, binR+1, binRR, dR, deltaDR, bR, deltaBR);
      }
      // Determine if it is an excess. Only keep if flagged to do so.
      if (!(fAllowDeficit)) { if (dC<=bC) continue; }

      // Get probabilities for observations in window and sidebands
      if (fDoErr) {
        //if (fNoisy) std::cout << "Doing number " << counter << ": bins " << binL << " to " << binR << std::endl;
        probC = PoissonConvGammaPval(dC,bC,deltaBC);
        if (fUseSidebands) {
          probL = PoissonConvGammaPval(dC,bC,deltaBC);
          probR = PoissonConvGammaPval(dC,bC,deltaBC);
        }
        if (probC < 0.0000001) {std::cout << "Error!" << std::endl;}
      } else {
        probC = PoissonPval(dC,bC);
        if (fUseSidebands) {
          probL = PoissonPval(dL,bL);
          probR = PoissonPval(dR,bR);
        }
      }

      //std::cout << "Comparing " << dC << " to " << bC << " got prob " << probC << std::endl;

      // Ignore cases where a significant discrepancy is observed in sidebands
      if (fUseSidebands) {
        if (probL <= 1e-3 || probR <= 1e-3)
          continue;
      }

      // Save information to use later in tomography plot
      fLowEdgesAllBumps.push_back(h_data.GetBinLowEdge(binL));
      fHighEdgesAllBumps.push_back(h_data.GetBinLowEdge(binR) + h_data.GetBinWidth(binR));
      fProbAllBumps.push_back(probC);
      if (probC < smallestPforWidth) {
        smallestPforWidth = probC;
        lowEdgeForWidth = h_bkg.GetBinLowEdge(binL);
        highEdgeForWidth = h_bkg.GetBinLowEdge(binR) + h_bkg.GetBinWidth(binR);
      }
    } // next center
    if (smallestPforWidth < fMostInterestingP) {
      fMostInterestingP = smallestPforWidth;
      fLowEdge = lowEdgeForWidth;
      fHighEdge = highEdgeForWidth;
    }
  } // next width

}

void MjjBumpHunter::FindBumpInCaseOfIncalculable(TH1D h_data, TH1D h_bkg, int firstBin, int lastBin) {

  vector<int> singlebinsinf;
  bool lastBinWasInf = false;
  bool allInfsConsecutive = true;
  for (int bin = firstBin; bin < lastBin+1; bin ++) {
    double D = h_data.GetBinContent(bin);
    double B = h_bkg.GetBinContent(bin);
    double thisbinpval = PoissonPval(D,B);
    if (thisbinpval==0 && D>B) {
      if (singlebinsinf.size()>0 && lastBinWasInf==false) allInfsConsecutive = false;
      singlebinsinf.push_back(bin);
      lastBinWasInf = true;
    } else {
      lastBinWasInf = false;
    }
  }
  
  if (singlebinsinf.size() > 0 && allInfsConsecutive) {
    fMostInterestingP = 0;
    fLowEdge = h_data.GetBinLowEdge(singlebinsinf.at(0));
    fHighEdge = h_data.GetBinLowEdge(singlebinsinf.at(singlebinsinf.size()-1))
                + h_data.GetBinWidth(singlebinsinf.at(singlebinsinf.size()-1));
  }

}

// ---------------------------------------------------------
void MjjBumpHunter::GetEffectiveBandContentsWithError(TH1D histData, TH1D histBkg, int firstBin, int lastBin, double& dataIntegral, double& dataError, double& bkgIntegral, double& bkgError) {

  dataIntegral = dataError = bkgIntegral = bkgError = 0;
  for (int bin=firstBin; bin<lastBin+1; bin++) {
    //std::cout << "In bin " << bin << " adding " << histData.GetBinContent(bin) << ", " << histData.GetBinError(bin) << ", " << histBkg.GetBinContent(bin) << ", " << fErrHist.GetBinContent(bin) << std::endl;
    dataIntegral += histData.GetBinContent(bin);
    dataError += histData.GetBinError(bin);
    bkgIntegral += histBkg.GetBinContent(bin);
    bkgError += fErrHist.GetBinContent(bin);
  }
  //std::cout << "Returning data=" << dataIntegral << ", bkg=" << bkgIntegral << std::endl;
  return; 
}

// ---------------------------------------------------------
vector<double> MjjBumpHunter::GetFurtherInformation() {

  vector<double> lowHighLimits;
  lowHighLimits.clear();
  lowHighLimits.push_back(fLowEdge);
  lowHighLimits.push_back(fHighEdge);
  return lowHighLimits;

}

// ---------------------------------------------------------

