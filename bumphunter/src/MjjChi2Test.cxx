// ---------------------------------------------------------

#include "inc/MjjChi2Test.h"

// ---------------------------------------------------------
MjjChi2Test::MjjChi2Test() 
{

}

// ---------------------------------------------------------
MjjChi2Test::~MjjChi2Test() 
{
  
}

// ---------------------------------------------------------
double MjjChi2Test::DoTest(MjjHistogram & dataHist, MjjHistogram & bkgHist, int firstBinToUse, int lastBinToUse)
{

  // bkgHist is the expectation (background), dataHist the observed data.
  // chi2 = sum_bins of (data - bkg)^2 / (sqrt(sqrt(data)^2+deltaBkg^2)^2
  TH1D normData(dataHist.GetEffectiveHistogram());
  TString dataname(Form("%s_bh_data",normData.GetName()));
  normData.SetName(dataname);

  TH1D normBkg(bkgHist.GetEffectiveHistogram());
  TString bkgname(Form("%s_bh_bkg",normBkg.GetName()));
  normBkg.SetName(bkgname);

  int nBinsData = normData.GetNbinsX();
  int nBinsBkg = normBkg.GetNbinsX();
  assert(nBinsBkg==nBinsData);

  int firstBin, lastBin;
  if (firstBinToUse<0 || firstBinToUse>nBinsData+1) firstBin = dataHist.GetFirstBinWithData();
  else firstBin = firstBinToUse;
  if(lastBinToUse<0 || lastBinToUse>nBinsData+1 || lastBinToUse < firstBinToUse) lastBin = dataHist.GetLastBinWithData();
  else lastBin = lastBinToUse;

  double answer = 0;
  for (int bin = firstBin; bin < lastBin+1; bin++) {

    if (fExcludeWindow) {
      if (bin >= fFirstBinToExclude && bin <= fLastBinToExclude) continue;
    }

    double d = normData.GetBinContent(bin);
    if (d==0) continue;
    double b = normBkg.GetBinContent(bin);
    double deltaB = normBkg.GetBinError(bin);

    // Definitions of chi2 itself
//    double term = (d - b) / sqrt(b); //Pearson's Chi is defined with sqrt(b), and only for populated bins.
    double term = (d - b) / sqrt(b+deltaB*deltaB); //give it an extra uncertainty in the background if it's there
    answer += (term*term);

//    if (fExcludeWindow) std::cout << "Chi2: Including bin " << bin << " with contribution " << term*term << " from d, b " << d << ", " << b << std::endl;


  }

  return answer;

}

// ---------------------------------------------------------

