// ---------------------------------------------------------

#include "inc/MjjLogLikelihoodTest.h"

// ---------------------------------------------------------
MjjLogLikelihoodTest::MjjLogLikelihoodTest() 
{

}

// ---------------------------------------------------------
MjjLogLikelihoodTest::~MjjLogLikelihoodTest() 
{
  
}

// ---------------------------------------------------------
double MjjLogLikelihoodTest::LogLikelihood(double d, double b)
{ 

  double thisVal = TMath::PoissonI(d,b);

  if (thisVal==0) { return -log(std::numeric_limits< double >::min()); }
  return -log(thisVal);

}

// ---------------------------------------------------------
double MjjLogLikelihoodTest::DoTest(MjjHistogram & dataHist, MjjHistogram & backgroundHist, int firstBinToUse, int lastBinToUse)
{

  TH1D normData(dataHist.GetEffectiveHistogram());
  TString dataname(Form("%s_bh_data",normData.GetName()));
  normData.SetName(dataname);

  TH1D normBkg(backgroundHist.GetEffectiveHistogram());
  TString bkgname(Form("%s_bh_bkg",normBkg.GetName()));
  normBkg.SetName(bkgname);

  double answer = 0;
  int nBinsData = normData.GetNbinsX();
  int nBinsBkg = normBkg.GetNbinsX();
  assert(nBinsBkg==nBinsData);

  int firstBin, lastBin;
  if (firstBinToUse<0 || firstBinToUse>nBinsData+1) firstBin = dataHist.GetFirstBinWithData();
  else firstBin = firstBinToUse;
  if(lastBinToUse<0 || lastBinToUse>nBinsData+1 || lastBinToUse < firstBinToUse) lastBin = dataHist.GetLastBinWithData();
  else lastBin = lastBinToUse;

  for (int bin = firstBin; bin < lastBin+1; bin++) {

    if (fExcludeWindow) {
      if (bin >= fFirstBinToExclude && bin <= fLastBinToExclude) continue;
    }
    double d = normData.GetBinContent(bin);
    double b = normBkg.GetBinContent(bin);
    answer+=LogLikelihood(d,b);
  }

  return answer;

}

// ---------------------------------------------------------

