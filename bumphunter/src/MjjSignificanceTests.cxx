// ---------------------------------------------------------

#include "inc/MjjSignificanceTests.h"

// ---------------------------------------------------------
MjjSignificanceTests::MjjSignificanceTests() 
{

}

// ---------------------------------------------------------
MjjSignificanceTests::~MjjSignificanceTests() 
{
  
}

// ---------------------------------------------------------
TH1D MjjSignificanceTests::GetRelativeDifference(MjjHistogram & dataHist, MjjHistogram & bkgHist, int firstBinToUse, int lastBinToUse)
{

  TH1D normData(dataHist.GetEffectiveHistogram());
  TString dataname(Form("%s_reldiff_data",normData.GetName()));
  normData.SetName(dataname);

  TH1D normBkg(bkgHist.GetEffectiveHistogram());
  TString bkgname(Form("%s_reldiff_bkg",normBkg.GetName()));
  normBkg.SetName(bkgname);

  int nBinsData = normData.GetNbinsX();
  int nBinsBkg = normBkg.GetNbinsX();
  assert(nBinsBkg==nBinsData);

  TH1D result(normBkg);
  TString name(Form("result_reldiff_%s",normBkg.GetName()));
  result.SetName(name);
  result.Clear();

  int firstBin, lastBin;
  if (firstBinToUse<0 || firstBinToUse>nBinsData+1) firstBin = dataHist.GetFirstBinWithData();
  else firstBin = firstBinToUse;
  if(lastBinToUse<0 || lastBinToUse>nBinsData+1 || lastBinToUse < firstBinToUse) lastBin = dataHist.GetLastBinWithData();
  else lastBin = lastBinToUse;

  for (int bin=firstBin; bin<lastBin+1; bin++) {
   
    double D = normData.GetBinContent(bin);
    double Derr = normData.GetBinError(bin);
    double B = normBkg.GetBinContent(bin);

    if (B != 0) { // Protect against divide-by-zero
      double frac = (D-B)/B;
      double fracErr = fabs(Derr/B);

      result.SetBinContent(bin,frac);
      result.SetBinError(bin,fracErr);

    }
    else {
      result.SetBinContent(bin,0);
      result.SetBinError(bin,0);
    }
  }

  return result;

}

// ---------------------------------------------------------
TH1D MjjSignificanceTests::GetRelativeDifferenceWithErrorOnBkg(MjjHistogram & dataHist, MjjHistogram & bkgHist, int firstBinToUse, int lastBinToUse) 
{

  TH1D normData(dataHist.GetEffectiveHistogram());
  TString dataname(Form("%s_reldiff_data_werr",normData.GetName()));
  normData.SetName(dataname);

  TH1D normBkg(bkgHist.GetEffectiveHistogram());
  TString bkgname(Form("%s_reldiff_bkg_werr",normBkg.GetName()));
  normBkg.SetName(bkgname);

  int nBinsData = normData.GetNbinsX();
  int nBinsBkg = normBkg.GetNbinsX();
  assert(nBinsBkg==nBinsData);

  TH1D result(normBkg);
  TString name(Form("result_reldiffwerr_%s",normBkg.GetName()));
  result.SetName(name);
  result.Clear();

  int firstBin, lastBin;
  if (firstBinToUse<0 || firstBinToUse>nBinsData+1) firstBin = dataHist.GetFirstBinWithData();
  else firstBin = firstBinToUse;
  if(lastBinToUse<0 || lastBinToUse>nBinsData+1 || lastBinToUse < firstBinToUse) lastBin = dataHist.GetLastBinWithData();
  else lastBin = lastBinToUse;

  for (int bin=firstBin; bin<lastBin+1; bin++) {

    double D = normData.GetBinContent(bin);
    double Derr = normData.GetBinError(bin);
    double B = normBkg.GetBinContent(bin);
    double Berr = normBkg.GetBinContent(bin);

    if (B != 0) { // Protect against divide-by-zero
      double frac = (D-B)/B;
      double fracErr = sqrt(pow(Derr/Berr,2)+pow(-(D*Berr)/(B*B),2));

      result.SetBinContent(bin,frac);
      result.SetBinError(bin,fracErr);

    }
    else {
      result.SetBinContent(bin,0);
      result.SetBinError(bin,0);
    }
  }

  return result;

}

// ---------------------------------------------------------
TH1D MjjSignificanceTests::GetSignificanceOfDifference(MjjHistogram & dataHist, MjjHistogram & bkgHist, int firstBinToUse, int lastBinToUse)
{

  TH1D normData(dataHist.GetEffectiveHistogram());
  TString dataname(Form("%s_sigdiff_data",normData.GetName()));
  normData.SetName(dataname);

  TH1D normBkg(bkgHist.GetEffectiveHistogram());
  TString bkgname(Form("%s_sigdiff_bkg",normBkg.GetName()));
  normBkg.SetName(bkgname);

  int nBinsData = normData.GetNbinsX();
  int nBinsBkg = normBkg.GetNbinsX();
  assert(nBinsBkg==nBinsData);

  TH1D result(normBkg);
  TString name(Form("result_sigdiff_%s",normBkg.GetName()));
  result.SetName(name);
  result.Clear();

  int firstBin, lastBin;
  if (firstBinToUse<0 || firstBinToUse>nBinsData+1) firstBin = dataHist.GetFirstBinWithData();
  else firstBin = firstBinToUse;
  if(lastBinToUse<0 || lastBinToUse>nBinsData+1 || lastBinToUse < firstBinToUse) lastBin = dataHist.GetLastBinWithData();
  else lastBin = lastBinToUse;

  for (int bin=firstBin; bin<lastBin+1; bin++) {

    double D = normData.GetBinContent(bin);
    double Derr = normData.GetBinError(bin);
    double B = normBkg.GetBinContent(bin);
    double Berr = normBkg.GetBinError(bin);

    if (B != 0) { // Matches other functions now, but since would still give value not sure why
      double frac = (D-B)/sqrt(Berr*Berr + Derr*Derr);
      double fracErr = 0.;

      result.SetBinContent(bin,frac);
      result.SetBinError(bin,fracErr);

    }
    else {
      result.SetBinContent(bin,0);
      result.SetBinError(bin,0);
    }
  }

  return result;

}

// ---------------------------------------------------------
TH1D MjjSignificanceTests::GetResidual(MjjHistogram & dataHist, MjjHistogram & bkgHist, int firstBinToUse, int lastBinToUse, TH1D * errHist)
{

  TH1D normData(dataHist.GetEffectiveHistogram());
  TString dataname(Form("%s_residual_data_werr",normData.GetName()));
  normData.SetName(dataname);

  TH1D normBkg(bkgHist.GetEffectiveHistogram());
  TString bkgname(Form("%s_residual_bkg_werr",normBkg.GetName()));
  normBkg.SetName(bkgname);

  int nBinsData = normData.GetNbinsX();
  int nBinsBkg = normBkg.GetNbinsX();
  assert(nBinsBkg==nBinsData);

  TH1D result(normBkg);
  TString name(Form("result_residual_%s",normBkg.GetName()));
  result.SetName(name);
  result.Clear();

  int firstBin, lastBin;
  if (firstBinToUse<0 || firstBinToUse>nBinsData+1) firstBin = dataHist.GetFirstBinWithData();
  else firstBin = firstBinToUse;
  if(lastBinToUse<0 || lastBinToUse>nBinsData+1 || lastBinToUse < firstBinToUse) lastBin = dataHist.GetLastBinWithData();
  else lastBin = lastBinToUse;

  for (int bin=0; bin<result.GetNbinsX()+2; bin++) {
    if (bin>firstBin-1 && bin<lastBin+1) {

      double D = normData.GetBinContent(bin);
      double B = normBkg.GetBinContent(bin);
      double bErr = normBkg.GetBinError(bin);

      if (D != 0) { // Don't want residual Getting plotted outside of bins that actually contain something
        double PVal;
        if (errHist){ PVal = PoissonConvGammaPval(D,B,bErr); }
        else PVal = PoissonPval(D,B);
        double frac;
        frac = ProbToSigma(PVal);
        // Trim to a reasonable size for display
        if (frac > 100) frac = 20;
        // use poissonPvalNonNegative
        if (frac<0.) frac = 0.; //if it's negative (very insignificant), then set it to zero.
        if (D < B) frac *= -1; //make it signed.
  
        result.SetBinContent(bin,frac);
        result.SetBinError(bin,0);
        
//        std::cout << "Residual result in bin " << bin << " is " << frac << " from d, b " << D << ", " << B << std::endl;


      }
      else {
        result.SetBinContent(bin,0);
        result.SetBinError(bin,0);
      }

    } else {
      result.SetBinContent(bin,0);
      result.SetBinError(bin,0);
    }
  }

  return result;

}

// ---------------------------------------------------------

