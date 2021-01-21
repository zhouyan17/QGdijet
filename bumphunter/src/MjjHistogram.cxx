// ---------------------------------------------------------

#include "inc/MjjHistogram.h"

// ---------------------------------------------------------
MjjHistogram::MjjHistogram(TH1D * inputHisto, bool getEffFromErrs) :
  fHistogram(*(TH1D*) inputHisto),
  fEffectiveHistogram(*(TH1D*) inputHisto),
  fWeightsHistogram(*(TH1D*) inputHisto),
  fRandomNumberGenerator(0)
{

  GetHistOutermostBinsWithData();

  TString basicname(Form("%s_mjjhist_basic",inputHisto->GetName()));
  fHistogram.SetName(basicname);
  TString effectivename(Form("%s_mjjhist_eff",inputHisto->GetName()));
  fEffectiveHistogram.SetName(effectivename);
  TString weightname(Form("%s_mjjhist_weights",inputHisto->GetName()));
  fWeightsHistogram.SetName(weightname);

  StoreEffectiveHistAndErrors(getEffFromErrs);
}

MjjHistogram::MjjHistogram(TH1D * inputHisto, int seed, bool getEffFromErrs) :
  fHistogram(*(TH1D*) inputHisto),
  fEffectiveHistogram(*(TH1D*) inputHisto),
  fWeightsHistogram(*(TH1D*) inputHisto),
  fRandomNumberGenerator(seed)
{

  GetHistOutermostBinsWithData();

  TString basicname(Form("%s_mjjhist_basic",inputHisto->GetName()));
  fHistogram.SetName(basicname);
  TString effectivename(Form("%s_mjjhist_eff",inputHisto->GetName()));
  fEffectiveHistogram.SetName(effectivename);
  TString weightname(Form("%s_mjjhist_weights",inputHisto->GetName()));
  fWeightsHistogram.SetName(weightname);

  StoreEffectiveHistAndErrors(getEffFromErrs);
}

// ---------------------------------------------------------
MjjHistogram::~MjjHistogram() 
{

}

// ---------------------------------------------------------
void MjjHistogram::GetHistOutermostBinsWithData(float epsilon) 
{
  int firstBin = 1; 
  while (fHistogram.GetBinContent(firstBin)==0 && firstBin < fHistogram.GetNbinsX()) {firstBin++; }
  int lastBin = fHistogram.GetNbinsX();
  if (epsilon < 0) while (fHistogram.GetBinContent(lastBin)==0 && lastBin > 1) lastBin--;
  else {
    while (fHistogram.GetBinContent(lastBin) < epsilon && lastBin > 1) { lastBin--; }
  }
  if (firstBin==fHistogram.GetNbinsX() && lastBin == 1) {
    std::cout << "No data in histogram! Resetting limits to first and last bin." << std::endl;
    firstBin = 1; lastBin = fHistogram.GetNbinsX();
  }
  fFirstBinWithData = firstBin;
  fLastBinWithData = lastBin;
  return;
}

// ---------------------------------------------------------
TH1D MjjHistogram::NormalizeByBinWidth()
{

  TH1D normalizedHistogram(fHistogram);
  string extra = "normalized";
  extra.append(fHistogram.GetName());
  normalizedHistogram.SetName(extra.c_str());
  normalizedHistogram.SetDefaultSumw2();

  //Take a spectrum that is in events per bin
  double width=0. , sigma=0. , dsigma=0. , events=0. , devents=0.; 
  int nBins = fHistogram.GetNbinsX()+2;
  for (int bin=0; bin < nBins; ++bin) {
    width= fHistogram.GetBinWidth(bin);
    events= fHistogram.GetBinContent(bin);
    devents= fHistogram.GetBinError(bin);
    sigma = events / width;
    dsigma = devents / width;
    normalizedHistogram.SetBinContent(bin,sigma);
    normalizedHistogram.SetBinError(bin,dsigma);
  }
  
  return normalizedHistogram;
}

// ---------------------------------------------------------
std::pair<int,int> MjjHistogram::GetIntervalContainingPercentage(double percentage) 
{

  double integral = fHistogram.Integral();
  if (integral==0) return std::make_pair(-1,-1);
  double thisPercentage=0;
  double thisInterval=0;
  int remember1=1;
  int remember2=1;
  double smallestInterval= std::numeric_limits<double>::max();
//  double rememberThisPercentage=0;
  for (int bin1=fFirstBinWithData; bin1<=fLastBinWithData; bin1 += 1) { 
    for (int bin2=bin1; bin2<=fLastBinWithData; bin2 += 1) {
      thisPercentage = (fHistogram.Integral(bin1,bin2))/integral;
      thisInterval = fHistogram.GetBinLowEdge(bin2) + fHistogram.GetBinWidth(bin2) - fHistogram.GetBinLowEdge(bin1);
      if (thisPercentage >= percentage) {
        if (thisInterval < smallestInterval) {
          remember1=bin1;
          remember2=bin2;
          smallestInterval=thisInterval;
          //rememberThisPercentage=thisPercentage;
        }
        break;
      }
    }
  }

  return std::make_pair(remember1,remember2);

}

// ---------------------------------------------------------

void MjjHistogram::StoreEffectiveHistAndErrors(bool getEffFromErrs)
{

  double Nobs, NobsErr, Neff, weff;
  for (int bin=0; bin<=fHistogram.GetNbinsX()+1; bin++) {
    Nobs = fHistogram.GetBinContent(bin);
    NobsErr = fHistogram.GetBinError(bin);
    Neff = (NobsErr!=0)?pow(Nobs/NobsErr,2):0;
    weff = (Neff!=0)? Nobs / Neff:1.;
    // Protect against future floating point errors
    if (AreIdentical(weff,double(1.0)) or !(getEffFromErrs)) {
      fEffectiveHistogram.SetBinContent(bin,Nobs);
      fEffectiveHistogram.SetBinError(bin,NobsErr);
      fWeightsHistogram.SetBinContent(bin,double(1.0));
      fWeightsHistogram.SetBinError(bin,double(0.0));
    } else {
      fEffectiveHistogram.SetBinContent(bin,Neff);
      fEffectiveHistogram.SetBinError(bin,sqrt(Neff));
      fWeightsHistogram.SetBinContent(bin,weff);
      fWeightsHistogram.SetBinError(bin,0.);
    }
  }

  return;
}

// ---------------------------------------------------------
void MjjHistogram::SetEffectiveFromBasicAndWeights(TH1D * weights) 
{
  double basicCont, basicErr, weight;
  for (int bin=0; bin<=fHistogram.GetNbinsX()+1; bin++) {
    fWeightsHistogram.SetBinContent(bin,weights->GetBinContent(bin));
    basicCont = fHistogram.GetBinContent(bin);
    basicErr = fHistogram.GetBinError(bin);
    weight = weights->GetBinContent(bin);
    // Protect against future floating point errors
    if (AreIdentical(weight,double(1.0))) {
      fEffectiveHistogram.SetBinContent(bin,basicCont);
      fEffectiveHistogram.SetBinError(bin,basicErr);
    } else {
      fEffectiveHistogram.SetBinContent(bin,basicCont/weight);
      fEffectiveHistogram.SetBinError(bin,basicErr/weight);
    }
  }
  return;
}

// ---------------------------------------------------------
void MjjHistogram::PoissonFluctuateBinByBin(TH1D * pseudoHist) 
{

  double effExp,pseudo,weight;
  for (int bin=0; bin<=fHistogram.GetNbinsX()+1; bin++) {
    effExp = fEffectiveHistogram.GetBinContent(bin);
//    effExpErr = fEffectiveHistogram.GetBinError(bin);
    weight = fWeightsHistogram.GetBinContent(bin);
    pseudo = fRandomNumberGenerator.PoissonD(effExp);
    // Protect against future floating point errors
    if (AreIdentical(weight,double(1.0))) {
      pseudoHist->SetBinContent(bin,pseudo);
      pseudoHist->SetBinError(bin,sqrt(pseudo));
    } else {
      pseudoHist->SetBinContent(bin,pseudo*weight);
      pseudoHist->SetBinError(bin,sqrt(pseudo)*weight);
    }
  }

  return;

}

// ---------------------------------------------------------
bool MjjHistogram::AreIdentical(double input, double compare) 
{

  if (fabs(compare - input) < std::numeric_limits<float>::epsilon()) return true;
  return false;

}

// ---------------------------------------------------------

