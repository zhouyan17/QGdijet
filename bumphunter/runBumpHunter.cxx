#include "TH1.h"
#include "TF1.h"
#include <iostream>
#include <vector>
#include <fstream>
#include <iomanip>
#include <math.h>
#include "inc/MjjHistogram.h"
#include "inc/MjjStatisticalTest.h"
#include "inc/MjjPseudoExperimenter.h"
#include "inc/MjjStatisticsBundle.h"
#include "inc/MjjSignificanceTests.h"
#include "inc/MjjBumpHunter.h"
#include "inc/MjjChi2Test.h"
#include "inc/MjjLogLikelihoodTest.h"

#include "TEnv.h"
#include "TFile.h"
#include "TH1D.h"
#include "TString.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TVector.h"
#include "TStopwatch.h"
#include "TLine.h"
#include "TLegend.h"
#include "TSystem.h"

using namespace std;
int main (int argc,char **argv)
{
  //------------------------------------------
  // Start counting time
  TStopwatch totaltime;
  totaltime.Start();

  float minBHMass = -1;
  float maxBHMass = -1;
  TString inputFileName = "" ;
  TString outPath = "";
  TString outFileName = "";
  TString dataHistoName = "" ;
  TString bkgHistoName = "";
  int nPseudoExpBH = 1e3;
  //------------------------------------------

  // Start reading input configuration
  int ip=1;
  while (ip<argc)  
  {
    if (string(argv[ip]).substr(0,2)=="--") 
    {
      //input file
      if (string(argv[ip])=="--inFile") 
      {
        if (ip+1<argc && string(argv[ip+1]).substr(0,2)!="--") 
        {
          inputFileName=argv[ip+1];
          ip+=2;
        } 
        else 
        {
          std::cout<<"\nno input file name inserted"<<std::endl; 
          break;
        }
      }
      // out path
      else if(string(argv[ip])=="--outPath")
      {
        if (ip+1<argc && string(argv[ip+1]).substr(0,2)!="--")
        {
          outPath=argv[ip+1];
          ip+=2;
        }
        else
        {
          std::cout<<"\nno name of plot to be created not included"<<std::endl;
          break ;
        }
      }
      //output file name
      else if (string(argv[ip])=="--outFileName") 
      {
        if (ip+1<argc && string(argv[ip+1]).substr(0,2)!="--") 
        {
          outFileName=argv[ip+1];
          ip+=2;
        } 
        else 
        {
          std::cout<<"\nno name of plot to be created not included"<<std::endl; 
          break ;
        }
      }
      //data histogram name
      else if (string(argv[ip])=="--dataHist") 
      {
        if (ip+1<argc && string(argv[ip+1]).substr(0,2)!="--") 
        {
          dataHistoName=argv[ip+1];
          ip+=2;
        } 
        else 
        {
          std::cout<<"\nno histogram name inserted"<<std::endl; 
          break;
        }
      }
      //background prediction histogram name
      else if (string(argv[ip])=="--bkgHist") 
      {
        if (ip+1<argc && string(argv[ip+1]).substr(0,2)!="--") 
        {
          bkgHistoName=argv[ip+1];
          ip+=2;
        } 
        else 
        {
          std::cout<<"\nno histogram name inserted"<<std::endl; 
          break;
        }
      }
      //minimum mass for Bump Hunter, if you don't want to examine all bins
      else if (string(argv[ip])=="--minBH") 
      {
        if (ip+1<argc && string(argv[ip+1]).substr(0,2)!="--") 
        {
          minBHMass=std::stof(string(argv[ip+1]));
          ip+=2;
        } 
        else 
        {
          std::cout<<"\nNo BumpHunter minimum value given "<<std::endl; 
          break ;
        }
      }
      //maximum mass for Bump Hunter, if you don't want to examine all bins
      else if (string(argv[ip])=="--maxBH")  
      {
        if (ip+1<argc && string(argv[ip+1]).substr(0,2)!="--") 
        {
          maxBHMass=std::stof(string(argv[ip+1]));
          ip+=2;
        } 
        else 
        {
          std::cout<<"\nNo BumpHunter maximum value given "<<std::endl; 
          break;
        }
      }
      //number of pseudoexperiments to use in BH (default 1000)
      else if (string(argv[ip])=="--nPseudoExpBH") 
      {
        if (ip+1<argc && string(argv[ip+1]).substr(0,2)!="--") 
        {
          nPseudoExpBH=std::stoi(string(argv[ip+1]));
          ip+=2;
        } 
        else 
        {
          std::cout<<"\nNumber of pseudoexperiments not specified "<<std::endl; 
          break ;
        }
      }
      //unknown command
      else 
      {
        std::cout<<"\nSearchPhase: command '"<<string(argv[ip])<<"' unknown"<<std::endl;
        if (ip+1<argc && string(argv[ip+1]).substr(0,2)!="--") ip+=2;
        else ip+=1;
      } 
    }//end if "--command"
    else 
    { //if command does not start with "--"
      std::cout << "\nSearchPhase: command '"<<string(argv[ip])<<"' unknown"<<std::endl;
      break;
    }//end if "--"
  }//end while loop

  //------------------------------------------
  // Retrieve histograms and set up for use
  // Open files
  //TFile * infile = TFile::Open(inputFileName,"READ");
  TFile * infile = new TFile(inputFileName,"READ");
  if(gSystem->AccessPathName(inputFileName))
  {
    std::cout<<inputFileName<<" doesn't exist."<<std::endl ;
    exit(0) ;
  }
  TH1::AddDirectory(kFALSE);

  // Get and store histograms
  TH1D * rawDataHisto = (TH1D*) infile->Get(dataHistoName);
  if(!rawDataHisto)
  {
    std::cout<<dataHistoName<<" doesn't exist."<<std::endl ;
    exit(0) ;
  }
  TH1D * rawBkgHisto = (TH1D*) infile->Get(bkgHistoName);
  if(!rawBkgHisto)
  {
    std::cout<<bkgHistoName<<" doesn't exist."<<std::endl ;
    exit(0) ;
  }

  // check out path
  if(gSystem->AccessPathName(outPath))
    gSystem->Exec("mkdir -p "+outPath) ;

  // Make histogram wrapper class that the code requires.
  // The "false" in the background histogram constructor just
  // tells the code not to assume its effective number of events
  // is related to the error bars (sorry for inconvenience)
  MjjHistogram dataHistogram(rawDataHisto, false);
  MjjHistogram backgroundHistogram(rawBkgHisto, false);
  
  //------------------------------------------
  // Set up range to use. You may not always want to bump-hunt
  // the entire set of filled bins (say if your background is poorly modeled
  // below some point.
  // This script will use the first and last bins with data to define the range
  // unless instructed otherwise.

  //int firstBin = dataHistogram.GetFirstBinWithData();
  int lastBin = dataHistogram.GetLastBinWithData();
  int firstBin = backgroundHistogram.GetFirstBinWithData();
  int firstBinBH; 
  int lastBinBH;

  if( minBHMass == -1) 
    firstBinBH = firstBin;
  else 
    firstBinBH = rawDataHisto->FindBin(minBHMass);
  if( maxBHMass == -1) 
    lastBinBH = lastBin;
  else 
    lastBinBH = rawDataHisto->FindBin(maxBHMass);

  std::cout << "Will bump hunt the spectrum in bins [" << firstBinBH << " - " << lastBinBH << "]\n" << "\tcorresponding to a range [" << rawDataHisto->GetBinLowEdge(firstBinBH) << " - "<< rawDataHisto->GetBinLowEdge(lastBinBH)+rawDataHisto->GetBinWidth(lastBinBH) << "]" << std::endl;

  //------------------------------------------
  // Set up tools for calculation.

  // Create bump hunter.
  MjjBumpHunter theBumpHunter;

  // Make statistical tests: Chi2 and Likelihood
  MjjChi2Test theChi2Test;
  MjjLogLikelihoodTest theLikelihoodTest;

  // Make a vector to store them.
  vector<MjjStatisticalTest*> theStatsTests;
  theStatsTests.push_back(&theLikelihoodTest);
  theStatsTests.push_back(&theChi2Test);
  theStatsTests.push_back(&theBumpHunter);


  // We will not consider one-bin bumps. Change this if
  // you want tigher or stricter minimum widths, but 2 is usually fine.
  theBumpHunter.SetMinBumpWidth(4);
  
  // The default maximum bin width is 1/2 the number of bins in the spectrum.
  // If you want another maximum, uncomment this and set it like so.
  // theBumpHunter.SetMaxBumpWidth(15)
  
  // The sideband option is rarely used in ATLAS but has been
  // included in BH implementations for a long time.
  // The default in the current code is "false" but you can
  // control it here
  theBumpHunter.SetUseSidebands(false);
  
  // Un-set any exclusions and bump hunt the spectrum
  theBumpHunter.SetUseWindowExclusion(false);
  double bumpHunterStat1 = theBumpHunter.DoTest(dataHistogram,backgroundHistogram,firstBinBH,lastBinBH);
  vector<double> bumpEdges1 = theBumpHunter.GetFurtherInformation();
  double lowEdgeOfBump1 = bumpEdges1.at(0);
  double highEdgeOfBump1 = bumpEdges1.at(1);
  std::cout << "BumpHunter results: stat = " << bumpHunterStat1 << std::endl;
  std::cout << "Low edge, high edge of bump: " << lowEdgeOfBump1 << " " << highEdgeOfBump1 << std::endl;
  // Get bump hunter tomography plot
  TGraphErrors bumpHunterTomography = theBumpHunter.GetBumpHunterTomography();

  // The standard operating mode of the BumpHunter is only interested in
  // excesses, since this is how we expect new physics to appear. If you
  // would like to consider both excesses and deficits when finding the
  // region of greatest discrepancy, activate this:
  // theBumpHunter.AllowDeficit(true);

  // The BumpHunter class contains a function you can call like so to
  // compute the test statistic comparing any data histogram to a background estimate:
  // theBumpHunter.DoTest(MjjHistogram & dataHist, MjjHistogram & bkgHist, int firstBinToUse=-1, int lastBinToUse=-1)
  // However, the test statistic is not the important quantity for the BumpHunter.
  // We need to compute a test statistic for each of a range of pseudoexperiments
  // and find the proportion with a more extreme value than that found in our data.

  // Make a pseudoexperimenter. I keep a toy-generating engine
  // as a separate class because I use it for a range of statistical tests
  MjjPseudoExperimenter thePseudinator;

  // Obtain estimate of BumpHunter p-value.
  // The MjjStatisticsBundle is just a struct containing the test statistic in data
  // and the test statistic in all of the histograms (they will be extracted below).
  // The higher the number of pseudoexperiments used, the more precisely known
  // the final p-value will be.
  // Sorry about putting the prediction first and observation second: it's counterintuitive
  // but historical at this point.
  vector<MjjStatisticsBundle> myStats = thePseudinator.GetPseudoExperimentStatsOnHistogram(backgroundHistogram,dataHistogram,theStatsTests,firstBinBH,lastBinBH,nPseudoExpBH);

  // This is the test statistic in data
  //double bumpHunterStat = myStats.originalStatistic;
  MjjStatisticsBundle theseLogLStats = myStats.at(0);
  MjjStatisticsBundle theseChi2Stats = myStats.at(1);
  MjjStatisticsBundle theseBHStats = myStats.at(2);

  std::cout<<"Absolute values of LogL: "<<theseLogLStats.originalStatistic<<std::endl; 
  std::cout<<"Absolute values of Chi2: "<<theseChi2Stats.originalStatistic<<std::endl; 
  std::cout<<"Absolute values of BumpHUnter: "<<theseBHStats.originalStatistic<<std::endl; 
  // This uses the test statistics in pseudoexperiments and the value in data to
  // compute a p-value (with stat uncertainties)
  std::pair<double,double> logLPValAndErr = GetFrequentistPValAndError(theseLogLStats.statisticsFromPseudoexperiments,theseLogLStats.originalStatistic);
  std::pair<double,double> chi2PValAndErr = GetFrequentistPValAndError(theseChi2Stats.statisticsFromPseudoexperiments,theseChi2Stats.originalStatistic);
  std::pair<double,double> bumpHunterPValAndErr = GetFrequentistPValAndError(theseBHStats.statisticsFromPseudoexperiments,theseBHStats.originalStatistic);
  std::cout << "Pval of logL: " << logLPValAndErr.first << std::endl;
  std::cout << "Pval of chi2: " << chi2PValAndErr.first << std::endl;
  std::cout << "Pval of bumpHunter are: " << bumpHunterPValAndErr.first << std::endl;
  
  // The BumpHunter stored the edges of the bins corresponding to the region of
  // greatest excess. If you ran only DoTest on a single spectrum, you
  // could retrieve it like so:
  // vector<double> bumpEdges = theBumpHunter.GetFurtherInformation();
  // double lowEdgeOfBump = bumpEdges.at(0);
  // double highEdgeOfBump = bumpEdges.at(1);
  // However, this contains the information corresponding to the last spectrum
  // examined (in this case, one of the pseudoexperiments).
  // The pseudoexperimenter saved the results in data, so we can access it like this:


  double lowEdgeOfBump = theseBHStats.originalFurtherInformation.at(0);
  double highEdgeOfBump = theseBHStats.originalFurtherInformation.at(1);
  double bumpHunterStat = theseBHStats.originalStatistic ;
  TVectorD bumpHunterStatLowHigh(3);
  bumpHunterStatLowHigh[0] = bumpHunterStat;
  bumpHunterStatLowHigh[1] = lowEdgeOfBump;
  bumpHunterStatLowHigh[2] = highEdgeOfBump;
  
  std::cout << "BumpHunter results: stat = " << theseBHStats.originalStatistic << std::endl;
  
  // Make histo of residual: we'll put this in the plot because it's helpful
  // for visualisation
  MjjSignificanceTests theTestMaker;
  TH1D residualHist = theTestMaker.GetResidual(dataHistogram, backgroundHistogram, firstBinBH, lastBinBH);

  //------------------------------------------
  // Print BH quantities for inspection/checking
  std::cout << "******************************************" << std::endl;
  std::cout << "*** Final values " << std::endl;
  std::cout << "*** BH p-value = " << bumpHunterPValAndErr.first << std::endl;
  std::cout << "*** BH test statistic value = " << theseBHStats.originalStatistic << std::endl;
  std::cout << "*** Selected most discrepant range = " << lowEdgeOfBump << " - " << highEdgeOfBump << std::endl;
  std::cout << "******************************************" << std::endl;
  //------------------------------------------
  infile->Close();

  // save to out file
  TFile* outFile = new TFile(outPath+outFileName, "recreate") ;

  // Save fit range
  TVectorD FitRange(2) ;
  FitRange[0] = rawBkgHisto->GetBinLowEdge(firstBinBH) ;
  FitRange[1] = rawBkgHisto->GetBinLowEdge(lastBinBH)+rawBkgHisto->GetBinWidth(lastBinBH);
  FitRange.Write("FitRange");

  rawDataHisto->SetName("basicData");
  rawDataHisto->Write();
  rawBkgHisto->SetName("basicBkg");
  rawBkgHisto->Write();
  residualHist.SetName("residualHist");
  residualHist.Write();

  theseLogLStats.statisticsFromPseudoexperimentsHist.SetName("logLikelihoodStatHistNullCase");
  theseLogLStats.statisticsFromPseudoexperimentsHist.Write();
  TVectorD statPValErrOfFitToData(3);
  statPValErrOfFitToData[0] = theseLogLStats.originalStatistic;
  statPValErrOfFitToData[1] = logLPValAndErr.first;
  statPValErrOfFitToData[2] = logLPValAndErr.second;
  statPValErrOfFitToData.Write("logLOfFitToData");

  theseChi2Stats.statisticsFromPseudoexperimentsHist.SetName("chi2StatHistNullCase");
  theseChi2Stats.statisticsFromPseudoexperimentsHist.Write();
  statPValErrOfFitToData[0] = theseChi2Stats.originalStatistic;
  statPValErrOfFitToData[1] = chi2PValAndErr.first;
  statPValErrOfFitToData[2] = chi2PValAndErr.second;
  statPValErrOfFitToData.Write("chi2OfFitToData");

  theseBHStats.statisticsFromPseudoexperimentsHist.SetName("bumpHunterStatHistNullCase");
  theseBHStats.statisticsFromPseudoexperimentsHist.Write();
  statPValErrOfFitToData[0] = theseBHStats.originalStatistic;
  statPValErrOfFitToData[1] = bumpHunterPValAndErr.first;
  statPValErrOfFitToData[2] = bumpHunterPValAndErr.second;
  statPValErrOfFitToData.Write("bumpHunterStatOfFitToData");
  bumpHunterTomography.SetName("bumpHunterTomographyFromPseudoexperiments");
  bumpHunterTomography.Write();

  bumpHunterStatLowHigh.Write("bumpHunterPLowHigh");
  outFile->Close() ;

  totaltime.Stop();
  std::cout << "Process ran in " << totaltime.CpuTime() << " seconds. " << std::endl;

  return 0;

}
