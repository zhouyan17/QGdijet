#ifndef MJJBUMPHUNTER_H
#define MJJBUMPHUNTER_H

/*!
 * \class MjjBumpHunter
 * \brief Implementation of a BumpHunter as defined in arXiv:1101.0390v2
 * \author Katherine Pachal
 * \date 2013
 * 
 * This class creates an object for determining a BumpHunter
 * test statistic to compare a measured spectrum to a theoretical
 * template. For a detailed description, see arXiv:1101.0390v2. 
 * The default BumpHunter uses a window width WC between 1 and NBins/2
 * and sidebands of width max(1,WC). It includes the lowest and highest
 * specified bins in the tests for the bump, ignoring sidebands if none exist
 * outside this range.
 */

// ---------------------------------------------------------


#include <iostream>
#include <math.h>
#include <algorithm>
#include <TMath.h>
#include <Math/SpecFuncMathCore.h>
#include <TRandom3.h>
#include <cmath>

#include "inc/MjjHistogram.h"
#include "inc/MjjStatisticalTest.h"
#include "inc/MjjStatisticsBundle.h"
#include "inc/MathFunctions.h"

// ---------------------------------------------------------
class MjjBumpHunter : public MjjStatisticalTest 
{

   public:

      /** \name Constructors and destructors */
      /** @{ */

      /**
       * The default constructor. */
      MjjBumpHunter();

      /**
       * The default destructor. */
      ~MjjBumpHunter();

      /** @} */
      /** \name Member functions (get) */
      /** @{ */

      /**
       * Retrieves the BumpHunter tomography plot.
       * This TGraph visually represents the probability
       * of a "bump" associated with each range of bins.
       * Width of the shown line (horizontal point error) 
       * marks the relevant bin range; height (vertical 
       * point location) represents probability.
       * @return The TGraphErrors tomography plot. */
      TGraphErrors GetBumpHunterTomography()
          { return fSpectrumTomography; };

      /** @} */

      /** \name Member functions (set) */
      /** @{ */

      /**
       * Determine whether BumpHunter can allow the
       * most significant deviation to be a deficit.
       * @param doAllow Permits deficit if true. */
      void AllowDeficit(bool doAllow) 
          { fAllowDeficit = doAllow; };

      /**
       * Sets minimum width of bump in number of bins.
       * @param nbins Minimum bump width. */
      void SetMaxBumpWidth(int nbins)
          { fMaxBinsInBump = nbins; };

      /**
       * Sets maximum width of bump in number of bins.
       * @param nbins Maximum bump width. */
      void SetMinBumpWidth(int nbins)
          { fMinBinsInBump = std::max(nbins,1); };

      /**
       * Enables sidebands. If activated, some number
       * of bins on either side of the region of interest
       * will also be considered. A bump is ignored if
       * the sidebands also display a discrepancy.
       * @param yesOrNo Activates sidebands if true,
       * deactivates them otherwise. */
      void SetUseSidebands(bool yesOrNo)
          { fUseSidebands = yesOrNo; };
      /**
       * Sets width of sidebands in number of bins.
       * Does NOT automatically activate sidebands. 
       * @param nbins Desired sideband width. */
      void SetSidebandWidth(int nbins)
          { fNBinsInSideband = nbins; };

      /**
       * Sets number of bins by which bump hunter
       * window shifts in scan. Default is 1. In 
       * the case of a very large number of bins,
       * setting this higher will speed up running.
       * @param nbins Desired shift in number of bins */
      void SetNBinsInWindowShift(int nbins)
          { fNBinsInShift = nbins; };

      void SetUseError(TH1D sumInQuad)
          { fErrHist = sumInQuad; fDoErr = true;};

      void SetNoisy(bool isnoisy)
          { fNoisy = isnoisy; };

      /** @} */
      /** \name Member functions (overloaded from MjjStatisticalTest) */
      /** @{ */

      /**
       * Returns the value of the BumpHunter test
       * statistic for the input data and background.
       * This is NOT a frequentist p-value. It is
       * the probability as defined by BumpHunter for the 
       * least probable range of bins in the data spectrum.
       * @param dataHist The observed data spectrum
       * @param bkgHist The expected data spectrum
       * @param firstBinToUse The first bin to be considered.
       * If unspecified, bump hunt will begin at first bin
       * with nonzero data.
       * @param lastBinToUse The last bin to be considered.
       * If unspecified, bump hunt will begin at last bin
       * with nonzero data.
       * @return The probability associated with the least
       * likely interval.   */
      double DoTest(MjjHistogram & dataHist, MjjHistogram & bkgHist, int firstBinToUse=-1, int lastBinToUse=-1);

      /**
       * @return Low edge and high edge of most significant
       * bump, in that order, as vector of size two. */
      vector<double> GetFurtherInformation();

      /** @} */

   private:

      /**
       * Calculate the effective data and background content
       * of a window between firstBin and lastBin inclusive.
       * Although bkgError, dataError, bkgIntegral, and 
       * dataIntegral are parameters, they are not implemented
       * and are a relic of an older BumpHunter.
       * @param histData The effective data histogram
       * @param histBkg The effective background histogram
       * @param firstBin The first bin to be included in the calculation
       * @param lastBin The last bin to be included in the calculation
       * @param dataIntegral Anything. Will be set to 0 anyway.
       * @param dataError Anything. Will be set to 0 anyway.
       * @param bkgIntegral Anything. Will be set to 0 anyway.
       * @param bkgError Anything. Will be set to 0 anyway.  */
      void GetEffectiveBandContentsWithError(TH1D histData, TH1D histBkg, int firstBin, int lastBin, double& dataIntegral, double& dataError, double& bkgIntegral, double& bkgError);

      void DoTestCore(TH1D h_data, TH1D h_bkg, int minWidth, int maxWidth, int firstBin, int lastBin);

      void FindBumpInCaseOfIncalculable(TH1D h_data, TH1D h_bkg, int firstBin, int lastBin);

     /** 
       * Controls whether only excesses are considered
       * or if deficits are also taken into account.
       * Default is false. */
      bool fAllowDeficit;

      /**
       * Minimum width of bump to consider, in number of bins.
       * Default is 1. */
      int fMinBinsInBump;

      /**
       * Maximum width of bump to consider, in number of bins.
       * Default is N/2 where N is total number of bins. */
      int fMaxBinsInBump;

      /**
       * Number of bins by which the search window is shifted in each test.
       * Default is max(1,WC/2) */
      int fNBinsInShift;

      /**
       * Determines whether sidebands are used in search. */
      bool fUseSidebands;

      /**
       * Number of bins included in sidebands (symmetric). Default is max(1,WindowWidth/2). */
      int fNBinsInSideband;

      /**
       * Records whether to treat the searched spectrum as weighted.
       * This affects whether or not errors on data bins are taken into account in 
       * Poisson p-value calculation. */
      bool fIsSpectrumWeighted;

      /**
       * Variables for storing bump information during calculation. */
      double fLowEdge, fHighEdge; //most interesting bump limits
      double fMostInterestingP;

      /**
       * Vectors for storing the low edge, high edge, and probability
       * of each possible bump location. Used in tomography plot. */
      vector<double> fLowEdgesAllBumps, fHighEdgesAllBumps, fProbAllBumps; 

      /**
       * Plot holding location and statistic of bump for each
       * possible window during latest bump search. */
      TGraphErrors fSpectrumTomography;

      bool fDoErr;

      TH1D fErrHist;

      bool fNoisy;

};

#endif
