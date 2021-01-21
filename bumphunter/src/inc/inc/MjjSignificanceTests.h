#ifndef MJJSIGNIFICANCETESTS_H
#define MJJSIGNIFICANCETESTS_H

/*!
 * \class MjjSignificanceTests
 * \brief A class for assessing the significance of differences
 * between histograms.
 * \author Katherine Pachal
 * \date 2013
 * 
 * This class has three constituent methods for determining the 
 * significance of the difference between two histograms. All three 
 * methods return a TH1D whose bins have the same spacing and structure
 * as the input histograms. Each bin of the output represents the 
 * significance of the difference between the values of that bin in the inputs.
 */

// ---------------------------------------------------------

#include <iostream>
#include <math.h>
#include <TMath.h>
#include <TH1D.h>
#include <limits>

#include "inc/MjjHistogram.h"
#include "inc/MjjStatisticalTest.h"
#include "inc/MathFunctions.h"

// ---------------------------------------------------------
class MjjSignificanceTests 
{

   public:

      /** \name Constructors and destructors */
      /** @{ */

      /**
       * The default constructor. */
      MjjSignificanceTests();

      /**
       * The default destructor. */
      ~MjjSignificanceTests();
  
      /** @} */
      /** \name Member functions (get) */
      /** @{ */

      /**
       * Expresses significance of difference as (data - background)/background
       * with errors of |(dataError)/bkg|. Uses effective data entries.
       * @param dataHist The observed MjjHistogram
       * @param bkgHist The predicted or theoretical MjjHistogram.
       * @return Histogram with relative difference in each bin. */
      TH1D GetRelativeDifference(MjjHistogram & dataHist, MjjHistogram & bkgHist, int firstBinToUse=-1, int lastBinToUse = -1);

      /**
       * Similar to GetRelativeDifference() except error on background is also
       * used in finding bin errors. Bin error is
       * sqrt((Derr/Berr)^2+((D*Berr)/(B^2))^2)
       * @param dataHist The observed MjjHistogram
       * @param bkgHist The predicted or theoretical MjjHistogram.
       * @return Histogram with relative difference in each bin. */
      TH1D GetRelativeDifferenceWithErrorOnBkg(MjjHistogram & dataHist, MjjHistogram & bkgHist, int firstBinToUse = -1, int lastBinToUse = -1);

      /**
       * Expresses significance of difference as 
       * (data - bkg)/sqrt(dataErr^2 + bkgErr^2) with no bin errors.
       * Uses effective data entries.
       * @param dataHist The observed MjjHistogram
       * @param bkgHist The predicted or theoretical MjjHistogram.
       * @return Histogram with significance of difference in each bin. */
      TH1D GetSignificanceOfDifference(MjjHistogram & dataHist, MjjHistogram & bkgHist, int firstBinToUse=-1, int lastBinToUse = -1);

      /**
       * Expresses significance of difference as a residual as defined 
       * in arXiv:1111.2062. In each bin, p-value found as the Poisson
       * probability for getting a value at least as extreme, and this
       * is converted into a z-value in number of sigmas. For values
       * below zero the result is considered extremely probable and 
       * bin content is plotted empty. For values above zero, bin
       * content is set to the number of sigmas and its sign is given
       * by whether data or fit is higher. Uses effective data entries.
       * @param dataHist The observed MjjHistogram
       * @param bkgHist The predicted or theoretical MjjHistogram.
       * @return Histogram with residual in each bin. */
      TH1D GetResidual(MjjHistogram & dataHist, MjjHistogram & bkgHist, int firstBinToUse=-1, int lastBinToUse=-1, TH1D * errHist = NULL);

      /** @} */
  
   private:

  
};
// ---------------------------------------------------------

#endif
