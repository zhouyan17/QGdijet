#ifndef MJJCHI2TEST_H
#define MJJCHI2TEST_H

/*!
 * \class MjjChi2Test
 * \brief A class for calculating the Chi2 between two MjjHistograms.
 * \author Kate Pachal
 * \date 2013
 * 
 * This class is an implementation of Pearson's chi-squared
 * test. It is an MjjStatisticalTest daughter class.
 */

// ---------------------------------------------------------

#include <math.h>
#include <iostream>
#include <assert.h>

#include "inc/MjjHistogram.h"
#include "inc/MjjStatisticalTest.h"

// ---------------------------------------------------------
class MjjChi2Test : public MjjStatisticalTest
{

   public:

      /** \name Constructors and destructors */
      /** @{ */

      /**
       * The default constructor. */
      MjjChi2Test();

      /**
       * The default destructor. */
      ~MjjChi2Test();

      /** @} */
      /** \name Member functions (overloaded from MjjStatisticalTest) */
      /** @{ */
 
      /**
       * Returns the value of the Pearson's Chi2 test
       * statistic for the input data and background
       * within the specified range. Uses the effective
       * histograms to ensure weighted data handled correctly.
       * @param dataHist The observed data spectrum
       * @param bkgHist The expected data spectrum
       * @param firstBinToUse The first bin to be considered.
       * If unspecified, bump hunt will begin at first bin
       * with nonzero data.
       * @param lastBinToUse The last bin to be considered.
       * If unspecified, bump hunt will begin at last bin
       * with nonzero data.
       * @return The Chi2 of the comparison between the 
       * two histograms.   */
      double DoTest(MjjHistogram & dataHist, MjjHistogram & bkgHist, int firstBinToUse=-1, int lastBinToUse=-1);
  
   private:
 
 
};

// ---------------------------------------------------------

#endif
