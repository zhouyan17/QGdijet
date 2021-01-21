#ifndef MJJLOGLIKELIHOODTEST_H
#define MJJLOGLIKELIHOODTEST_H

/*!
 * \class MjjLogLikelihoodTest
 * \brief A definition of the log likelihood for use anywhere
 * \author Katherine Pachal
 * \date 2013
 * 
 * This class implements a log likelihood which can
 * either do a simple calculation for one bin
 * or return a value for two MjjHistograms. When given 
 * MjjHistograms it handles the weights appropriately
 * using the effective statistical information stored in
 * the inputs.
 */

// ---------------------------------------------------------

#include <math.h>
#include <iostream>
#include <assert.h>
#include <float.h>

#include "inc/MjjHistogram.h"
#include "inc/MjjStatisticalTest.h"

#include <TMath.h>

// ---------------------------------------------------------
class MjjLogLikelihoodTest : public MjjStatisticalTest
{

   public:

      /** \name Constructors and destructors */
      /** @{ */

      /**
       * The default constructor. */
      MjjLogLikelihoodTest();

      /**
       * The default destructor. */
      ~MjjLogLikelihoodTest();

      /** @} */
      /** \name Member functions (overloaded from MjjStatisticalTest) */
      /** @{ */

      /**
       * Returns the value of the log likelihood 
       * for the input data and background
       * within the specified range. Uses the effective
       * histograms to ensure weighted data handled correctly.
       * @param dataHist The observed data spectrum
       * @param backgroundHist The expected data spectrum
       * @param firstBinToUse The first bin to be considered.
       * If unspecified, bump hunt will begin at first bin
       * with nonzero data.
       * @param lastBinToUse The last bin to be considered.
       * If unspecified, bump hunt will begin at last bin
       * with nonzero data.
       * @return The log likelihood for the comparison between 
       * two histograms.   */
      double DoTest(MjjHistogram & dataHist, MjjHistogram & backgroundHist, int firstBinToUse=-1, int lastBinToUse=-1);

      /**
       * Returns the value of the log likelihood 
       * for the input single bin values of data
       * and background. Please remember if using this
       * to handle weights appropriately yourself.
       * @param d Observed value in bin.
       * @param b Predicted value in bin.
       * @return Log of poisson likelihood of observation
       * given prediction.  */
      double LogLikelihood(double d, double b);  
  
   private:
 
};

#endif
