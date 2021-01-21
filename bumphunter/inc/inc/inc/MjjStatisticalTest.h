#ifndef MJJSTATISTICALTEST_H
#define MJJSTATISTICALTEST_H

/*!
 * \class MjjStatisticalTest
 * \brief A skeleton class on which to build various statistical tests
 * \author Katherine Pachal
 * \date 2013
 * 
 * This is the base class for the statistical tests used
 * throughout the package. There are two functions to overload:
 * DoTest returns a double with a number representing some comparison
 * between the input data and background histograms. GetFurtherInformation
 * can be used to pass anything else that the specific test computes.
 */

// ---------------------------------------------------------

#include "inc/MjjHistogram.h"

// ---------------------------------------------------------
class MjjStatisticalTest
{

   public:

      /** \name Constructors and destructors */
      /** @{ */

      /**
       * The default constructor. */
      MjjStatisticalTest();

      /**
       * The default destructor. */
      virtual ~MjjStatisticalTest();

      /** @} */
      /** \name Member functions (set) */
      /** @{ */

      void SetUseWindowExclusion(bool doExclusion)
          { fExcludeWindow = doExclusion; };
 
      void SetWindowToExclude(int lowBin, int highBin)
          { fExcludeWindow = true;
            fFirstBinToExclude = lowBin;
            fLastBinToExclude = highBin; };
 
      /** @} */
      /** \name Member functions (virtual) */
      /** @{ */

      /**
       * Method for obtaining a comparison between data and 
       * background. Pure virtual in this class. In daughter classes,
       * use this to define the statistical test.
       * @param h_data The observed spectrum.
       * @param h_background The predicted spectrum. 
       * @param firstBinToUse The first bin to be included in 
       * the comparison.
       * @param lastBinToUse The last bin to be included in
       * the comparison. 
       * @return The value specified in the daughter class to 
       * represent the comparison between data and prediction. */
      virtual double DoTest(MjjHistogram & h_data, MjjHistogram & h_background, int firstBinToUse, int lastBinToUse) = 0;

      /**
       * Method which allows daughter classes to return 
       * any other useful information from the test: for example,
       * the limits of the bumps found in the BumpHunter.
       * This is NOT pure virtual: if writing a derived class,
       * it does not need to contain this. In that case it just
       * returns and empty vector.
       * @return Vector of further information. */
      virtual vector<double> GetFurtherInformation();

      /** @} */

   protected:

      bool fExcludeWindow;

      int fFirstBinToExclude;

      int fLastBinToExclude;
  
};
// ---------------------------------------------------------

#endif
