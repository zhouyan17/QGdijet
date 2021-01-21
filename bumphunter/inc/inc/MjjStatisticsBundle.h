#ifndef MJJSTATISTICSBUNDLE_H
#define MJJSTATISTICSBUNDLE_H

/*!
 * \struct MjjStatisticsBundle
 * \brief A structure for holding information from statistical tests.
 * \author Katherine Pachal
 * \date 2013
 * 
 * This struct is for passing all the interesting output
 * from the MjjStatisticalTest class. For a given test, this 
 * contains the original test statistic between data and prediction,
 * the further information (if any) from that comparison, the
 * test statistic and further information from each pseudoexperiment,
 * and a histogram showing all the test statistics from the pseudo-
 * experiments together.
 */

// ---------------------------------------------------------

#include <iostream>
#include <string>
#include <sstream>
#include <TH1D.h>
#include <vector>
#include <TGraphErrors.h>

// ---------------------------------------------------------
struct MjjStatisticsBundle {

      /**
       * A histogram of the pseudoexperiment test statistics. */
      TH1D statisticsFromPseudoexperimentsHist;

      /**
       * The test statistic in the observed histogram. */
      double originalStatistic;

      /**
       * The test statistic in each pseudoexperiment. */
      vector<double> statisticsFromPseudoexperiments;

      /**
       * The further information in the observed histogram. */
      vector<double> originalFurtherInformation;

      /**
       * The further information in the pseudoexperiments. */
      vector<vector<double> > furtherInformationFromPseudoexperiments;

};
// ---------------------------------------------------------

#endif
