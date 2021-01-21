#ifndef MJJPSEUDOEXPERIMENTER_H
#define MJJPSEUDOEXPERIMENTER_H

/*!
 * \class MjjPseudoExperimenter
 * \brief A class for performing pseudoexperiments and collecting statistics.
 * \author Katherine Pachal
 * \date 2013
 * 
 * This class takes an input spectrum and one or more test
 * statistics, generates pseudoexperiments from the spectrum, and 
 * returns a collection of the test statistic(s) corresponding
 * to them. It is set up to handle weighted data and will take the
 * information from the weights histogram and effective histogram
 * in the MjjHistogram passed to it, so please ensure these are
 * correct.
 */

// ---------------------------------------------------------

#include <iostream>
#include <math.h>
#include <TMath.h>
#include <TRandom3.h>
#include <cmath>
#include <limits>
#include <TGraphErrors.h>

#include "inc/MjjHistogram.h"
#include "inc/MjjStatisticalTest.h"
#include "inc/MjjStatisticsBundle.h"
#include "inc/MathFunctions.h"

// ---------------------------------------------------------
class MjjPseudoExperimenter {

   public:

      /** \name Constructors and destructors */
      /** @{ */

      /**
       * The default constructor. */
      MjjPseudoExperimenter();

      /**
       * Constructor including an MjjBATModel for pseudoexperiment
       * generation which accounts for nuisance parameters. */
      //MjjPseudoExperimenter(MjjBATModel * m, unsigned int seed = 1234);

      /**
       * The default destructor. */
      ~MjjPseudoExperimenter();

      /** @} */
      /** \name Member functions (get) */
      /** @{ */

      /**
       * This performs a designated number of pseudoexperiments
       * taking the input histogram as a theoretical prediction
       * or hypothesis. Each pseudoexperiment is compared to this
       * theoretical prediction using the statistical test passed
       * as a parameter within the specified bin range. This test
       * is distinguished by the fact that the template or theory
       * does not change and each pseudoexperiment is compared to
       * it in exactly the same way that the observation is. This
       * means that the calculation of the fraction of pseudo-
       * experiments more extreme than the data corresponds to a
       * frequentist p-value where the null hypothesis is given
       * by the template histogram.
       * @param templateHist The prediction or hypothesis
       * @param observedHist The observed histogram
       * @param theStatTest The statistical test to be used in the
       * comparison
       * @param firstBinToUse First bin to include in comparison
       * @param lastBinToUse Last bin to include in comparison
       * @param nExperiments Number of pseudoexperiments to use.
       * @return An MjjStatisticsBundle with the test statistic
       * value in data and in the pseudoexperiments, plus any
       * further information from the statistical test. */
      MjjStatisticsBundle GetPseudoExperimentStatsOnHistogram(MjjHistogram & templateHist, 
               MjjHistogram & observedHist, MjjStatisticalTest * theStatTest, 
               int firstBinToUse=1, int lastBinToUse=1e3, int nExperiments=1000);

      /**
       * Equivalent to the version of 
       * GetPseudoExperimentStatsOnHistogram() that uses 
       * a single MjjStatisticalTest. This version simply performs
       * each statistical test in the passed vector on the
       * pseduoexperiments. Note this means the same pseudo-
       * experiments are used in each test. The returned vector
       * of MjjStatisticsBundles correspond to the different
       * statistics tests in the same order in which they were
       * in the passed vector.
       * @param templateHist The prediction or hypothesis
       * @param observedHist The observed histogram
       * @param theStatTests The statistical tests to be used in 
       * the comparison
       * @param firstBinToUse First bin to include in comparison
       * @param lastBinToUse Last bin to include in comparison
       * @param nExperiments Number of pseudoexperiments to use.
       * @return A vector of MjjStatisticsBundles with the test 
       * statistic values in data and in the pseudoexperiments, 
       * plus any further information from the statistical test. 
       * Same ordering as the input vector of tests.  */
      vector<MjjStatisticsBundle> GetPseudoExperimentStatsOnHistogram(MjjHistogram & templateHist, 
               MjjHistogram & observedHist, vector<MjjStatisticalTest*> theStatTests, 
               int firstBinToUse=1, int lastBinToUse=1e3, int nExperiments=1000);

      /**
       * This performs a designated number of pseudoexperiments
       * taking the given function fit to the observed histogram
       * as the theoretical prediction or hypothesis. Each
       * pseudoexperiment is generated from this prediction, but
       * each is individually re-fitted with the function. Thus 
       * the test statistics calculated correspond to the
       * difference between the observed data and its fit (in the
       * original statistic) or to the difference between each 
       * pseudoexperiment and its own fit. This can be informative,
       * but note that because the prediction to which they are 
       * being compared is not the same, the fraction of pseudo-
       * experiments more extreme than the data does not really
       * correspond to a frequentist p-value by strict definition.
       * @param functionToFit The fit function to use in generating
       * the theoretical predictions
       * @param observedHist The observed histogram
       * @param theStatTest The statistical test to be used in the
       * comparison
       * @param firstBinToUse First bin to include in comparison
       * @param lastBinToUse Last bin to include in comparison
       * @param nExperiments Number of pseudoexperiments to use.
       * @return An MjjStatisticsBundle with the test statistic
       * value in data and in the pseudoexperiments, plus any
       * further information from the statistical test. */
      //MjjStatisticsBundle GetPseudoExperimentStats_Refitting(MjjFitter * theFitter, MjjFitFunction * functionToFit, MjjHistogram & observedHist, MjjStatisticalTest * theStatTest, int firstBinToUse=1, int lastBinToUse=1e3, int nExperiments=1000);

      /**
       * Equivalent to the version of 
       * GetPseudoExperimentStats_Refitting() that uses
       * a single MjjStatisticalTest. This version simply performs
       * each statistical test in the passed vector on the
       * pseduoexperiments. Note this means the same pseudo-
       * experiments are used in each test. The returned vector
       * of MjjStatisticsBundles correspond to the different
       * statistics tests in the same order in which they were
       * in the passed vector.
       * @param functionToFit The fit function to use in generating
       * the theoretical predictions
       * @param observedHist The observed histogram
       * @param theStatTests The statistical tests to be used in 
       * the comparison
       * @param firstBinToUse First bin to include in comparison
       * @param lastBinToUse Last bin to include in comparison
       * @param nExperiments Number of pseudoexperiments to use.
       * @return A vector of MjjStatisticsBundles with the test 
       * statistic values in data and in the pseudoexperiments, 
       * plus any further information from the statistical test. 
       * Same ordering as the input vector of tests.  */
      //vector<MjjStatisticsBundle> GetPseudoExperimentStats_Refitting(MjjFitter * theFitter, MjjFitFunction * functionToFit, MjjHistogram & observedHist, vector<MjjStatisticalTest*> theStatTests, int firstBinToUse=1, int lastBinToUse=1e3, int nExperiments=1000);

      /**
       * In the case that GetPseudoExperimentStats_Refitting was
       * used, the parameters of the fit function after the fit
       * will have been stored for each pseudoexperiment.
       * @return A histogram of the values after the fit for
       * each parameter in the fit function. Vector of hists
       * is in same order as function parameters. */
      vector<TH1D> GetParameterHistsFromLatestPseudoexeperiments() 
              { return fParameterHistsFromLatestPseudoexperiments; };

      /**
       * In the case that GetPseudoExperimentStats_Refitting was
       * used, the percentage of fits which succeeded will
       * have been calculated and stored.
       * @return Fraction of fits to pseudoexperiments which 
       * succeeded. */
      double GetFitSuccessRateFromLatestPseudoexperiments() 
              { return fFitSuccessRate; };

      std::vector<TH1D> GetToyHists() { return toyHists; };

      std::vector<TH1D> GetFitHists() { return fitHists; };
  
      std::vector<std::vector<double> > GetToyNPs() { return toyNPs; };

   private:

      /** 
       * Transform vector of statistics into a histogram with
       * appropriate axis ranges and granularity. 
       * @param statistics Vector of statistics to be plotted.
       * @return Histogram of statistics. */
//      TH1D MakeHistoFromStats(vector<double> statistics);

      /**
       * The histograms of fit parameters from the latest
       * set of pseudoexperiments with fits. 
       * Empty vector if none have been run. */
      vector<TH1D> fParameterHistsFromLatestPseudoexperiments;

      TH1D * GetToy(MjjHistogram * templateHist);

      /**
       * Rate of success of fits in most recent call to
       * GetPseudoExperimentStats_Refitting */
      double fFitSuccessRate;

      bool fUseWindowExclusion;
  
      int fFirstBinInWindow;

      int fLastBinInWindow;
  
      //MjjBATModel * model;
  
      //MjjBATAnalysisFacility * analysisFacility;
  
      std::vector<TH1D> toyHists;

      std::vector<std::vector<double> > toyNPs;
  
      std::vector<TH1D> fitHists;

      /**
       * A random number generator for the priors.
       * Should no longer be required in BAT 1.0.0. */
      TRandom3 fRandomNumberGenerator;

};
// ---------------------------------------------------------

#endif
