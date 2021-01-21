#ifndef MATHFUNCTIONS_H
#define MATHFUNCTIONS_H

/*!
 * \file
 * \brief A library of functions for common math problems
 * \author Kate Pachal
 * \date 2013
 * 
 * This library provides methods for getting
 * getting statistical information out of vectors of 
 * values and out of TGraphs.
 */

// ---------------------------------------------------------

#include <TROOT.h>
#include <TMath.h>
#include <TGraph.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <assert.h>
#include <string>
#include <TH1D.h>

using std::vector;

// ---------------------------------------------------------

      /** 
       * Get mean of a set of values stored as an std::vector<double>.
       * @param input The input vector of values.
       * @return The mean. */
      double GetMean(vector<double> input);

      /** 
       * Get RMS of a set of values stored as an std::vector<double>..
       * @param input The input vector of values. 
       * @return The RMS.*/
      double GetRMS(vector<double> input);

      /** 
       * Get median and 1- and 2-sigma values from a vector
       * of input doubles. Returned values are entries in vector
       * which fall nearest below the specified quantile.
       * @param input The input vector of values. 
       * @return A vector of doubles corresponding to the values
       * at [-2sigma, -1sigma, 0sigma, 1sigma, 2sigma]*/
      vector<double> GetCenterAndSigmaDeviations(vector<double> input);

      /** 
       * Get p-value and error on p-value for an observation given
       * its value in data and a collection of values in pseudo-
       * experiments. The Frequentist p-value is the fraction of
       * cases giving an equal or more extreme value of the observation.
       * @param pseudoExpStatistics A vector of test statistic values in 
       * pseudo-data.
       * @param observedStatistic The test statistc in real data.
       * @return A pair of doubles: Frequenstist p-value and its error. */
      std::pair<double,double> GetFrequentistPValAndError(vector<double> pseudoExpStatistics, double observedStatistic);

      /** 
       * Convert probability as p-value in [0,1] to
       * probability in number of sigmas. Negative return value
       * corresponds to very low significance/high probability.
       * Increasingly positive return values correspond to increasingly
       * interesting results/low probability.
       * @param prob The p-value. 
       * @return The probability in standard deviations.*/
      double ProbToSigma(double prob);

      /** 
       * Reverse of ProbToSigma. Convert probability
       * in number of sigmas to a p-value.
       * @param sigma The probability in standard deviations.
       * @return The p-value. */
      double SigmaToProb(double sigma);

      /** 
       * Calculate Poisson p-value for observed events d and
       * predicted events b. This is probability of observing
       * a value at least as extreme as d given expectation b
       * and Poisson distributed statistics. Equivalent to
       * \f$\sum_{n=D}^{\infty} TMath::PoissonI(d,b) \f$
       * @param d The observed number of events.
       * @param b The expected number of events. 
       * @return The Poisson p-value for d given b. */
      double PoissonPval(const double& d, const double& b);

      /**
       * Calculate a Poisson p-value for observing d events
       * when the predicted value of the mean parameter b is itself uncertain.
       * Assumes that the mean parameter is distributed as a gamma density,
       * with the expectation (mean) of the gamma function = E
       * and the variance of the gamma density equal to err^2.
       * Implementation stolen from Will Buttinger
       * @param d The observed number of events.
       * @param b The expected number of events. 
       * @param bErr The uncertainty on the expected number of events.
       * @return The probability for d marginalised over possible values of b,
       * using a Gamma function prior. */
      double PoissonConvGammaPval(const double& d, const double& b, const double& bErr);

      /**
       * Find the integral of a gaussian with the specified mean and 
       * width between the two specified x values.
       * @param mean The mean of the gaussian
       * @param sigma The width of the gaussian
       * @param x1 The lower bound of the region to integrate
       * @param x2 The upper bound of the region to integrate
       * @return The integral of the specified gaussian over the specified range. */
      double IntegrateGaussian(const double& mean, const double& sigma, const double& x1, const double& x2);

      /**
       * Find the integral under the normal distribution centered 
       * at zero between the two specified x values.
       * @param x1 The lower bound of the region to integrate
       * @param x2 The upper bound of the region to integrate
       * @return The integral of the normal distribution over the specified range. */
      double integrateNormalDistribution(const double& x1, const double& x2);

      /** 
       * Get integral of a TGraph. Space to integrate defined
       * by connecting horizontally nearest points with
       * straight lines (trapezoid approximation of integral).
       * @param graph The input TGraph.
       * @return The integral of the TGraph. */
      double GetIntegralTGraph(TGraph graph);

      /** 
       * Get x-value below which a fraction CL of the
       * contents of a TGraph lie. Extrapolation is linear
       * between the two nearest points in the graph.
       * @param graph The input TGraph.
       * @param CL The quantile of the TGraph to calculate.
       * @return The x value corresponding to
       * a quantile CL of the TGraph. */
      double GetQuantileTGraph(TGraph graph, double CL);

      /**
       * Tests if two doubles are identical within float 
       * precision. Called during generation of effective and
       * weighted histograms to protect against future
       * rounding errors.. */
      bool AreIdentical(double input, double compare);
 
      /** 
       * Transform vector of statistics into a histogram with
       * appropriate axis ranges and granularity. 
       * @param statistics Vector of statistics to be plotted.
       * @return Histogram of statistics. */
      TH1D MakeHistoFromStats(vector<double> statistics);


// ---------------------------------------------------------

#endif
