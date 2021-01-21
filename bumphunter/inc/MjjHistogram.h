#ifndef MJJHISTOGRAM_H
#define MJJHISTOGRAM_H

/*!
 * \class MjjHistogram
 * \brief A class for handling data histograms in various forms
 * \author Katherine Pachal
 * \date 2013
 *
 * This class is mainly a wrapper for a TH1D, but comes
 * with several other parallel histograms for effective statistics
 * and weights. There are also numerous functions for handling
 * the histograms. This is the standard base class for all
 * histogram use in the package.
 */

// ---------------------------------------------------------

#include <iostream>
#include <math.h>
#include <TMath.h>
#include <cmath>

#include "TH1D.h"
#include <TRandom3.h>
#include "TString.h"
#include "TFile.h"
#include <string>
#include <vector>

using std::string;
using std::vector;

// ---------------------------------------------------------
class MjjHistogram : public TObject
{

   public:

      /** \name Constructors and destructors */
      /** @{ */

      /**
       * The standard constructor, taking a TH1D as input.
       * @param inputHisto This will be the basic histogram. */
      MjjHistogram(TH1D * inputHisto, bool getEffFromErrs = true);

      /**
       * A variant on the constructor taking a TH1D as input and
       * a seed with which to initialize the random number generator.
       * @param inputHisto This will be the basic histogram.
       * @param seed The random number generator will be initialized
       * with this seed. */
      MjjHistogram(TH1D * inputHisto, int seed, bool getEffFromErrs = true);


      /**
       * An empty constructor. */
      MjjHistogram() {}

      /**
       * The default destructor. */
      ~MjjHistogram();

      /** @} */
      /** \name Member functions (get) */
      /** @{ */

      /**
       * @return The basic data histogram. */
      TH1D GetHistogram()
          { return fHistogram; };

      /**
       * The effective histogram has a value in each bin
       * equal to the equivalent statistical power of the
       * events in that bin in the basic histogram.
       * @return The effective histogram. */
      TH1D GetEffectiveHistogram()
          { return fEffectiveHistogram; };

      /**
       * The weights histogram has a value in each bin
       * equal to the effective weight of each event in
       * that bin in the effective histogram. It can
       * be used to convert between basic and effective
       * spectra.
       * @return The weights histogram. */
      TH1D GetWeightsHistogram()
          { return fWeightsHistogram; };

      /**
       * This histogram is normalized by bin width
       * compared to the basic data histogram.
       * @return The normalized histogram. */
      TH1D GetNormalizedHistogram()
          { return NormalizeByBinWidth(); };

      /**
       * @return First non-zero bin in basic data histogram. */
      int GetFirstBinWithData()
          { return fFirstBinWithData; };

      /**
       * @return Last non-zero bin in basic data histogram. */
      int GetLastBinWithData(float epsilon = -1)
          { if ((epsilon) < 0) return fLastBinWithData;
            else {
              int saveval = fLastBinWithData;
              GetHistOutermostBinsWithData(epsilon);
              int returnval = fLastBinWithData;
              fLastBinWithData = saveval;
              return returnval;
            }
          };

      /**
       * Finds the smallest set of bins containing the specified
       * fraction of the basic histogram contents.
       * @param percentage The fraction of the histogram to be found.
       * @return pair corresponding to bin numbers of
       * low and high edge (inclusive) of that range. */
      std::pair<int,int> GetIntervalContainingPercentage(double percentage);

      /** @} */
      /** \name Member functions (set) */
      /** @{ */

      /**
       * When an MjjHistogram is created, the effective histogram
       * is automatically generated from the error bars on the
       * bins. In some cases, the user may wish to specify a different
       * set of weights for a histogram, such as when Monte Carlo is
       * being used. Here, an input weight histogram is used to
       * overwrite the automatic one.
       * @param weights The desired weights histogram.  */
      void SetEffectiveFromBasicAndWeights(TH1D * weights);

      /** @} */
      /** \name Member functions (miscellaneous methods) */
      /** @{ */

      /**
       * Function for generating a pseudoexperiment from the MjjHistogram.
       * Pseudoexperiments are generated using the effective
       * histogram and then scaled up with the weights histogram
       * to match the scale of the basic histogram.
       * @return pseudoHist Histogram to be filled with pseudo-data. */
      void PoissonFluctuateBinByBin(TH1D * pseudoHist);

      void PrintSeed(TString desc="")
      {
        std::cout<<"\nSeed ("<<desc<<") : "<<fRandomNumberGenerator.GetSeed()<<std::endl;
      };

      /** @} */

   private:

      /**
       * Called on creation of an MjjHistogram by a TH1.
       * Calculates the range of data and stores it in
       * fFirstBinWithdAta, fLastBinWithData. */
      void GetHistOutermostBinsWithData(float epsilon = -1);

      /**
       * Produces a histogram with each bin value equal
       * to the bin value in the basic histogram divided
       * by the bin width.
       * @return The normalized histogram. */
      TH1D NormalizeByBinWidth();

      /**
       * Called internally when the MjjHistogram is created.
       * Calculates and stores effective and weights histograms
       * from the bin content and bin errors of the basic
       * histogram. */
      void StoreEffectiveHistAndErrors(bool getEffFromErrs);

      /**
       * Tests if two doubles are identical within float
       * precision. Called during generation of effective and
       * weighted histograms to protect against future
       * rounding errors.. */
      bool AreIdentical(double input, double compare);

      /**
       * The basic histogram at the center of the class. */
      TH1D fHistogram;

      /**
       * The effective histogram representing the actual
       * statistical power in each bin of the basic histogram. */
      TH1D fEffectiveHistogram;

      /**
       * The histogram of weights. Each bin contains the effective
       * weight of the events in that bin of the effective histogram. */
      TH1D fWeightsHistogram;

      /**
       * A random number generator for the pseudoexperiments. */
      TRandom3 fRandomNumberGenerator;

      /**
       * The first bin in the basic histogram with nonzero entries. */
      int fFirstBinWithData;

      /**
       * The last bin in the basic histogram with nonzero entries. */
      int fLastBinWithData;

};
// ---------------------------------------------------------

#endif
