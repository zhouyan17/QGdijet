// ---------------------------------------------------------

#include "inc/MathFunctions.h"

// ---------------------------------------------------------
double GetMean(vector<double> input) 
{

  int size = input.size();
  if (size==0) {
    std::cout << "ERROR: Can't GetMean of zero-size vector." << std::endl;
    assert(false);
  }
  double sum = 0.0;
  for (int i=0; i<size; ++i) {
    sum += input.at(i);
  }
  return sum / (double)(size);

}

// ---------------------------------------------------------
double GetRMS(vector<double> input) 
{

  int size = input.size();
  if (size==0) {
    std::cout << "ERROR: Can't GetRMS of zero-size vector." << std::endl;
    assert(false);
  }
  double mean=GetMean(input);
  double x2=0;
  for (int i=0; i<size; ++i) {
    x2 += (input.at(i)-mean)*(input.at(i)-mean);
  }
  x2 /= (double)size;
  return sqrt(x2);

}

// ---------------------------------------------------------
vector<double> GetCenterAndSigmaDeviations(vector<double> input) 
{

  std::sort(input.begin(),input.end());
  unsigned int nVals = input.size();
  vector<double> statVals;
  double wantEvents;
  int bestEvent;
  double quantiles [5] = {0.02275,0.1587,0.5,0.8413,0.9772};
  for (int i=0; i<5; i++) {
    wantEvents = nVals*quantiles[i];
    bestEvent = (int) wantEvents;
    statVals.push_back(input.at(bestEvent));
  }

  return statVals;

}

// ---------------------------------------------------------
std::pair<double,double> GetFrequentistPValAndError(vector<double> pseudoExpStatistics, double observedStatistic) 
{

  int trials = pseudoExpStatistics.size();
  int successes=0;
  for (int i=0; i<trials; i++){
    if (pseudoExpStatistics.at(i) > observedStatistic) successes++;
  }
  double pVal = (double)successes/(double)trials;
  double pValErr = sqrt(pVal*(1.-pVal))/sqrt(trials);

  return std::make_pair(pVal,pValErr);

}

// ---------------------------------------------------------
double ProbToSigma(double prob) 
{

  assert (prob >=0 && prob <=1);
  /* p = 0.5 - 0.5*erf(s/sqrt(2.0))
     2p = 1 - erf(s/sqrt(2.0))
     erf(s/sqrt(2.0)) = 1-2p
     s/sqrt(2.0) = erfInverse(1-2p)
     s = sqrt(2.0) * erfInverse(1-2p)
   */
  double valtouse = 1-2*prob;
  // if (pval > 0 && pval < 1)
  if (valtouse>-1 && valtouse<1) return sqrt(2.0)*TMath::ErfInverse(valtouse);
  else if (valtouse==1) return std::numeric_limits<double>::max();
  else return -1*std::numeric_limits<double>::max();

}

// ---------------------------------------------------------
double SigmaToProb(double sigma) 
{
  return 0.5*(1 - TMath::Erf(sigma/sqrt(2.0)));
}

// ---------------------------------------------------------
double PoissonPval(const double& d, const double& b) 
{
  //We want to find the integral of poisson probability, to have d >= observed d, given we expect b.  We pass here the observed d.
  //Be careful, the mathematica Gamma[a,b] is equal to TMath::Gamma(a) * ( 1 - TMath::Gamma(a,b) )
  double answer = 1;

  if (d >= b) { //sum upwards
    // p[d_, b_] := b^d/d!*Exp[-b]
    // Sum[p[n, b], {n, d, Infinity}]
    //  answer : (Gamma[d] - Gamma[d, b])/Gamma[d]
    //      double gammad = TMath::Gamma(d);
    //      double gammadb = TMath::Gamma(d) * (1 - TMath::Gamma(d,b));
    //      answer = (gammad - gammadb)/gammad;
    // after simplification:
    answer = TMath::Gamma(d,b);
  }
  else {
    answer = 1. - TMath::Gamma(d+1,b);
  }
  return answer;
}


// ---------------------------------------------------------
double PoissonConvGammaPval(const double& d, const double& b, const double& bErr) {
   //parameters (a,b) of gamma density can be written in terms of the expectation and variance:
   //Gammma(x, a,b);   - note this is not equal to gamma(x) or gamma(a,b), which are different functions

   double beta = b/(bErr*bErr); // = E/V
   double alpha = b*beta; // = E^2/V

   double pval = 0.;

   // If the error is small enough don't bother with all this nonsense
   if (alpha > 100*d) {
      pval = PoissonPval(d,b);
   } else {
      // Use recursive formula to solve:
      // (nData>nMC): pval = 1 - sum(x=0->nData-1, Integral(y=0->inf, Poisson(x,y)Gamma(y,a,b) dy))
      // (nData<=nMC): pval = sum(x=0->nData, Integral(y=0->inf, Poisson(x,y)Gamma(y,a,b) dy))
      // i.e. integrating out the unknown parameter y
      // Recursive formula: P(n;A,B) = P(n-1;A,B) (A+n-1)/(n*(1+B))
      unsigned int stop = d;
      if (d > b) stop--;
      double sum = 0;

      if(alpha>100) {
         /// NB: must work in log-scale or the math fails!
         double logProb = alpha * log(beta/(1+beta));
         sum = exp(logProb); // P(n=0)
         for (unsigned u=1; u<stop+1; u++) {
            logProb += log((alpha+u-1)/(u*(1+beta)));
            sum += exp(logProb);
         }
       } else {
         // Now, working in log scale similarly messes things up.
         double p0 = pow(beta/(1+beta),alpha);
         double pLast = p0;
         sum = p0;
         for (unsigned k=1; k<stop+1; ++k) {
            double p = pLast * (alpha+k-1) / (k*(1+beta));
            sum += p;pLast = p;
        }
      }
      pval = (d > b) ?  1-sum : sum;
   } 

   return pval;
}

// ---------------------------------------------------------
double IntegrateGaussian(const double& mean, const double& sigma, const double& x1, const double& x2) 
{
  assert (sigma > 0);
  return 0.5*(TMath::Erf((x2-mean)/(sqrt(2.)*sigma)) - TMath::Erf((x1-mean)/(sqrt(2.)*sigma)));
}

// ---------------------------------------------------------
double IntegrateNormalDistribution(const double& x1, const double& x2) 
{
  return IntegrateGaussian(0, 1, x1, x2);
}


// ---------------------------------------------------------
double GetIntegralTGraph(TGraph graph) 
{

  // Sort graph by x-val of points
  graph.Sort();
  // Integrate by trapzoids, point-by-point
  double ptLeftX,ptLeftY;
  double ptRightX,ptRightY;
  int points = graph.GetN();
  double integral = 0;
  double yval, width, area;
  for (int point=0; point<points-1; point++) {
    // Center of slice
    graph.GetPoint(point,ptLeftX,ptLeftY);
    graph.GetPoint(point+1,ptRightX,ptRightY);
    width = ptRightX - ptLeftX;
    yval = (ptLeftY + ptRightY)/2.;
    area = width*yval;
    integral += area;
  }

  return integral;

}

// ---------------------------------------------------------
double GetQuantileTGraph(TGraph graph, double CL) 
{

  // Sorts graph and gives us total area
  double total = GetIntegralTGraph(graph);

  double ptLeftX=0,ptLeftY=0,newPtLeftX,newPtLeftY;
  double ptRightX=0,ptRightY=0,newPtRightX,newPtRightY;
  double width, yval;
  double x,area;
  double upToHere = 0;
  int points = graph.GetN();
  for (int point=0; point<points-1; point++) {
    // Center of slice
    graph.GetPoint(point,newPtLeftX,newPtLeftY);
    graph.GetPoint(point+1,newPtRightX,newPtRightY);
    width = newPtRightX - newPtLeftX;
    yval = (newPtLeftY + newPtRightY)/2.;
    area = width*yval;
    if ((upToHere + area)/total > CL) break;
    else {
      upToHere += area;
      ptLeftX = newPtLeftX;
      ptLeftY = newPtLeftY;
      ptRightX = newPtRightX;
      ptRightY = newPtRightY;
    }
  }

  // Now we have pts to left and right. Find exact xval between them that gives us closest CL.
  double slope = (ptRightY-ptLeftY)/(ptRightX-ptLeftX);
  double onLeft = ptLeftX; double onRight = ptRightX;
  double thisQuantile = 0;
  int count = 0;
  while (true && count<100) {
    x = (onLeft + onRight)/2.;
    //y = slope*(x-ptLeftX) + ptLeftY;
    area = upToHere + (slope*(x-ptLeftX)/2. + ptLeftY)*(x-ptLeftX);
    thisQuantile = area/total;
    // Stop loop when quantile close enough to CL
    if (fabs(thisQuantile - CL) < 0.00001) break;
    // Not close enough? then:
    count++;
    if (thisQuantile > CL) onRight = x;
    else onLeft = x;
  }

  return x;

}

bool AreIdentical(double input, double compare) 
{

  if (fabs(compare - input) < std::numeric_limits<float>::epsilon()) return true;
  return false;

}

// ---------------------------------------------------------
TH1D MakeHistoFromStats(vector<double> statistics)
{

  double min=-1, max=-1;

  unsigned int nentries = statistics.size();
  int nBins = (int) nentries/10.;

  double maxVal = -(std::numeric_limits<double>::max());
  double minVal = std::numeric_limits<double>::max();
  for (unsigned int i=0; i<nentries; i++) { 
    double thisVal = statistics.at(i);
    if (thisVal > maxVal)  maxVal = thisVal;
    if (thisVal < minVal)  minVal = thisVal;
  }
  double range = maxVal - minVal;

  if (min==-1.) min = minVal-0.05*range;
  if (max==-1.) max = maxVal+0.05*range;

  TH1D statPlot("statPlot","",nBins,min,max);
  for (unsigned int i=0; i<nentries; i++) {
    statPlot.Fill(statistics.at(i));
  }

  return statPlot;

}
// ---------------------------------------------------------
