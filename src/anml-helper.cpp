#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double rcpp_median(std::vector<double> x)
{
  // checks
  if (x.size() == 0) {
    throw std::range_error("median of vector with length 0");
    return -1;
  }

  // calculate median
  const int n = x.size();
  std::sort(x.begin(), x.end());

  if (n == 1) // even
  {
    return x[0];
  }
  else if (n % 2 == 0) // even
  {
    return (x[n / 2] + x[n / 2 - 1]) / 2;
  } else // odd
  {
    return x[n / 2];
  }
}

// [[Rcpp::export]]
double anml_helper(NumericVector x, NumericVector ref, NumericVector mad,
                          int iterations, double tolerance)
{
  const int N = x.size();
  double sf = 1;
  std::vector<double> ratios;

  // loop iterations
  for (int j = 0; j < iterations; j++)
  {
    ratios.clear();
    // loop analytes: calculate ratios
    for (int i = 0; i < N; i++)
    {
      const double logval = log2(sf * x[i]);
      const double logref = log2(ref[i]);
      const double b = 2 * mad[i];
      bool in_range = logval < logref + b && logval > logref - b;
      if (in_range)
      {
        // add to ratios
        if (x[i] == 0) stop("division by zero");
        ratios.push_back(ref[i] / x[i]);
      }
    }

    // return early if within tolerance or no valid ratios
    if (ratios.size() > 0)
    {
      double sf_candiate = rcpp_median(ratios);
      if (abs(sf_candiate - sf) < tolerance)
      {
        return sf_candiate;
      }
      else
      {
        sf = sf_candiate;
      }
    }
    else
    {
      return 1;
    }
  }

  return sf;
}


