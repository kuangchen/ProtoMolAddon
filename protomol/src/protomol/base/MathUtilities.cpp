#include <protomol/base/MathUtilities.h>
#include <protomol/base/Timer.h>

using namespace std;

namespace ProtoMol {
  namespace Rand{
    bool isSeeded = false;
  };

//____ erf
  Real myerf(Real x) {
    // This function returns the error function ERF(x) with fractional
    // error everywhere less than 1.2*10-7.
    // Adapted from the book "Numerical Recipes"

    Real erfcc = 0.0;

    if (x < 9.0) {
      Real z = fabs(x);
      Real t = 1.0 / (1.0 + 0.5 * z);
      erfcc = t * exp(-z * z - 1.26551223 + t *
                      (1.00002368 + t *
                       (.37409196 + t *
                        (.09678418 + t *
                         (-.18628806 + t *
                          (.27886807 + t *
                           (-1.13520398 + t *
                            (1.48851587 + t *
                             (-.82215223 + t * .17087277)))))))));
      if (x < 0.0)
        erfcc = 2.0 - erfcc;
    }

    return 1.0 - erfcc;
  }

//____ poly
  const Real E1 = 0.254829592;      // Polynomial Constants used in
  const Real E2 = -0.284496736;      // Evaluation of the complementary
  const Real E3 = 1.421413741;      // Error function.
  const Real E4 = -1.453152027;      // Approximation used is that of
  const Real E5 = 1.061405429;      // Abramowitz & Stegun p299.
  const Real PP = 0.3275911;

  Real poly5(Real ar) {
    Real t = 1.0 / (1.0 + PP * ar);
    return (t) * (E1 + (t) * (E2 + (t) * (E3 + (t) * (E4 + (t) * E5))));
  }

//____ lowerboundExp2
  Real lowerboundExp2(Real r) {
    Real res = 1.0;
    while (res < r)
      res *= 2.0;

    while (res > r)
      res /= 2.0;

    return res;
  }

//____ Radians <-> Degrees

  vector<Real> rtod(const vector<Real> &rad) {
    vector<Real> deg(rad.size());
    for (unsigned int i = 0; i < rad.size(); ++i)
      deg[i] = rtod(rad[i]);

    return deg;
  }

  vector<Real> dtor(const vector<Real> &deg) {
    vector<Real> rad(deg.size());
    for (unsigned int i = 0; i < deg.size(); ++i)
      rad[i] = dtor(deg[i]);

    return rad;
  }

//____ sincos
//____
//____ Copied from lsys by Jonathan P. Leech.
//____
//____ Computes sine and cosine of angle similar to the
//____ sincos(Real, Real*, Real*) function. Checkes if
//____ the results are close to -1, 0, 1 and rounds
//____ appropriately.
//____

  void sincos(Real alpha, Real &sinAlpha, Real &cosAlpha) {
    sinAlpha = sin(alpha);
    cosAlpha = cos(alpha);

    if (cosAlpha > 1 - Constant::EPSILON) {
      cosAlpha = 1.;
      sinAlpha = 0.;
    } else if (cosAlpha < -1 + Constant::EPSILON) {
      cosAlpha = -1.;
      sinAlpha = 0.;
    }

    if (sinAlpha > 1 - Constant::EPSILON) {
      cosAlpha = 0.;
      sinAlpha = 1.;
    } else if (sinAlpha < -1 + Constant::EPSILON) {
      cosAlpha = 0.;
      sinAlpha = -1.;
    }
  }

//____ power

  Real power(Real x, int n) {
    // Simple and fast way to compute the int-power
    if (n == 0)
      return 1.0;
    if (n < 0)
      return 1.0 / power(x, -n);

    Real z = 1;
    int m = n;
    Real y = x;

    while (m > 0) {
      if (m & 1)
        z *= y;

      y *= y;
      m >>= 1;
    }

    return z;
  }

//____ norm()
  Real norm(Real a, Real b) {
    Real absa = fabs(a);
    Real absb = fabs(b);

    if (absa > absb)
      return absa *sqrt(1.0 + power<2>(absb / absa));
    else
      return absb == 0.0 ? 0.0 : absb *sqrt(1.0 + power<2>(absa / absb));
  }

//____ randomNumber()
  Real randomNumber(unsigned int seed) {
    return randomNumber(seed, 0);
  }

//____ randomNumber()
  Real randomNumber(unsigned int seed, unsigned int randomTypeIn) {
    static int randomType = 0;

    //new random number?
    if (randomTypeIn || randomType) {

      if (!Rand::isSeeded) {
        Rand::isSeeded = true;
        Random::Instance().seed( seed );

        randomType = randomTypeIn;
      }

      return Random::Instance().rand();
    } else {
#ifdef _WIN32
      if (!Rand::isSeeded) {
        srand(seed);

        Rand::isSeeded = true;
      }

      return double (rand()) / double (RAND_MAX);
#else
      if (!Rand::isSeeded) {
        srand48((long)seed);

        Rand::isSeeded = true;
      }

      return drand48();
#endif
    }
  }

//____ randomGaussian()
//____  This section generates a Gaussian random
//____  deviate of 0.0 mean and standard deviation RFD for
//____  each of the three spatial dimensions.
//____  The algorithm is a "sum of uniform deviates algorithm"
//____  which may be found in Abramowitz and Stegun,
//____  "Handbook of Mathematical Functions", pg 952.
  Real randomGaussian(Real sdv, unsigned int seed) {
    Real rnd = 0.0;
    Real sdv2 = 2 * sdv;
    for (int i = 0; i < sdv2; ++i)
      rnd += randomNumber(seed);

    rnd -= sdv;
    return rnd;
  }

//____ randomGaussianNumber()
  Real randomGaussianNumber(unsigned int seed) {
    static bool iset = false;

    static Real gset = 0;

    Real fac = 0.,
         r = 2.,     //  (r >= 1.523e-8) ensures (abs result<6)
         v1 = 0.,
         v2 = 0.;

    if (!iset) { // we do not have an extra result ready,
      while (r >= 1.0 || r < 1.523e-8) { // make sure we are within unit circle
        v1 = 2.0 * randomNumber(seed) - 1.0;

        v2 = 2.0 * randomNumber(seed) - 1.0;

        r = v1 * v1 + v2 * v2;
      }

      fac = sqrt(-2.0 * log(r) / r);

      // now make the Box-Muller transformation to get two normally
      // distributed random numbers. Save one and return the other.
      gset = v1 * fac;

      iset = true;

      return v2 * fac;
    } else { // use previously computed value
      iset = false;

      return gset;
    }
  }

//____ randomGaussianNumber()
  Real randomGaussianNumber(Real mean, Real stdev, unsigned int iseed) {
    //------------------------------------------------------------
    //  FUNCTION: gauss
    //
    //  INPUTS:
    //   mean - average value of gaussian
    //   stdev - standard deviation
    //   iseed - seed for the random number generator
    //
    //  PURPOSE:
    //   The gauss function generates a number from a
    //   near-Gaussian distribution using the Box-Muller method.
    //
    //  REFERENCE:
    //   Box et al., Ann. Math. Stat., v. 29, p.610-611, 1958.
    //------------------------------------------------------------
    // function call flag
    static bool flag = false;

    // independent random variables from the same uniform density function
    static Real zeta1 = 0.0;
    static Real zeta2 = 0.0;

    // random numbers
    Real sqg1 = 0.0;
    Real sqg2 = 0.0;

    // value of gaussian
    Real gaussian = 0.0;

    // if the function call flag is equal to one generate two
    // random numbers, save one and calculate the gaussian else
    // calculate the gaussian
    if (flag) {
      // set function call flag to zero
      flag = false;

      // generate two random numbers on the interval (0,1)
      sqg1 = randomNumber(iseed);
      sqg2 = randomNumber(iseed);

      // generate random variables
      zeta1 = sqrt(-2.0 * log(sqg1)) * cos(2 * M_PI * sqg2);
      zeta2 = sqrt(-2.0 * log(sqg1)) * sin(2 * M_PI * sqg2);

      // calculate gaussian
      gaussian = mean + stdev * zeta1;
    } else {
      // set function call flag to one
      flag = true;

      // calculate gaussian
      gaussian = mean + stdev * zeta2;
    }

    return gaussian;
  }

//____ getTimerSeed()
  int getTimerSeed() {
    Real currentTime = Timer::getCurrentTime().getRealTime();
    return static_cast<int>((currentTime * 10000 -
                             floor(currentTime * 10000)) * 100000);
  }

//____ splitRangeQuadratic()
  void splitRangeQuadratic(unsigned int p, unsigned int from, unsigned int to,
                           vector<PairUInt> &fromRange,
                           vector<PairUInt> &toRange) {
    fromRange.clear();
    toRange.clear();

    const int size = to - from;

    if (p > 2) {
      for (unsigned int i = 0; i < p; i++) {
        fromRange.push_back
          (PairUInt(static_cast<unsigned int>(from + (size * i) / p),
            static_cast<unsigned int>(from + (size * (i + 1.0)) / p)));
        toRange.push_back
          (PairUInt(static_cast<unsigned int>(from + (size * i) / p),
            static_cast<unsigned int>(from + (size * (i + 1.0)) / p)));
      }

      for (unsigned int i = 0; i < p; i++)
        for (unsigned int j = i + 1; j < p; j++) {
          fromRange.push_back
            (PairUInt(static_cast<unsigned int>(from + (size * i) / p),
              static_cast<unsigned int>(from + (size * (i + 1.0)) / p)));
          toRange.push_back
            (PairUInt(static_cast<unsigned int>(from + (size * j) / p),
              static_cast<unsigned int>(from + (size * (j + 1.0)) / p)));
        }

    } else
      for (unsigned int i = 0; i < p; i++)
        for (unsigned int j = i; j < p; j++) {
          fromRange.push_back
            (PairUInt(static_cast<unsigned int>(from + (size * i) / p),
              static_cast<unsigned int>(from + (size * (i + 1.0)) / p)));
          toRange.push_back
            (PairUInt(static_cast<unsigned int>(from + (size * j) / p),
              static_cast<unsigned int>(from + (size * (j + 1.0)) / p)));
        }

  }

//____ splitRangeArea()
  void splitRangeArea(unsigned int p, unsigned int from, unsigned int to,
                      vector<PairUInt> &fromRange,
                      vector<PairUInt> &toRange) {
    fromRange.clear();
    toRange.clear();

    const int size = to - from;
    int first = 0;
    int second = size;

    for (unsigned int i = 0; i < p; i++) {
      second =
        static_cast<int>(size - 0.5 *
                         sqrt(power<2>(1.0 + size * 2.0) -
                              (4.0 * (i + 1.0) * size * (size + 1.0)) / p)
                         + 0.5);
      fromRange.push_back(PairUInt(from + first, from + second));
      toRange.push_back(PairUInt(from + first, to));
      first = second;
    }
  }
}
