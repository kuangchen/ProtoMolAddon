/*  -*- c++ -*-  */
#ifndef MATHUTILITIES_H
#define MATHUTILITIES_H

#include <protomol/base/PMConstants.h>
#include <protomol/type/SimpleTypes.h>
#include <protomol/base/Random.h>

#include <cmath>
#include <algorithm>
#include <queue>
#include <vector>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#ifdef WIN32
#include <float.h>
#endif

#ifdef HAVE_NO_LIMITS
#ifndef WIN32
#include <values.h>
#endif
#else /* HAVE_NO_LIMITS */
#include <limits>
#endif /* HAVE_NO_LIMITS */
namespace ProtoMol {

#ifdef _MSC_VER
  //___________________________________________________________________ isnan()

  inline int isnan(Real x) {
    bool b1 = (x < 0.0);
    bool b2 = (x >= 0.0);
    return !(b1 || b2);
  }
#endif

  //___________________________________________________________________ square()
  inline Real square(Real x) {
    return x * x;
  }

  //___________________________________________________________________ power()
#if defined (NO_GLOBAL_TEMPLATE_SPECIALIZATION)
  namespace Private {
    template<int _power>
    struct Power {
      static Real power(Real x) {
        return square(Power < _power / 2 >
            ::power(x)) *Power < _power - 2 *(_power / 2) > ::power(x);
      }
    };

    template<>
    struct Power<1> {
      static Real power(Real x) {return x;}
    };

    template<>
    struct Power < -1 > {
      static Real power(Real x) {return 1 / x;}
    };

    template<>
    struct Power<0> {
      static Real power(Real) {return 1;}
    };
  }

  template<int _power> inline Real power(Real x) {
    return Private::Power<_power>::power(x);
  }


#else
  template<int _power> inline Real power(Real x) {
    return square(power < _power / 2 >
      (x)) * power < _power - 2 * (_power / 2) > (x);
  }

  template<> inline Real power<0>(Real) {
    return 1;
  }
  template<> inline Real power<1>(Real x) {
    return x;
  }
  template<> inline Real power < -1 > (Real x){return 1. / x;}
#endif

  //___________________________________________________________________ power()
  Real power(Real x, int n);

  //_____________________________________________________________________ fact()
  inline int fact(int n) {
    int a = 1;
    for (int i = 2; i <= n; ++i)
      a *= i;

    return a;
  }

  //__________________________________________________________ BinaryExponent
  namespace Private {
    template<bool cmp, int A> struct BinaryExponentHelper {
      enum {N = 1};
    };

    template<int A> struct BinaryExponentHelper<false, A> {
      enum {N = 1 + BinaryExponentHelper < A <= 1, A / 2 > ::N};
    };
  }

  /**
   * BinaryExponent<int>::N gives the number of required bits
   */
  template<int P> struct BinaryExponent {
    enum {N = Private::BinaryExponentHelper < P <= 1, P / 2 > ::N};
  };

  //_____________________________________________________________ lowerboundExp2
  Real lowerboundExp2(Real r);
  //_______________________________________________________ Radians <-> Degrees
  inline Real dtor(Real degree) {
    return degree * M_PI / 180.0;
  }

  //________________________________________________________ Radians <-> Degrees
  inline Real rtod(Real rad) {
    return rad * 180.0 / M_PI;
  }
  //__________________________________________________________ Radians <-> Degre
  std::vector<Real> rtod(const std::vector<Real> &rad);
  //________________________________________________________ Radians <-> Degrees
  std::vector<Real> dtor(const std::vector<Real> &deg);

  //___________________________________________________________________ sincos()
  /**
   * Computes sine and cosine of angle similar to the
   * sincos(Real, Real*, Real*) function. Checkes if
   * the results are close to -1, 0, 1 and rounds
   * appropriately.
   */
  void sincos(Real alpha, Real &sinAlpha, Real &cosAlpha);

  //_____________________________________________________________________ erf()
  /**
   * Polynomial approximation of the error function with fractional error
   * everywhere less than 1.2*10-7.
   */
  Real myerf(Real x);
#ifdef _WIN32
  inline Real erf(Real x) {
    return myerf(x);
  }
  inline Real erfc(Real x) {
    return 1.0 - myerf(x);
  }
#endif

  //_____________________________________________________________________ poly()
  /// Polynomial constants used in evaluation of the complementary error
  /// function
  Real poly5(Real ar);

  //___________________________________________________________________ equal()
  inline bool equal(Real x, Real y) {
    return fabs(x - y) < Constant::EPSILON;
  }

  //____________________________________________________________________ equal()
  inline bool equal(Real x, Real y, Real epsilon) {
    return fabs(x - y) < epsilon;
  }

  //_____________________________________________________________________ sign()
  template<class T>
  inline int sign(T a) {
    if (a < 0)
      return -1;
    else if (a > 0)
      return 1;
    else
      return 0;
  }

  //_____________________________________________________________________ sign()
  inline Real sign(Real a, Real b) {
    return b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a);
  }


  //_____________________________________________________________________ norm()
  /// 2-norm
  Real norm(Real a, Real b);


  //_____________________________________________________________ randomNumber()
  Real randomNumber(unsigned int seed = 1234);
  Real randomNumber(unsigned int seed, unsigned int randomtype);
  //___________________________________________________________ randomGaussian()
  Real randomGaussian(Real sdv, unsigned int seed = 1234);
  //_____________________________________________________ randomGaussianNumber()
  Real randomGaussianNumber(unsigned int seed = 1234);
  //_____________________________________________________ randomGaussianNumber()
  /// Generates a number from a near-Gaussian distribution using the Box-Muller
  ///  method
  Real randomGaussianNumber(Real mean, Real stdev, unsigned int seed = 1234);
  //_____________________________________________________________ getTimerSeed()
  /// Computes a seed based on the actual time
  int getTimerSeed();

  //______________________________________________________ splitRangeQuadratic()
  /**
   * Splits indices of a diagonal quadratic (sub) matrix (from,from)-(to,to)
   * into p(p-1)/2 quadratic sub matrices.
   */
  void splitRangeQuadratic(unsigned int p,
                           unsigned int from,
                           unsigned int to,
                           std::vector<PairUInt> &fromRange,
                           std::vector<PairUInt> &toRange);
  //___________________________________________________________ splitRangeArea()
  /**
   * Splits indices of a diagonal quadratic (sub) matrix (from,from)-(to,to)
   * into p sub matrices of equal area
   */
  void splitRangeArea(unsigned int p,
                      unsigned int from,
                      unsigned int to,
                      std::vector<PairUInt> &fromRange,
                      std::vector<PairUInt> &toRange);

  //_____________________________________________________________________ max()
  /// return the larger of the numbers a and b
  using std::max;

  //_____________________________________________________________________ max()
  /// return the largest of the numbers a, b, and c
  template<class T>
  inline T max(const T &a, const T &b, const T &c) {
    return (a > b) ? ((a > c) ? a : c) : ((b > c) ? b : c);
  }

  //_____________________________________________________________________ max()
  /// return the largest of the numbers a, b, c, and d
  template<class T>
  inline T max(const T &a, const T &b, const T &c, const T &d) {
    return (a > b) ? ((a > c) ? (a > d ? a : d) : (c > d ? c : d)) :
      ((b > c) ?  (b > d ? b : d) : (c > d ? c : d));
  }

  //_____________________________________________________________________ min()
  /// return the smallest of the numbers a and b
  using std::min;

  //_____________________________________________________________________ min()
  /// return the smallest of the numbers a, b, and c
  template<class T>
  inline T min(const T &a, const T &b, const T &c) {
    return (a < b) ? ((a < c) ? a : c) : ((b < c) ? b : c);
  }

  //_____________________________________________________________________ min()
  /// return the smallest of the numbers a, b, c, and d
  template<class T>
  inline T min(const T &a, const T &b, const T &c, const T &d) {
    return (a < b) ? ((a < c) ? (a < d ? a : d) : (c < d ? c : d)) :
      ((b < c) ?  (b < d ? b : d) : (c < d ? c : d));
  }

  //  ________________________________________________________  Class WindowAve
  //  I wrote this class to be able to easily keep track of the average of   //
  //  a moving window of values.  By default, I save the most 25 recent      //
  //  values and when needed, can quickly compute the average.  Not sure     //
  //  this merits a new class, but I thought it was a good idea.  :)  --ssh  //

  class WindowAve {
private:
    std::queue<Real> data;
    unsigned int maxItems;
    Real sumItems;

public:
    WindowAve(unsigned int max = 25) : maxItems(max), sumItems(0.0) {}
    void clear();
    void add(Real newValue);
    Real getAve();
  };

  inline void WindowAve::clear() {
    while (!data.empty()) data.pop();
  }

  inline void WindowAve::add(Real newValue) {
    Real oldValue = 0.0;

    //  If we have the maximum allowed items, then subtract out the oldest
    //  value.
    if (data.size() == maxItems) {
      oldValue = data.front();
      data.pop();
    }

    data.push(newValue);
    sumItems += (newValue - oldValue);
  }

  inline Real WindowAve::getAve() {
    return data.size() > 0 ? sumItems /
           data.size() : Constant::REAL_NAN;
  }
}
#endif /* MATHUTILITIES_H*/
