/*  -*- c++ -*-  */
#ifndef REAL_H
#define REAL_H

// The standard Real number type, settable by the user.
//
// Define double as default type
#if !defined (USE_REAL_IS_FLOAT) && !defined (USE_REAL_IS_DOUBLE)
#define USE_REAL_IS_DOUBLE
#endif

// The following constants are defined for compilation in Windows
#ifdef _MSC_VER
#define rint(X) floor(X + 0.5)
#define BITSPERBYTE 8
#define BITS(type)  (BITSPERBYTE * (long)sizeof(type))
#define HIBITI      (1U << (BITS(int) - 1))
#define MAXINT          ((int)(~HIBITI))
#define MAXDOUBLE   1.7976931348623157e+308
#define MINDOUBLE       2.2250738585072014e-308
#define M_2_SQRTPI     1.12837916709551257390  /* 2/sqrt(pi) */
#define M_SQRT2        1.41421356237309504880  /* sqrt(2) */
#define M_SQRT1_2      0.70710678118654752440  /* 1/sqrt(2) */


#define M_PI    3.14159265358979323846
#define M_PI_2  1.57079632679489661923  // pi/2
#define M_2PI   6.283185307179586476926

#endif

namespace ProtoMol {
  /// The standard Real number type, settable by the user: double
#ifdef USE_REAL_IS_DOUBLE
  typedef double Real;
#endif

  /// The standard Real number type, settable by the user: float
#ifdef USE_REAL_IS_FLOAT
  typedef float Real;
#endif

  inline Real RealAbs(Real x) {return x < 0 ? -x : x;}
}
#endif /*REAL_H*/

