/*  -*- c++ -*-  */
#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <protomol/type/Real.h>
#include <string>

namespace ProtoMol {
  /**
   * Repository namespace of constants.
   */
  namespace Constant {
    extern const char* PROTOMOL_HR;
    extern const char* PRINTINDENT;
    extern const unsigned int PRINTMAXWIDTH;

    extern const int FASTDELTAMAX;
    extern const Real SQRTCOULOMBCONSTANT;
    extern const Real PRESSUREFACTOR;
    extern const Real BOLTZMANN;
    extern const Real PDBVELSCALINGFACTOR;

    extern const Real MAXREAL;
    extern const Real MINREAL;
    extern const Real REAL_INFINITY;
    extern const Real REAL_NAN;
    extern const int MAX_INT;
    extern const int MAX_INT_2;
    extern const Real EPSILON;
    extern const Real TINY;

    extern const Real TIMEFACTOR;
    extern const Real INV_TIMEFACTOR;
    extern const Real PERIODIC_BOUNDARY_TOLERANCE;

    extern const Real EPS_GOURAUD_THRESHOLD;
    extern const Real EPS_SMOOTH_LINE_FACTOR;

    // Constants for GROMACS-PROTOMOL conversion

    extern const Real ANGSTROM_NM;
    extern const Real NM_ANGSTROM;
    extern const Real INV_ANGSTROM_NM;
    extern const Real INV_NM_ANGSTROM;
    extern const Real KCAL_KJ;
    extern const Real KJ_KCAL;
    extern const Real FS_PS;
    extern const Real PS_FS;

    namespace SI {
      extern const Real C;
      extern const Real COULOMB_FACTOR;
      extern const Real ELECTRON_CHARGE;
      extern const Real LENGTH_AA;
      extern const Real AVOGADRO;
      extern const Real AMU;
      extern const Real KCAL;
      extern const Real TIME_FS;
      extern const Real BOLTZMANN;
    }
    /// Returns actual limits/values of the numerical constants
    std::string print();
  }
}

#endif
