/* -*- c++ -*- */
#ifndef LENNARDJONESPARAMETERS_H
#define LENNARDJONESPARAMETERS_H

#include <protomol/type/Real.h>

namespace ProtoMol {
  //________________________________________ LennardJonesParameters

  /// The Lennard-Jones parameters for one pair of atom types.
  struct LennardJonesParameters {
  public:
    LennardJonesParameters() : A(0.0), B(0.0), A14(0.0), B14(0.0) {};
    LennardJonesParameters(Real a, Real b) : A(a), B(b), A14(0.0), B14(0.0) {};
    LennardJonesParameters(Real a, Real b, Real a14,
                           Real b14) : A(a), B(b), A14(a14), B14(b14) {};

    Real A, B;
    ///< The non-1-4 parameters for these atom types.

    Real A14, B14;
    ///< The 1-4 parameters for these atom types.
  };
}
#endif /* not LENNARDJONESPARAMETERS_H */
