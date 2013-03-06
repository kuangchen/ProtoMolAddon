/*  -*- c++ -*-  */
#ifndef ANGLE_H
#define ANGLE_H

#include <protomol/type/Real.h>

namespace ProtoMol {
  //________________________________________ Angle

  /// This class contains the information for one Angle.
  class Angle {
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    Angle() : atom1(-1), atom2(-1), atom3(-1), restAngle(0.0),
      forceConstant(0.0), ureyBradleyConstant(0.0), ureyBradleyRestLength(0.0),
      iSGmodifierIndex(0), DeltaK(0.0), DeltaTheta0(0.0), Delta_ubK(0.0),
      Delta_ubR0(0.0) {}

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    int atom1;
    ///< The first atom in this angle (one end).

    int atom2;
    ///< The second atom in this angle (the vertex).

    int atom3;
    ///< The third atom in this angle (the other end).

    Real restAngle;
    ///< The rest angle ($\theta_0$ value) for this angle.

    Real forceConstant;
    ///< The force constant ($k_\theta$ value) for this angle.

    Real ureyBradleyConstant;
    ///< The Urey-Bradley constant ($k_{ub}$ value) for this angle.

    Real ureyBradleyRestLength;
    ///< The Urey-Bradley rest length ($r_{ub}$ value) for this angle.


    int iSGmodifierIndex;
    ///< the index # of this angle type in the ModifierISG angle structure.

    /**
     * difference in spring constant and rest angle for two different
     * angle identities.  This is needed to compute the chemical potential
     * (or free energy) difference between the identities.
     */
    Real DeltaK, DeltaTheta0, Delta_ubK, Delta_ubR0;
  };
}
#endif /* ANGLE_H */
