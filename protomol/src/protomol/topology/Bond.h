/*  -*- c++ -*-  */
#ifndef BOND_H
#define BOND_H

#include <protomol/type/Real.h>

namespace ProtoMol {
  //________________________________________ Bond
  /**
     This class contains the information for one bond between two atoms.
   */
  class Bond {
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Typedef & const
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    struct Constraint {
      Constraint() : atom1(-1), atom2(-1), restLength(0.0) {}
      Constraint(int a, int b, Real c) : atom1(a), atom2(b), restLength(c) {}
      int atom1;
      int atom2;
      Real restLength;
    };
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    Bond() : atom1(-1), atom2(-1), springConstant(0.0), restLength(0.0),
      iSGmodifierIndex(0), DeltaK(0.0), DeltaR0(0.0) {}
    Bond(int a, int b, Real c,
         Real d) : atom1(a), atom2(b), springConstant(c), restLength(d) {}
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    int atom1;
    ///< The first atom involved in this bond.

    int atom2;
    /// <The second atom involved in this bond.

    Real springConstant;
    ///< The spring constant ($k$ value) for this bond.

    Real restLength;
    ///< The rest length ($r_0$ value) for this bond.

    int iSGmodifierIndex;
    // /<the index # of this bond type in the ModifierISG bond structure

    /**
     *  difference in spring constant and rest length for two different
     * bond identities.  This is needed to compute the chemical potential
     * (or free energy) difference between the identities.
     */
    Real DeltaK, DeltaR0;
  };
}
#endif /* BOND_H */
