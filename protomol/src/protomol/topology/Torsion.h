/*  -*- c++ -*-  */
#ifndef TORSION_H
#define TORSION_H

#include <vector>
#include <protomol/type/Real.h>

namespace ProtoMol {
  //________________________________________ Torsion

  /**
   * This class contains the information for one dihedral or improper.
   * Interactions with multiple force terms should be written multiple
   * times, one for each term.
   */
  class Torsion {
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    Torsion() : atom1(-1), atom2(-1), atom3(-1), atom4(-1) {
      forceConstant.clear();
      phaseShift.clear();
      periodicity.clear();

      DeltaK.clear();
      DeltaPhase.clear();
    }


    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    int atom1;         ///< The first atom in this interaction.
    int atom2;         ///< The second atom in this interaction.
    int atom3;         ///< The third atom in this interaction.
    int atom4;         ///< The fourth atom in this interaction.

    int multiplicity;                ///< The number of terms describing the
                                     ///< potential energy of this torsion.
    std::vector<Real> forceConstant; ///< The force constant ($k$ value) for
                                     ///< this interaction.
    std::vector<Real> phaseShift;    ///< The phase shift ($\delta$ value) for
                                     ///< this interaction, in radians.
    std::vector<int> periodicity;    ///< The periodicity ($n$ value) for this
                                     ///< interaction.


    int iSGmodifierIndex;              ///< The index # of this torsion type
                                       ///< in the ModifierISG angle structure.

    /**
     * difference in spring constant and rest angle for two different
     * torsion identities.  This is needed to compute the chemical potential
     *  (or free energy) difference between the identities.
     */
    std::vector<Real> DeltaK, DeltaPhase;
  };
}
#endif /* TORSION_H */
