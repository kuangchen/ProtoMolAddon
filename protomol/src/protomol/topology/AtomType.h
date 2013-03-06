/*  -*- c++ -*-  */
#ifndef ATOMTYPE_H
#define ATOMTYPE_H

#include <string>

#include <protomol/type/Real.h>

namespace ProtoMol {
  enum Hbonded {NO = 0, PH = 1, PA = 2};

  struct SCPISMAtomTypeParameters {
    Real sqrt_alpha;        ///< JI: sqrt of alpha i used in CoulombSCPISM
    Real alpha;             ///< JI: alpha i used in CoulombSCPISM
    Real g_i;               ///< TC: used for polar H+, in addendum
    Hbonded isHbonded;      ///< TC: type of H-bond, none, polar, or acceptor
    Real A_i;               ///< JI: A_i in (2) on SCPISM.impl
    Real B_i;               ///< JI: B_i in (2) on SCPISM.impl
    Real C_i;               ///< JI: C_i in (2) on SCPISM.impl    
    //Real eta;               ///< CS: Constant for Born radius
  };

  //________________________________________ AtomType
  /**
   * This class contains information common to one type of atom.  This
   * currently only includes the name and mass of the atom type.  The
   * Lennard-Jones parameters are stored in the LennardJonesParameterTable
   * structure.
   */
  struct AtomType {
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    AtomType() : name(""), mass(0.0), charge(0.0), symbolName(""), mySCPISM_T(0)
    {}

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    std::string name;       ///< The name of this atom type.
    Real mass;              ///< The mass of this atom type.
    Real charge;            ///< The charge of this atom type.
    std::string symbolName; ///< The symbol untity name of this atom type.

    // Van der Waal Radius, used for implicit solvents.
    // Note: In fixing GB, we realized that the Van der Waal radius is
    // a property of the Atom, not the AtomType, and have created a 
    // field in the Atom type appropriately.  As SCPISM still depends
    // on this field, we will not remove it until SCPISM is fixed.
    Real vdwR;              ///< The van del waals radius, used for implicit solvents

    //LJ parameters for openMM and Amber force fields.
    Real sigma, sigma14, epsilon, epsilon14;

    SCPISMAtomTypeParameters *mySCPISM_T;
  };
}
#endif /* ATOMTYPE_H */
