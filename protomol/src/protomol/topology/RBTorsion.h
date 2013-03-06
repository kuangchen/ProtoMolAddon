/*  -*- c++ -*-  */
#ifndef RBTORSION_H
#define RBTORSION_H

#include <vector>
#include <protomol/type/Real.h>

namespace ProtoMol {
  //________________________________________ RBTorsion

  class RBTorsion {
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    RBTorsion() : atom1(-1), atom2(-1), atom3(-1), atom4(-1),
        C0(0), C1(0), C2(0), C3(0), C4(0), C5(0), Offset(0) {}


    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    int atom1;         ///< The first atom in this interaction.
    int atom2;         ///< The second atom in this interaction.
    int atom3;         ///< The third atom in this interaction.
    int atom4;         ///< The fourth atom in this interaction.

    Real C0; 
    Real C1; 
    Real C2; 
    Real C3; 
    Real C4; 
    Real C5; 
    
    Real Offset;

  };

}
#endif /* RBTORSION_H */
