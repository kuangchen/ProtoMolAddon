/*  -*- c++ -*-  */
#ifndef ATOM_H
#define ATOM_H

#include <protomol/type/Real.h>
#include <string>
#include <vector>

namespace ProtoMol {
  //________________________________________ Atom

  struct SCPISMAtomParameters {
	  SCPISMAtomParameters() {}

    //Pre force initialization
    void preForce(){
      bornRadius = zeta;
      energySum = true;
      D_s = 0.0;
    }

    Real sqrtalphaSCPISM;

    ///< For SCPISM calculations
    //Updated version CRS 01/24/09

    Real bornRadius, D_s;
    Real zeta, eta;
    bool energySum;
  };

  //add variables to store GBSA paramaters for each atom
  struct GBSAAtomParameters {

    GBSAAtomParameters() {
      bornRadiusDerivatives = NULL;
      Lvalues = NULL;
      Uvalues = NULL;
      distij = NULL;
      btDerv1 = NULL;
      psiDerv1 = NULL;

      fij = NULL;
    }

    ~GBSAAtomParameters() {

      if (bornRadiusDerivatives != NULL) delete [] bornRadiusDerivatives;

      if (Lvalues != NULL) delete [] Lvalues;
      if (Uvalues != NULL) delete [] Uvalues;
      if (distij != NULL) delete [] distij;
      if (btDerv1 != NULL) delete [] btDerv1;
      if (psiDerv1 != NULL) delete [] psiDerv1;

      if (fij != NULL) delete [] fij;

    }

    void SetSpaceForBornRadiusDerivatives(int sz) {
      if (bornRadiusDerivatives == NULL) {
         bornRadiusDerivatives = new Real[sz];
      }

      btDerv1 = new Real[sz];
      psiDerv1 = new Real[sz];

      fij = new Real[sz];
    }

    void SetSpaceLvalues(int sz) {
       Lvalues = new Real[sz];
    }
    void SetSpaceUvalues(int sz) {
       Uvalues = new Real[sz];
    }

    void SetSpaceDistij(int sz) {
       distij = new Real[sz];
    }

    //Pre force initialization
    void preForce() {
      burialTerm = 0.0;
      doneCalculateBornRadius = false;
      doSelfForceTerm = false;
      doneACEPotential = false;
      havePartialGBForceTerms = false;
      havePartialGBHessianTerms = false;
    }

    Real bornRad, burialTerm;
    Real scalingFactor;
    Real offsetRadius;

    //int doneBornRadii;

    Real PsiValue;

    //Need to define an array which will store the derivatives of born radius of atom
    //i w.r.t. all j
    Real *bornRadiusDerivatives;

    Real *Lvalues;
    Real *Uvalues;

    Real *distij;

    Real *btDerv1;

    Real *psiDerv1;

    Real *fij;

    //Defining a flag which will be false if Born Radius value has not been calculated
    //from the burialTerm/PsiValue. It will be set to true first time Born Radius is
    //estimated for THIS atom. It has to be set to false in preForce() routine above.
    bool doneCalculateBornRadius;

    //flag to control calculation of i-i term
    bool doSelfForceTerm;

    // For optimization purposes, we're memoizing the calculation of the forces
    // between i and all of the other atoms.  We use a flag to determine if the
    // force needs to be computed.
    bool havePartialGBForceTerms;
    Real partialGBForceTerms;
    
    //the same for the Hessian calculations
    bool havePartialGBHessianTerms;
    Real partialGBHessianTerms_term1;
    Real partialGBHessianTerms_term2;

    //flag for setting ACE potential
    bool doneACEPotential;

    // Variable for storing Van der Waal radius used in GB calculations.
    // As part of the incremental fix for GB, we realized the Van der Waal radii
    // should be associated with the Atom, not the Atom Type. The Van der Wall radius
    // in AtomType will not be removed until SCPISM has been updated to refer to
    // this field for the Van der Waal Radii.
    Real vanDerWaalRadius;


  };

  /**
   * This class defines the information for one atom.  It contains the type
   * of the atom, its charge (scaled by a constant factor), and the next
   * atom in this atom's cell list.  All other pieces of atom information
   * are either in the atom type or the coordinate vectors (for positions,
   * velocities, and so on).
   */
  struct Atom {
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Atom() : type(-1), scaledCharge(0.0), scaledMass(0.0), hvyAtom(-1),
      atomNum(-1), cellListNext(-1), molecule(-1), mybonds(std::vector<int>()),
             mySCPISM_A(0), myGBSA_T(0) {}

    // A default constructor for the atom class

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    int type;
    ///< The type number of the atom (an index into the array of atomTypes).

    Real scaledCharge;
    ///< The charge of this atom, scaled by the square root of the Coulomb
    ///< constant.

    Real scaledMass;
    ///< The mass of this atom (can be modified, e.g., iSG)

    int hvyAtom;
    ///< The size of the Group for Heavy Atom Come First ordering
    ///< for hydrogen hvyAtom=0, else, hvyAtom = 1+#of attached hydrogens

    int atomNum;
    ///< Original order number of the atom, used to undo Heavy Atom Come
    ///< First ordering

    mutable int cellListNext;
    ///< The index of the next atom in this atom's cell list, or -1 if this
    ///< atom is the last in its list.

    int molecule;
    ///< The ID# of the molecule to which this atom belongs

    std::vector<int> mybonds;
    ///< this vector holds the bond IDs for those bonds made with this atom

    std::string name;
    std::string residue_name;
    int residue_seq;
    SCPISMAtomParameters *mySCPISM_A;

    GBSAAtomParameters *myGBSA_T;

  };
}
#endif /* ATOM_H */

