/* -*- c++ -*- */
#ifndef GBBORNBURIALTERM_H
#define GBBORNBURIALTERM_H

#include <protomol/topology/GenericTopology.h>
#include <protomol/type/ScalarStructure.h>
#include <protomol/topology/ExclusionTable.h>
#include <protomol/type/Real.h>
#include <protomol/type/Vector3D.h>
#include <protomol/config/Parameter.h>
#include <protomol/force/born/BornSwitch.h>
#include <string>

#include <protomol/base/Report.h>

using namespace ProtoMol::Report;

namespace ProtoMol {

  class GBBornBurialTerm {

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // GB Born radii calculation.
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


  public:
    enum {DIST_R2 = 1};
    enum {CUTOFF = 0};

   public:
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Constructors, destructors, assignment
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    GBBornBurialTerm(){}

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // New methods of class GBBornBurialTerm
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:

    void operator()(Real &energy, Real &force,
                    Real distSquared, Real rDistSquared, const Vector3D &diff,
                    const GenericTopology *topo,
                    int atom1, int atom2, ExclusionClass excl) const {

    Real one = (Real) 1.0 ;
    Real two = (Real) 2.0 ;
    Real four = (Real) 4.0;

     if (!topo->doGBSAOpenMM) {
       report << error <<"GBBornBurialTerm require GBSA set"<<endr;
     }

     energy = 0;

     force = 0;



      // If either molecule belongs to a water, do nothing.
      // Won't happen in most simulations, but could in the 
      // case of comparing forces.
      if (topo->molecules[topo->atoms[atom1].molecule].water ||
              topo->molecules[topo->atoms[atom2].molecule].water)
                return;

     //r_ij
     Real dist = sqrt(distSquared);


      //store pairwise distances in some data structure for future use
      topo->atoms[atom1].myGBSA_T->distij[atom2] = dist;
      topo->atoms[atom2].myGBSA_T->distij[atom1] = dist;

      Real radius_i = topo->atoms[atom1].myGBSA_T->vanDerWaalRadius;
      Real radius_j = topo->atoms[atom2].myGBSA_T->vanDerWaalRadius;

     //offset radii ({\tilde{\rho}_{j}}) Equation (4) in the writeup
     Real offsetRadius_i = radius_i - topo->atoms[atom1].myGBSA_T->offsetRadius;
     Real offsetRadius_j = radius_j - topo->atoms[atom2].myGBSA_T->offsetRadius;

     //Scaling factors
     Real S_i = topo->atoms[atom1].myGBSA_T->scalingFactor;
     Real S_j = topo->atoms[atom2].myGBSA_T->scalingFactor;


     Real Lij, Uij, Cij;

     //Equations (6-8)
     if (offsetRadius_i >=  dist + S_j*offsetRadius_j) {
        Lij = one;
        Uij = one;

     }else {
        Lij = (offsetRadius_i > fabs(dist - S_j*offsetRadius_j)) ? offsetRadius_i : fabs(dist - S_j*offsetRadius_j);

        Uij = dist + S_j*offsetRadius_j;
        
     }
     if (offsetRadius_i < offsetRadius_j*S_j- dist) {
       Cij = two*(one/offsetRadius_i - one/Lij);
     }else {
       Cij = 0;
     }

     //store Lvalues Uvalues in an array so that we can use later
     topo->atoms[atom1].myGBSA_T->Lvalues[atom2] = Lij;
     topo->atoms[atom1].myGBSA_T->Uvalues[atom2] = Uij;


     Real invLij = one/Lij;
     Real invUij = one/Uij;

     Real invLij2 = invLij*invLij;
     Real invUij2 = invUij*invUij;

     Real ratio_i = log(Lij/Uij);


     //add this to the burial term for atom i (atom1) (see Equation (5) in the writeup)
     Real term_i = (invLij - invUij) + (dist/four)*(invUij2 - invLij2) + (one/(two*dist))*(Real)ratio_i + ((S_j*S_j*offsetRadius_j*offsetRadius_j)/(four*dist))*(invLij2 - invUij2) + Cij;
     topo->atoms[atom1].myGBSA_T->burialTerm += term_i;


    // These are required for jith term. Note that Lij not necessary equal to Lji. Same is true for Uij and Cij.
    // These quantities may not be symmetric. I need these for calculating R_j.
    Real Lji, Uji, Cji;

    if (offsetRadius_j >=  dist + S_i*offsetRadius_i) {
        Lji = one;
        Uji = one;
     }else {
        Lji = (offsetRadius_j > abs(dist - S_i*offsetRadius_i)) ? offsetRadius_j : abs(dist - S_i*offsetRadius_i);

        Uji = dist + S_i*offsetRadius_i;
     }

     //Store calculated L and U values
     topo->atoms[atom2].myGBSA_T->Lvalues[atom1] = Lji;
     topo->atoms[atom2].myGBSA_T->Uvalues[atom1] = Uji;

     if (offsetRadius_j < offsetRadius_i*S_i - dist) {
       Cji = two*(one/offsetRadius_j - one/Lji);

     }else {
       Cji = 0;
     }

     Real invLji = one/Lji;
     Real invUji = one/Uji;

     Real invLji2 = invLji*invLji;
     Real invUji2 = invUji*invUji;
 
     Real ratio_j = log(Lji/Uji);

     //add this to the burial term for atom j (atom2) (see Equation (5) in the writeup)
     Real term_j = (invLji - invUji) + (dist/four)*(invUji2 - invLji2) + (one/(two*dist))*(Real)ratio_j + ((S_i*S_i*offsetRadius_i*offsetRadius_i)/(four*dist))*(invLji2 - invUji2) + Cji;
     topo->atoms[atom2].myGBSA_T->burialTerm += term_j;

    
        
   }

   static void accumulateEnergy(ScalarStructure *energies, Real energy) {
      //(*energies)[ScalarStructure::COULOMB] += energy;
    }

    static Real getEnergy(const ScalarStructure *energies) {
      return (*energies)[ScalarStructure::COULOMB];
    }

    // Parsing
    static std::string getId() {return keyword;}
    static unsigned int getParameterSize() {return 0;}
    void getParameters(std::vector<Parameter> &) const {}

    static GBBornBurialTerm make(const std::vector<Value> &) {
      return GBBornBurialTerm();
    }
 

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // My data members
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    static const std::string keyword;

 

  };

}

#endif /* GBBornBurialTerm */
