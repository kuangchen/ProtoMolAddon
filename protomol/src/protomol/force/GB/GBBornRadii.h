/* -*- c++ -*- */
#ifndef GBBORNRADII_H
#define GBBORNRADII_H


#include <protomol/topology/GenericTopology.h>
#include <protomol/type/ScalarStructure.h>
#include <protomol/topology/ExclusionTable.h>
#include <protomol/config/Parameter.h>
#include <protomol/base/MathUtilities.h>

#include <protomol/base/Report.h>

using namespace ProtoMol::Report;

namespace ProtoMol {

  class GBBornRadii {

  public:

    enum {DIST_R2 = 1};
    enum {CUTOFF = 0};

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Constructors, destructors, assignment
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    GBBornRadii(){}

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // New methods of class GBBornRadii
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  public:

    void operator()(Real &energy, Real &force,
                    Real distSquared, Real rDistSquared, const Vector3D &diff,
                    const GenericTopology *topo,
                    int atom1, int atom2, ExclusionClass excl) const {

      //r_{ij}
      //Real dist = sqrt(distSquared);


      Real radius_i = topo->atoms[atom1].myGBSA_T->vanDerWaalRadius;
      Real radius_j = topo->atoms[atom2].myGBSA_T->vanDerWaalRadius;

     //offset radii ({\tilde{\rho}_{j}})
     Real offsetRadius_i = radius_i - topo->atoms[atom1].myGBSA_T->offsetRadius;
     Real offsetRadius_j = radius_j - topo->atoms[atom2].myGBSA_T->offsetRadius;


     //Scaling factors
     //Real S_i = topo->atoms[atom1].myGBSA_T->scalingFactor;
     //Real S_j = topo->atoms[atom2].myGBSA_T->scalingFactor;

     Real psi_i, psi_j;

     Real tanhparam_i, tanhparam_j;

      if (!topo->atoms[atom1].myGBSA_T->doneCalculateBornRadius)  {
        //calculate born radius
        //set flag to true

        //Equation (3) and (5)
        topo->atoms[atom1].myGBSA_T->PsiValue = 0.5*topo->atoms[atom1].myGBSA_T->burialTerm*offsetRadius_i;
        psi_i = topo->atoms[atom1].myGBSA_T->PsiValue;

        //part of Equation (1)
        tanhparam_i = topo->alphaObc*psi_i - topo->betaObc*psi_i*psi_i + topo->gammaObc*psi_i*psi_i*psi_i;
        //Second part of Equation (1)
        Real invBornRad_i = (1/offsetRadius_i) - (1/radius_i)*tanh(tanhparam_i);
        topo->atoms[atom1].myGBSA_T->bornRad = 1/invBornRad_i;
        topo->atoms[atom1].myGBSA_T->doneCalculateBornRadius = true;

      }

      if (!topo->atoms[atom2].myGBSA_T->doneCalculateBornRadius)  {
        //calculate born radius
        //set flag to true

         //Equation (3) and (5)
         topo->atoms[atom2].myGBSA_T->PsiValue = 0.5*topo->atoms[atom2].myGBSA_T->burialTerm*offsetRadius_j;
        psi_j = topo->atoms[atom2].myGBSA_T->PsiValue;
 
        //part of Equation (1)
        tanhparam_j = topo->alphaObc*psi_j - topo->betaObc*psi_j*psi_j + topo->gammaObc*psi_j*psi_j*psi_j;
        //Second part of Equation (1)
         Real invBornRad_j = (1/offsetRadius_j) - (1/radius_j)*tanh(tanhparam_j);
         topo->atoms[atom2].myGBSA_T->bornRad = 1/invBornRad_j;
         topo->atoms[atom2].myGBSA_T->doneCalculateBornRadius = true;

      }


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

    static GBBornRadii make(const std::vector<Value> &) {
      return GBBornRadii();
    }


  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // My data members
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    static const std::string keyword;


     

  };

}

#endif
