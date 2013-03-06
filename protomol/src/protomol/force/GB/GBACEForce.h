/* -*- c++ -*- */
#ifndef GBACEFORCE_H
#define GBACEFORCE_H

#include <protomol/topology/GenericTopology.h>
#include <protomol/type/ScalarStructure.h>
#include <protomol/topology/ExclusionTable.h>
#include <protomol/config/Parameter.h>
#include <protomol/base/MathUtilities.h>

#include <protomol/base/Report.h>

#define PI 3.14169
#define SIXPOW 6

using namespace ProtoMol::Report;

namespace ProtoMol {

  class GBACEForce {

  public:

    enum {DIST_R2 = 1};
    enum {CUTOFF = 0};

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Constructors, destructors, assignment
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   GBACEForce() : sigma(0), rho_s(0) {}
   GBACEForce(Real s, Real r_s) : sigma(s), rho_s(r_s) {}

   //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class GBForce
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    void operator()(Real &energy, Real &force, Real distSquared,
                    Real rDistSquared, const Vector3D &,
                    const GenericTopology *topo, int atom1, int atom2,
                    ExclusionClass excl) const {


      //r_{ij}
      Real dist = sqrt(distSquared);

      Real radius_i = topo->atoms[atom1].myGBSA_T->vanDerWaalRadius; 
      Real radius_j = topo->atoms[atom2].myGBSA_T->vanDerWaalRadius;

     //offset radii ({\tilde{\rho}_{j}})
     Real offsetRadius_i = radius_i - topo->atoms[atom1].myGBSA_T->offsetRadius;
     Real offsetRadius_j = radius_j - topo->atoms[atom2].myGBSA_T->offsetRadius;


     //Scaling factors
     Real S_i = topo->atoms[atom1].myGBSA_T->scalingFactor;
     Real S_j = topo->atoms[atom2].myGBSA_T->scalingFactor;

     Real psi_i, psi_j;

     Real tanhparam_i, tanhparam_j;

     psi_i = topo->atoms[atom1].myGBSA_T->PsiValue;
     psi_j = topo->atoms[atom2].myGBSA_T->PsiValue;
     tanhparam_i = topo->alphaObc*psi_i - topo->betaObc*psi_i*psi_i + topo->gammaObc*psi_i*psi_i*psi_i;
     tanhparam_j = topo->alphaObc*psi_j - topo->betaObc*psi_j*psi_j + topo->gammaObc*psi_j*psi_j*psi_j;


      Real bornRad_i = topo->atoms[atom1].myGBSA_T->bornRad;
      Real bornRad_j = topo->atoms[atom2].myGBSA_T->bornRad;

      Real ratio_i, ratio_j;

      //Check Equation (14) for ACE nonpolar solvation potential
      if (!topo->atoms[atom1].myGBSA_T->doneACEPotential) {
         Real p1 = radius_i/topo->atoms[atom1].myGBSA_T->bornRad;
         ratio_i = power(p1,SIXPOW);
         energy += 4*PI*sigma*(radius_i + rho_s)*(radius_i + rho_s)*ratio_i;
         topo->atoms[atom1].myGBSA_T->doneACEPotential = true;
      }

      if (!topo->atoms[atom2].myGBSA_T->doneACEPotential) {
         Real p2 = radius_j/topo->atoms[atom2].myGBSA_T->bornRad;
         ratio_j = power(p2,SIXPOW);
         energy += 4*PI*sigma*(radius_j + rho_s)*(radius_j + rho_s)*ratio_j;
         topo->atoms[atom2].myGBSA_T->doneACEPotential = true;
      }

      Real c1 = 24*PI*sigma;

      Real c_i = (radius_i + rho_s)*(radius_i + rho_s)*pow(radius_i,SIXPOW);
      Real c_j = (radius_j + rho_s)*(radius_j + rho_s)*pow(radius_j,SIXPOW);


      Real Lij = topo->atoms[atom1].myGBSA_T->Lvalues[atom2];
      Real Uij = topo->atoms[atom1].myGBSA_T->Uvalues[atom2];

      Real Lji = topo->atoms[atom2].myGBSA_T->Lvalues[atom1];
      Real Uji = topo->atoms[atom2].myGBSA_T->Uvalues[atom1];

      //Derivatives for calculation of the derivative of the born radii
      Real dLijdrij, dUijdrij, dCijdrij;

      //Check Equations (11-13)
      if (offsetRadius_i <= (dist - S_j*offsetRadius_j)) dLijdrij = 1;
      else dLijdrij = 0;

      if (offsetRadius_i < (dist + S_j*offsetRadius_j)) dUijdrij = 1;
      else dUijdrij = 0;

      if (offsetRadius_i <= (S_j*offsetRadius_j - dist)) dCijdrij = 2*(1/Lij)*(1/Lij)*dLijdrij;
      else dCijdrij = 0;


     //Derivatives for calculation of the derivative of the born radii
      Real dLjidrij, dUjidrij, dCjidrij;

      if (offsetRadius_j <= (dist - S_i*offsetRadius_i)) dLjidrij = 1;
      else dLjidrij = 0;

      if (offsetRadius_j < (dist + S_i*offsetRadius_i)) dUjidrij = 1;
      else dUjidrij = 0;

      if (offsetRadius_j <= (S_i*offsetRadius_i - dist)) dCjidrij = 2*(1/Lji)*(1/Lji)*dLjidrij;
      else dCjidrij = 0;

      //tanhk angle
      Real tanh_i = tanh(tanhparam_i);

      //tanhi angle
      Real tanh_j = tanh(tanhparam_j);

     //tanhk derivative
     Real tanhparam_derv_i = (topo->alphaObc - 2*topo->betaObc*psi_i + 3*topo->gammaObc*psi_i*psi_i);

     //tanhi derivative
     Real tanhparam_derv_j = (topo->alphaObc - 2*topo->betaObc*psi_j + 3*topo->gammaObc*psi_j*psi_j);

     Real S_i_term = (S_i*offsetRadius_i)/dist;
     Real S_j_term = (S_j*offsetRadius_j)/dist;

     //Check Equation (10) for the derivative of the burial term
     Real dBTidrij = -0.5*dLijdrij*(1/(Lij*Lij)) + 0.5*dUijdrij*(1/(Uij*Uij)) + 0.125*((1/(Uij*Uij)) - (1/(Lij*Lij))) + 0.125*dist*((2/(Lij*Lij*Lij))*dLijdrij - (2/(Uij*Uij*Uij))*dUijdrij) - 0.25*(1/(dist*dist))*log(Lij/Uij) + (Uij/(4*dist*Lij))*((1/Uij)*dLijdrij - (Lij/(Uij*Uij))*dUijdrij) - 0.125*power(S_j_term,2)*((1/(Lij*Lij)) - (1/(Uij*Uij))) + 0.25*((S_j*S_j*offsetRadius_j*offsetRadius_j)/(dist*Uij*Uij*Uij))*dUijdrij - 0.25*((S_j*S_j*offsetRadius_j*offsetRadius_j)/(dist*Lij*Lij*Lij))*dLijdrij + dCijdrij;

     //Check Equation (9)
     Real dRidrij = power(bornRad_i,2)*offsetRadius_i*(1-tanh_i*tanh_i)*tanhparam_derv_i*(1/radius_i)*dBTidrij;

     topo->atoms[atom1].myGBSA_T->btDerv1[atom2] = dBTidrij;
     topo->atoms[atom1].myGBSA_T->bornRadiusDerivatives[atom2] = dRidrij;


     Real dBTjdrji = -0.5*dLjidrij*(1/(Lji*Lji)) + 0.5*dUjidrij*(1/(Uji*Uji)) + 0.125*((1/(Uji*Uji)) - (1/(Lji*Lji))) + 0.125*dist*((2/(Lji*Lji*Lji))*dLjidrij - (2/(Uji*Uji*Uji))*dUjidrij) - 0.25*(1/(dist*dist))*log(Lji/Uji) + (Uji/(4*dist*Lji))*((1/Uji)*dLjidrij - (Lji/(Uji*Uji))*dUjidrij) - 0.125*power(S_i_term,2)*((1/(Lji*Lji)) - (1/(Uji*Uji))) + 0.25*((S_i*S_i*offsetRadius_i*offsetRadius_i)/(dist*Uji*Uji*Uji))*dUjidrij - 0.25*((S_i*S_i*offsetRadius_i*offsetRadius_i)/(dist*Lji*Lji*Lji))*dLjidrij + dCjidrij;

     Real dRjdrji = power(bornRad_j,2)*offsetRadius_j*(1-tanh_j*tanh_j)*tanhparam_derv_j*(1/radius_j)*dBTjdrji;

     topo->atoms[atom2].myGBSA_T->btDerv1[atom1] = dBTjdrji;
     topo->atoms[atom2].myGBSA_T->bornRadiusDerivatives[atom1] = dRjdrji;
     //Check Equation (15)
     force += c1*(c_i*(1/power(bornRad_i,7))*dRidrij*(1/dist) + c_j*(1/power(bornRad_j,7))*dRjdrji*(1/dist));

  }

   static void accumulateEnergy(ScalarStructure *energies, Real energy) {
      (*energies)[ScalarStructure::COULOMB] += energy;
    }

    static Real getEnergy(const ScalarStructure *energies) {
      return (*energies)[ScalarStructure::COULOMB];
    }

    // Parsing
    static std::string getId() {return keyword;}
    static unsigned int getParameterSize() {return 2;}
    void getParameters(std::vector<Parameter> &) const;

    static GBACEForce make(const std::vector<Value> &);


  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // My data members
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    static const std::string keyword;

  private:
    Real sigma;
    Real rho_s;

  };
}

#endif
