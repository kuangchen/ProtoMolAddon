/* -*- c++ -*- */
#ifndef GBFORCE_H
#define GBFORCE_H

#include <protomol/topology/GenericTopology.h>
#include <protomol/type/ScalarStructure.h>
#include <protomol/topology/ExclusionTable.h>
#include <protomol/config/Parameter.h>

#include <protomol/base/Report.h>

using namespace ProtoMol::Report;

namespace ProtoMol {

   class GBForce {

  public:

    enum {DIST_R2 = 1};
    enum {CUTOFF = 0};

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Constructors, destructors, assignment
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   GBForce() : soluteDielec(1.0), solventDielec(80.0) {}
   GBForce(Real solute_d, Real solvent_d) : soluteDielec(solute_d), solventDielec(solvent_d) {}

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class GBForce
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    void operator()(Real &energy, Real &force, Real distSquared,
                    Real rDistSquared, const Vector3D &,
                    const GenericTopology *topo, int atom1, int atom2,
                    ExclusionClass excl) const {

      if(!topo->doGBSAOpenMM)
        report << error << "GBForce requires GB  parameters." << endr;

      //r_{ij}
      Real dist = sqrt(distSquared);

      Real radius_i = topo->atoms[atom1].myGBSA_T->vanDerWaalRadius;
      Real radius_j = topo->atoms[atom2].myGBSA_T->vanDerWaalRadius;

      //offset radii ({\tilde{\rho}_{j}})
      Real offsetRadius_i = radius_i - topo->atoms[atom1].myGBSA_T->offsetRadius;
      Real offsetRadius_j = radius_j - topo->atoms[atom2].myGBSA_T->offsetRadius;

      Real bornRad_i, bornRad_j;
      Real psi_i, psi_j;
      Real tanhparam_i, tanhparam_j;


     psi_i = topo->atoms[atom1].myGBSA_T->PsiValue;
     tanhparam_i = topo->alphaObc*psi_i - topo->betaObc*psi_i*psi_i + topo->gammaObc*psi_i*psi_i*psi_i;
     psi_j = topo->atoms[atom2].myGBSA_T->PsiValue;
     tanhparam_j = topo->alphaObc*psi_j - topo->betaObc*psi_j*psi_j + topo->gammaObc*psi_j*psi_j*psi_j;


      bornRad_i = topo->atoms[atom1].myGBSA_T->bornRad;
      bornRad_j = topo->atoms[atom2].myGBSA_T->bornRad;

      //Equation (17)
      Real expterm = ((dist*dist)/(4*bornRad_i*bornRad_j));
      Real fGB_ij = sqrt(dist*dist + bornRad_i*bornRad_j*exp(-expterm));

      Real scaledCharge_i = topo->atoms[atom1].scaledCharge;
      Real scaledCharge_j = topo->atoms[atom2].scaledCharge;

      //Equation (16)
      energy = -(scaledCharge_i*scaledCharge_j)*(1/fGB_ij)*((1/soluteDielec) - (1/solventDielec));


      //self terms (Equation (18))
      if (!topo->atoms[atom1].myGBSA_T->doSelfForceTerm) {
         energy -= 0.5 * (scaledCharge_i*scaledCharge_i)*(1/bornRad_i)*((1/soluteDielec) - (1/solventDielec));
         topo->atoms[atom1].myGBSA_T->doSelfForceTerm = true;
      }

      if (!topo->atoms[atom2].myGBSA_T->doSelfForceTerm) {
         energy -= 0.5 * (scaledCharge_j*scaledCharge_j)*(1/bornRad_j)*((1/soluteDielec) - (1/solventDielec));
         topo->atoms[atom2].myGBSA_T->doSelfForceTerm = true;
      }




    //Scaling factors
     Real S_i = topo->atoms[atom1].myGBSA_T->scalingFactor;
     Real S_j = topo->atoms[atom2].myGBSA_T->scalingFactor;


      //It would be nice if we can precalculate this terms.
      Real Lij, Uij;

     if (offsetRadius_i >=  dist + S_j*offsetRadius_j) {
        Lij = 1;
        Uij = 1;
     }else {
        Lij =(offsetRadius_i > fabs(dist - S_j*offsetRadius_j)) ? offsetRadius_i : fabs(dist - S_j*offsetRadius_j);
        Uij = dist + S_j*offsetRadius_j;
     }

     Real invLij = 1/Lij;
     //Real invUij = 1/Uij;

      //Derivatives for calculation of the derivative of the born radii
      Real dLijdrij, dUijdrij, dCijdrij;
      if (offsetRadius_i <= (dist - S_j*offsetRadius_j)) dLijdrij = 1;
      else dLijdrij = 0;

      if (offsetRadius_i < (dist + S_j*offsetRadius_j)) dUijdrij = 1;
      else dUijdrij = 0;

      if (offsetRadius_i <= (S_j*offsetRadius_j - dist)) dCijdrij = 2*invLij*invLij*dLijdrij;
      else dCijdrij = 0;

      //Lik
      Real Lji = topo->atoms[atom2].myGBSA_T->Lvalues[atom1];
      //Uki
      Real Uji = topo->atoms[atom2].myGBSA_T->Uvalues[atom1];

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


     Real dBTidrij = -0.5*dLijdrij*(1/(Lij*Lij)) + 0.5*dUijdrij*(1/(Uij*Uij)) + 0.125*((1/(Uij*Uij)) - (1/(Lij*Lij))) + 0.125*dist*((2/(Lij*Lij*Lij))*dLijdrij - (2/(Uij*Uij*Uij))*dUijdrij) - 0.25*(1/(dist*dist))*log(Lij/Uij) + (Uij/(4*dist*Lij))*((1/Uij)*dLijdrij - (Lij/(Uij*Uij))*dUijdrij) - 0.125*power(S_j_term,2)*((1/(Lij*Lij)) - (1/(Uij*Uij))) + 0.25*((S_j*S_j*offsetRadius_j*offsetRadius_j)/(dist*Uij*Uij*Uij))*dUijdrij - 0.25*((S_j*S_j*offsetRadius_j*offsetRadius_j)/(dist*Lij*Lij*Lij))*dLijdrij + dCijdrij;
     Real dRidrij = power(bornRad_i,2)*offsetRadius_i*(1-tanh_i*tanh_i)*tanhparam_derv_i*(1/radius_i)*dBTidrij;

     Real dBTjdrji = -0.5*dLjidrij*(1/(Lji*Lji)) + 0.5*dUjidrij*(1/(Uji*Uji)) + 0.125*((1/(Uji*Uji)) - (1/(Lji*Lji))) + 0.125*dist*((2/(Lji*Lji*Lji))*dLjidrij - (2/(Uji*Uji*Uji))*dUjidrij) - 0.25*(1/(dist*dist))*log(Lji/Uji) + (Uji/(4*dist*Lji))*((1/Uji)*dLjidrij - (Lji/(Uji*Uji))*dUjidrij) - 0.125*power(S_i_term,2)*((1/(Lji*Lji)) - (1/(Uji*Uji))) + 0.25*((S_i*S_i*offsetRadius_i*offsetRadius_i)/(dist*Uji*Uji*Uji))*dUjidrij - 0.25*((S_i*S_i*offsetRadius_i*offsetRadius_i)/(dist*Lji*Lji*Lji))*dLjidrij + dCjidrij;

     Real dRjdrji = power(bornRad_j,2)*offsetRadius_j*(1-tanh_j*tanh_j)*tanhparam_derv_j*(1/radius_j)*dBTjdrji;

     Real expterm_ij = (dist*dist)/(4*bornRad_i*bornRad_j);
     Real fGBij = sqrt(dist*dist + bornRad_i*bornRad_j*exp(-expterm_ij));

     //force due to pairwise i-j term
     //Check Equation (19-20). Next line only finds the pairwise terms.
     force -= scaledCharge_i*scaledCharge_j*(1/(fGBij*fGBij))*0.5*(1/fGBij)*((2*dist - 0.5*exp(-expterm_ij)*dist) + exp(-expterm_ij)*dRidrij*(bornRad_j + (dist*dist)/(4*bornRad_i)) + exp(-expterm_ij)*dRjdrji*(bornRad_i + (dist*dist)/(4*bornRad_j)))*(1/dist);

     // calculation of self (i-i and j-j) terms
     force -= 0.5 * scaledCharge_i * scaledCharge_i * dRidrij / (bornRad_i * bornRad_i) / dist;
       
     force -= 0.5 * scaledCharge_j * scaledCharge_j * dRjdrji / (bornRad_j * bornRad_j) / dist;
     //
            
     //new N^2 handeling of interaction with other atoms
     if (!topo->atoms[atom1].myGBSA_T->havePartialGBForceTerms) {
        topo->atoms[atom1].myGBSA_T->partialGBForceTerms = Force_i_term(topo, atom1);
        topo->atoms[atom1].myGBSA_T->havePartialGBForceTerms = true;
     }
      
     force -= (topo->atoms[atom1].myGBSA_T->partialGBForceTerms - Force_i_j_term(topo, atom1, atom2)) *(dRidrij/dist);
      
     if (!topo->atoms[atom2].myGBSA_T->havePartialGBForceTerms) {
        topo->atoms[atom2].myGBSA_T->partialGBForceTerms = Force_i_term(topo, atom2);
        topo->atoms[atom2].myGBSA_T->havePartialGBForceTerms = true;
     }
      
     force -= (topo->atoms[atom2].myGBSA_T->partialGBForceTerms - Force_i_j_term(topo, atom2, atom1))*(dRjdrji/dist);
     //end
      
     force *= ((1/soluteDielec) - (1/solventDielec));

   }

     //estimate the force term for the sum over k,l where k=i,j and l {\neq} j if 
     //k=i and l {\neq}i if k=j
     Real Force_i_term(const GenericTopology *topo, int atom1) const{
       
       Real scaledCharge_i, scaledCharge_l;
       
       Real bornRad_i, bornRad_l;
       
       Real filGB, expterm;
       
       bornRad_i = topo->atoms[atom1].myGBSA_T->bornRad;
       
       scaledCharge_i = topo->atoms[atom1].scaledCharge;
       
       Real force = 0;
       
       Real ril;
       
       unsigned int sz = topo->atoms.size();
       
       for (unsigned int l = 0; l<sz ; l++) {
         
         if (l == (unsigned)atom1) ril = 0;
         else ril = topo->atoms[atom1].myGBSA_T->distij[l];
         
         bornRad_l = topo->atoms[l].myGBSA_T->bornRad;
         scaledCharge_l = topo->atoms[l].scaledCharge;
         
         expterm = (ril*ril)/(4.0*bornRad_i*bornRad_l);
         filGB = sqrt(ril*ril + bornRad_i*bornRad_l*exp(-expterm));
         
         if (l != (unsigned)atom1) {
           //Equation (24)
           force += scaledCharge_i*scaledCharge_l*(1/(filGB*filGB))*0.5*(1/filGB)*exp(-expterm)*(bornRad_l + (ril*ril)/(4.0*bornRad_i));
         }
       }
       
       return force; 
       
     }
     
     
     //estimate the force term for the sum over k,l where k=i,j and l {\neq} j if 
     //k=i and l {\neq}i if k=j
     // modified to calculate the force only between two atoms.  will be used in
     // optimization
     Real Force_i_j_term(const GenericTopology *topo, int atom1, int atom2) const{
       
       Real scaledCharge_i, scaledCharge_l;
       
       Real bornRad_i, bornRad_l;
       
       Real filGB, expterm;
       
       bornRad_i = topo->atoms[atom1].myGBSA_T->bornRad;
       
       scaledCharge_i = topo->atoms[atom1].scaledCharge;
       
       Real force = 0;
       
       Real ril;
       
       unsigned int l = (unsigned) atom2;

       ril = topo->atoms[atom1].myGBSA_T->distij[l];
       
       bornRad_l = topo->atoms[l].myGBSA_T->bornRad;
       scaledCharge_l = topo->atoms[l].scaledCharge;
       
       expterm = (ril*ril)/(4.0*bornRad_i*bornRad_l);
       filGB = sqrt(ril*ril + bornRad_i*bornRad_l*exp(-expterm));
       
       force += scaledCharge_i*scaledCharge_l*(1/(filGB*filGB))*0.5*(1/filGB)*exp(-expterm)*(bornRad_l + (ril*ril)/(4.0*bornRad_i));
       
       return force; 
       
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

    static GBForce make(const std::vector<Value> &);

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // My data members
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    static const std::string keyword;

  private:
    Real soluteDielec;
    Real solventDielec;


   };

}

#endif
