/* -*- c++ -*- */
#ifndef BORNSELFFORCE_H
#define BORNSELFFORCE_H

#include <protomol/topology/GenericTopology.h>
#include <protomol/type/ScalarStructure.h>
#include <protomol/topology/ExclusionTable.h>
#include <protomol/config/Parameter.h>
#include <protomol/force/born/BornSwitch.h>
#include <string>

#include <protomol/base/Report.h>

using namespace ProtoMol::Report;

namespace ProtoMol {
  //____ BornSelfForce
  class BornSelfForce {
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Born self energy and force.
    // Equation numbers from "Notes on SCPISM" by C.R.Sweet, based on Hassan's papers.
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    enum {DIST_R2 = 1};
    enum {CUTOFF = 0};

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Constructors, destructors, assignment
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    BornSelfForce() : myBornSwitch(3), myDielecConst(80.0) {};
    BornSelfForce(int bsw, Real dielec) : myBornSwitch(bsw), myDielecConst(dielec) {};

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class BornSelfForce
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    void operator()(Real &energy, Real &force, Real distSquared,
                    Real rDistSquared, const Vector3D &,
                    const GenericTopology *topo, int atom1, int atom2,
                    ExclusionClass excl) const {

      //SCPISM valid?
      if(!topo->doSCPISM)
        report << error << "BornSelfForce requires SCPISM parameters." << endr;
      //Stuff from the old born radii calculations      
      Real dist = sqrt(distSquared);
      Real rDist = 1.0 / dist;
      BornSwitch bSw(myBornSwitch);  // Always quartic now  6/1/2008
      Real f_ij = bSw.switchValue(dist);
      Real fp_ij = bSw.switchDerivative(dist);

      // Atom 1 variables
      int type1 = topo->atoms[atom1].type;      // Atom 2 variables
      int type2 = topo->atoms[atom2].type;

      //**********************************
      // Eq. (5)
      Real C_i = topo->atomTypes[type1].mySCPISM_T->C_i;
      Real eta_i = topo->atoms[atom1].mySCPISM_A->eta;
      Real dRs_i_j = eta_i * (fp_ij - f_ij * C_i) * exp(-C_i * dist);
      //**********************************
      // If atom 1 is a polar H+, accumulate the derivative dR.
      if (topo->atomTypes[type1].mySCPISM_T->isHbonded == PH &&
          topo->atomTypes[type2].mySCPISM_T->isHbonded == PA &&
          (topo->atoms[atom1].residue_seq !=
           topo->atoms[atom2].residue_seq)) { // Polar
        Real E_i = 0.80; // Currently all polar H+ have this value
        Real g_i;
        if (topo->atoms[atom1].name == "HN" &&
            topo->atoms[atom2].name == "O")
          g_i = -0.378;
        else
          g_i = topo->atomTypes[type1].mySCPISM_T->g_i;
        Real g_j = topo->atomTypes[type2].mySCPISM_T->g_i;
        //topo->atoms[atom1].mySCPISM->polarFrac += g_i * g_j * f_ij * exp(
        //  -E_i * dist);
        //**********************************
        // Eq. (7)
        dRs_i_j +=  g_i * g_j *
          (fp_ij - f_ij * E_i) * exp(-E_i * dist);
        //**********************************
      }

      // TC: I am commenting in which equations
      //     in my documentation correspond to which terms.
      Real D = myDielecConst; // bulk solvent dielectric '-D'
      // IMPORTANT: Eq. 1 assumes redundant pairs for the Sasa Fraction.
      // i.e. (atom1-atom2) and (atom2-atom1) would both contribute
      // ProtoMol does not compute redundant pairs.
      // Therefore, we have to multiply by 2.

      // More variable definitions
      Real bornRadius = topo->atoms[atom1].mySCPISM_A->bornRadius;
      Real rBornRadius = 1.0 / bornRadius;
      Real K = (D - 1.0) / 2;
      Real alpha_i = topo->atomTypes[type1].mySCPISM_T->alpha;
      //Screening
      Real Ds_r = topo->atoms[atom1].mySCPISM_A->D_s;
      if(Ds_r == 0.0){
        // Eq. (9)
        Ds_r = (1.0 + D) / (1 + K * exp(-alpha_i * bornRadius)) - 1.0;
        topo->atoms[atom1].mySCPISM_A->D_s = Ds_r;
      }
      Real rDs_r = 1.0 / Ds_r;
      // Eq. (12)
      Real dDs_r = alpha_i * (1.0 / (1.0 + D)) * (1 + Ds_r) * (D - Ds_r);
      Real q_i = topo->atoms[atom1].scaledCharge;

      // I term of Eq. (10)
      force = (-0.5 * q_i * q_i * rBornRadius * rBornRadius *
        (1.0 - rDs_r - bornRadius * rDs_r * rDs_r * dDs_r) * dRs_i_j) * rDist;

      // J term
      //**********************************
      // Eq. (5)
      Real C_j = topo->atomTypes[type2].mySCPISM_T->C_i;
      Real eta_j = topo->atoms[atom2].mySCPISM_A->eta;
      Real dRs_j_i = eta_j * (fp_ij - f_ij * C_j) * exp(-C_j * dist);
      //**********************************
      // If atom 2 is a polar H+, accumulate polar
      // fraction and derivative
      if (topo->atomTypes[type2].mySCPISM_T->isHbonded == PH &&
          topo->atomTypes[type1].mySCPISM_T->isHbonded == PA &&
          (topo->atoms[atom2].residue_seq !=
           topo->atoms[atom1].residue_seq)) {
        Real E_i = 0.80; // Currently all polar H+ have this value
        Real g_i;
        if (topo->atoms[atom2].name == "HN" && topo->atoms[atom1].name == "O")
          g_i = -0.378;
        else g_i = topo->atomTypes[type2].mySCPISM_T->g_i;
        Real g_j = topo->atomTypes[type1].mySCPISM_T->g_i;
        //**********************************
        // Eq. (7)
        dRs_j_i += g_i * g_j *
         (fp_ij - f_ij * E_i) * exp(-E_i * dist);
        //**********************************
      }


      Real bornRadius_j = topo->atoms[atom2].mySCPISM_A->bornRadius;
      Real rBornRadius_j = 1.0 / bornRadius_j;
      Real alpha_j = topo->atomTypes[type2].mySCPISM_T->alpha;
      //Screening
      Real Ds_r_j = topo->atoms[atom2].mySCPISM_A->D_s;
      if(Ds_r_j == 0.0){        
        // Eq. (9)
        Ds_r_j = (1.0 + D) / (1 + K * exp(-alpha_j * bornRadius_j)) - 1.0;
        topo->atoms[atom2].mySCPISM_A->D_s = Ds_r_j;
      }
      //Real Ds_r_j =
      //  (1.0 + D) / (1 + K * exp(-alpha_j * bornRadius_j)) - 1.0;
      Real rDs_r_j = 1.0 / Ds_r_j;

      // Eq. (12)
      Real dDs_r_j =
        alpha_j * (1.0 / (1.0 + D)) * (1 + Ds_r_j) * (D - Ds_r_j);
      Real q_j = topo->atoms[atom2].scaledCharge;

      // Add the J term to the force value, remainder of Eq. (10)
      force += (-0.5 * q_j * q_j * rBornRadius_j * rBornRadius_j *
        (1.0 - rDs_r_j - bornRadius_j * rDs_r_j * rDs_r_j * dDs_r_j) * dRs_j_i) * rDist;

      //**********************************************
      // Eq. (8)
      if(topo->atoms[atom1].mySCPISM_A->energySum){
        energy = 0.5 * q_i * q_i * rBornRadius * (rDs_r - 1);
        topo->atoms[atom1].mySCPISM_A->energySum = false;
      }
      if(topo->atoms[atom2].mySCPISM_A->energySum){
        energy += 0.5 * q_j * q_j * rBornRadius_j * (rDs_r_j - 1);
        topo->atoms[atom2].mySCPISM_A->energySum = false;
      }
      //**********************************************
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

    static BornSelfForce make(const std::vector<Value> &);

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // My data members
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    static const std::string keyword;

  private:
    int myBornSwitch;
    Real myDielecConst;
  };


}
#endif /* BORNSELFFORCE_H */
