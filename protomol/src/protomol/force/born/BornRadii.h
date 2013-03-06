/* -*- c++ -*- */
#ifndef BORNRADII_H
#define BORNRADII_H

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

  //____ BornRadii
  class BornRadii {
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Born radii calculation.
  // Equation numbers from "Notes on SCPISM" by C.R.Sweet, based on Hassan's papers.
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  public:
    enum {DIST_R2 = 1};
    enum {CUTOFF = 0};

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Constructors, destructors, assignment
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    BornRadii() : myBornSwitch(3){};
    BornRadii(int bsw) : myBornSwitch(bsw){};

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // New methods of class BornRadii
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    void operator()(Real &energy, Real &force,
                    Real distSquared, Real rDistSquared, const Vector3D &diff,
                    const GenericTopology *topo,
                    int atom1, int atom2, ExclusionClass excl) const {
      
      //SCPISM valid?
      if(!topo->doSCPISM)
        report << error << "BornRadii requires SCPISM parameters." << endr;
      // If either molecule belongs to a water, do nothing.
      // Won't happen in most simulations, but could in the 
      // case of comparing forces.
      if (topo->molecules[topo->atoms[atom1].molecule].water ||
	      topo->molecules[topo->atoms[atom2].molecule].water)
	        return;

      // Distance and shift function
      Real dist = sqrt(distSquared);

      BornSwitch bSw(myBornSwitch);  // Always quartic now  6/1/2008

      Real f_ij = bSw.switchValue(dist);

      // Atom 1 variables
      int type1 = topo->atoms[atom1].type;
      int type2 = topo->atoms[atom2].type;

      //**********************************
      // Part of Eq. (1)
      Real C_i = topo->atomTypes[type1].mySCPISM_T->C_i;
      Real eta_i = topo->atoms[atom1].mySCPISM_A->eta;
      topo->atoms[atom1].mySCPISM_A->bornRadius += eta_i * f_ij * exp(-C_i * dist);
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
        // Eq. (2)
        topo->atoms[atom1].mySCPISM_A->bornRadius += g_i * g_j * f_ij * exp(-E_i * dist);
      }

      //**********************************
      // Part of Eq. (1)      
      Real C_j = topo->atomTypes[type2].mySCPISM_T->C_i;      
      Real eta_j = topo->atoms[atom2].mySCPISM_A->eta;
      topo->atoms[atom2].mySCPISM_A->bornRadius += eta_j * f_ij * exp(-C_j * dist);
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
        //Eq. (2)
        topo->atoms[atom2].mySCPISM_A->bornRadius += g_i * g_j * f_ij * exp(-E_i * dist);

      }
    }

    static void accumulateEnergy(ScalarStructure *energies, Real energy) {
      // Nothing to accumulate here
    }

    static Real getEnergy(const ScalarStructure *energies) {return 0;}

    // Parsing
    static std::string getId() {return keyword;}
    void getParameters(std::vector<Parameter> &) const;
    static unsigned int getParameterSize() {return 1;}

    static BornRadii make(const std::vector<Value> &);

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // My data members
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    static const std::string keyword;

  private:
    int myBornSwitch;
  };

}
#endif /* BORNRADII_H */
