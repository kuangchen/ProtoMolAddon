/* -*- c++ -*- */
#ifndef ENERGYSTRUCTURE_H
#define ENERGYSTRUCTURE_H

#include <protomol/base/Proxy.h>
#include <protomol/type/Vector3D.h>
#include <protomol/base/PMConstants.h>

#include <iostream>
using namespace std;

namespace ProtoMol {
  //____________________________________________________________ ScalarStructure

  /**
   * Container holding energies and all kind of scalar values of interest.
   * The values are kept in array of fixed size and in parallel environment
   * the values are reduced from FIRST to LASTREDUCE.
   */
  class ScalarStructure : public Proxy {
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //  Index type of the array of relevant scalars
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    /// Index of the array of relevant scalars
    enum Index {
      FIRST = 0,       // Only internal purpose
      COULOMB = FIRST,
      LENNARDJONES,
      BOND,
      ANGLE,
      DIHEDRAL,
      IMPROPER,
      OTHER,
      VIRIALXX,
      VIRIALXY,
      VIRIALXZ,
      VIRIALYX,
      VIRIALYY,
      VIRIALYZ,
      VIRIALZX,
      VIRIALZY,
      VIRIALZZ,
      MOLVIRIALXX,
      MOLVIRIALXY,
      MOLVIRIALXZ,
      MOLVIRIALYX,
      MOLVIRIALYY,
      MOLVIRIALYZ,
      MOLVIRIALZX,
      MOLVIRIALZY,
      MOLVIRIALZZ,
      COULOMB_DELTAMU,      ///< needed for iSG and oSG simulations
      LENNARDJONES_DELTAMU, ///< needed for iSG and oSG simulations
      BOND_DELTAMU,         ///< needed for iSG and oSG simulations
      ANGLE_DELTAMU,        ///< needed for iSG and oSG simulations
      DIHEDRAL_DELTAMU,     ///< needed for iSG and oSG simulations
      IMPROPER_DELTAMU,     ///< needed for iSG and oSG simulations
      CQFLUCTUATION,        ///< needed for iSG and oSG simulations
      DELTATIME,            ///< needed for iSG and oSG simulations
      LAMBDA_TEMPERATURE,   ///< needed for iSG and oSG simulations
      INTEGRATOR,
      LASTREDUCE,           ///< Last value to be reduced in parallel
                            ///< environment, only internal purpose
      SHADOW = LASTREDUCE,
      LAST             // Only internal purpose
    };
    // Index of relevant scalares

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //  Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    ScalarStructure();

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class ScalarStructure
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    Real &operator[](Index i) {return myTable[static_cast<int>(i)];}
    Real operator[](Index i) const {return myTable[static_cast<int>(i)];}

    /// Clear all of the energies in this structure.
    void clear();

    ScalarStructure &intoAdd(const ScalarStructure &e);
    ScalarStructure &intoAssign(const ScalarStructure &e);
    ScalarStructure &intoSubtract(const ScalarStructure &e);

    /// The potential energy.
    Real potentialEnergy() const;

    /// The atomic pressure.
    Real pressure(Real volume) const;

    /// The Molecular pressure.
    Real molecularPressure(Real volume) const;

    /// The chemical potential difference
    Real deltaMu() const;

    /// Add viral term for one force pair
    void addVirial(const Vector3D &force12, const Vector3D &diff);

    ///  Add molecular viral term for one force pair.
    void addMolVirial(const Vector3D &force12, const Vector3D &diff);

    ///  Add viral term for one force pair and molecular viral term for one
    ///  force pair.
    void addVirial(const Vector3D &force12,
                   const Vector3D &diff,
                   const Vector3D &comDiff);

    /// test if molecular virial tensor desired
    bool molecularVirial() const {return myDoMolecularVirial;}
    bool molecularVirial(bool doMolecularVirial);
    /// test if virial tensor desired
    bool virial() const {return myDoVirial;}
    bool virial(bool doVirial);
    /// test if trajectory output desired
    bool trajectory() const {return myDoTrajectory;}
    bool trajectory(bool doTrajectory);
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  private:
    Real myTable[LAST - FIRST];
    // Table of all relevant scalars
    bool myDoVirial;
    bool myDoMolecularVirial;

    // flag telling when to write out the coordinate or velocity trajectory for
    // oSGMD (grand canonical or semigrand canonical) simulations.  For these
    // types of simulations it is best to output every 100 transformation
    // attempts, instead of every 100 timesteps -- TIM
    bool myDoTrajectory;
  };

  //___________________________________________________________________ INLINES
  inline Real ScalarStructure::potentialEnergy() const {
    return myTable[COULOMB] +
           myTable[LENNARDJONES] +
           myTable[BOND] +
           myTable[ANGLE] +
           myTable[DIHEDRAL] +
           myTable[IMPROPER] +
           myTable[OTHER];
  }

  inline Real ScalarStructure::deltaMu() const {
    return myTable[COULOMB_DELTAMU] +
           myTable[LENNARDJONES_DELTAMU] +
           myTable[BOND_DELTAMU] +
           myTable[ANGLE_DELTAMU] +
           myTable[DIHEDRAL_DELTAMU] +
           myTable[IMPROPER_DELTAMU];
  }

  inline Real ScalarStructure::pressure(Real volume) const {
    return (myTable[VIRIALXX] + myTable[VIRIALYY] +
            myTable[VIRIALZZ]) / 3.0 / volume * Constant::PRESSUREFACTOR;
  }

  inline Real ScalarStructure::molecularPressure(Real volume) const {
    return (myTable[MOLVIRIALXX] + myTable[MOLVIRIALYY] +
            myTable[MOLVIRIALZZ]) / 3.0 / volume * Constant::PRESSUREFACTOR;
  }
}
#endif /* ENERGYSTRUCTURE_H */
