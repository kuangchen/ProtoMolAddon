/* -*- c++ -*- */
#ifndef GENERICTOPOLOGY_H
#define GENERICTOPOLOGY_H

#include <protomol/topology/Atom.h>
#include <protomol/topology/AtomType.h>
#include <protomol/topology/Angle.h>
#include <protomol/topology/Bond.h>
#include <protomol/topology/Torsion.h>
#include <protomol/topology/ExclusionTable.h>
#include <protomol/topology/ExclusionType.h>
#include <protomol/type/Vector3DBlock.h>
#include <protomol/base/Makeable.h>
#include <protomol/topology/Molecule.h>
#include <protomol/base/Zap.h>

#include <protomol/topology/RBTorsion.h>

#include <protomol/topology/LennardJonesParameterTable.h>
#include <protomol/topology/BankLennardJonesParameterTable.h>

namespace ProtoMol {
  //________________________________________ GenericTopology

  /**
   * GenericTopology represents the topology of the system and contains all 
   * force field parameters, which are needed to compute the interactions.
   * GenericTopology is the base class of SemiGenericTopology and Topology,
   * hiding all implementation details from the applications and integrators.
   * GenericTopology itself does not know how the measure distances and how to
   * implement the cell list algorithm. SemiGenericTopology is in charge of
   * implementing the boundary conditions and all other services related to
   * distances. Topology implements the cell algorithm, which relates to
   * boundary conditions.
   * Topology<class BC, class CM> is the actual implementation used by forces,
   * since they are instantiated with the same template parameters they
   * can correctly access the needed services provided by Topology. Some forces
   * will get along with only accessing SemiGenericTopology or even 
   * GenericTopology. @n
   *
   * There are a couple of service depending on Topology<class BC, class CM>,
   * which are passed through such one can access them, but one may consider
   * the price of calling virtual methods. @n
   *
   * The design was motived by these facts:
   * - we do not what to mess with templates at application and integrator level
   * - we do not what to pay for the cost of calling virtual methods in
   *   the most inner loops, i.e., forces and potential evaluation
   * - measurement of distances must be non-virtual
   */

  ///Available force fields
  enum ForceFieldType {
    CHARMM,
    GROMACS
  };

  ///Available implicit solvents
  enum ImplicitSolventType {
    NONE,
    SCPISM,
    GBSA
  };
    
  class GenericTopology : public Makeable<GenericTopology> {
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    GenericTopology();
    GenericTopology(Real c, const ExclusionType &e);
    virtual ~GenericTopology() {
      if (doSCPISM) {
        for (unsigned int i = 0; i < atoms.size(); i++)
          zap(atoms[i].mySCPISM_A);

        for (unsigned int i = 0; i < atomTypes.size(); i++)
          zap(atomTypes[i].mySCPISM_T);
      }
    };

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class GenericTopology
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    /// Actual volume of the simulation box or the particles
    virtual Real getVolume(const Vector3DBlock &positions) const = 0;
    /// Volume by boundary conditions
    virtual Real getVolume() const = 0;
    /// rescales the volume of the system cell/minimal image
    virtual void rescaleVolume(Real fac) = 0;
    /// Tags the cell list as out of date
    virtual void uncacheCellList() = 0;
    /// Sets default parameters basaed on particle positions for cell manager
    /// and bounary conditions
    virtual std::vector<Parameter> getDefaults(const Vector3DBlock &positions)
    const = 0;
    /// computes the minimal bounding box of the positions taking into account
    /// of the boundary conditions
    virtual void getBoundingbox(const Vector3DBlock &positions,
                                Vector3D &min,
                                Vector3D &max) const = 0;
    /// returns the bounding box of the boundary conditions
    virtual void getBoundaryConditionsBox(Vector3D &min,
                                          Vector3D &max) const = 0;

    /// Perform a minimal-image subtraction.
    virtual Vector3D minimalDifference(const Vector3D &c1,
                                       const Vector3D &c2) const = 0;
    /// Find the position in the basis/original cell/image.
    virtual Vector3D minimalPosition(const Vector3D &c) const = 0;

    /**
     * checks whether the plain distances between all atoms on each molecule are
     * the same as with applying the boundary conditions. If true, there is no
     * molecule, which is wrapped around.
     */
    virtual bool checkMoleculePairDistances(const Vector3DBlock &positions)
    const = 0;
    /**
     * translates all positions of a each molecule such that the plain 
     * distances between them yields and that the center of mass is inside the 
     * minimal image/basis cell. It may not work for molecules which are bigger
     * then half of the simulation cell.
     */
    virtual void minimalImage(Vector3DBlock &positions) = 0;

    virtual std::string print(const Vector3DBlock *positions = NULL) const = 0;

    GenericTopology *make(const std::vector<Value> &values) const;

    static const std::string &getKeyword() {return keyword;}

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class Makable
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    virtual std::string getScope() const {return scope;}
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    std::vector<Atom>     atoms;
    std::vector<AtomType> atomTypes;
    std::vector<Bond>     bonds;
    std::vector<Angle>    angles;
    std::vector<Torsion>  dihedrals;
    std::vector<Torsion>  impropers;
    std::vector<Molecule> molecules;

    //Ryckert-Belleman potential
    std::vector<RBTorsion> rb_dihedrals;

    /// look-up table of exclusions
    ExclusionTable exclusions;
    /// type of exclusion on intra-molecular interactions
    ExclusionType exclude;
    /// Lennard-Jones table for all atoms
    LennardJonesParameterTable lennardJonesParameters;
    /// coulomb scaling factor for modified 1-4 exclusions
    Real coulombScalingFactor;

    ///LJ scaling factor for Gromacs/AMBER force fields
    Real LJScalingFactor;

    /// the actual time of the system, intended to by modified/incremented by
    /// an Observer (Modfier*)
    mutable Real time;
    mutable Vector3D min;
    mutable Vector3D max;

    static const std::string scope;
    static const std::string keyword;

    /// a bank of LennardJonesParameterTables needed only for iSGMD
    BankLennardJonesParameterTable isgLJParms;

    /// the # of molecules of each component
    std::vector<int> iSGNumMols;

    //Implicit Solvent flags
    ImplicitSolventType implicitSolvent;
    // Flag for SCPISM
    int doSCPISM;

    //flag for force field
    ForceFieldType forceFieldFlag;

    /**
     * the number of degrees of freedom in the system
     * this number is needed to make sure we are properly computing the 
     * temperature the number of degrees of freedom is = 3 * (# of atoms) - 
     * (# of constraints)
     */
    int degreesOfFreedom;

    /// list of bond constraints (RATTLE/SHAKE)
    std::vector<Bond::Constraint> bondRattleShakeConstraints;

    /// if the distances on the same molecule are minimal regardless boundary 
    /// conditions
    bool minimalMolecularDistances;

    //Parameters for GBSA with OpenMM
    int doGBSAOpenMM;
    int obcType;
    Real dielecOffset;
    Real alphaObc, betaObc, gammaObc;

  };
  //________________________________________ INLINES
}
#endif /* GENERICTOPOLOGY_H */
