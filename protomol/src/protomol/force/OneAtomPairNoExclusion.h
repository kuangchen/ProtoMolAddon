/* -*- c++ -*- */
#ifndef ONEATOMPAIRNOEXCLUSION_H
#define ONEATOMPAIRNOEXCLUSION_H

#include <protomol/topology/Topology.h>
#include <protomol/config/Parameter.h>
#include <protomol/type/ScalarStructure.h>
#include <protomol/force/OneAtomContraints.h>

namespace ProtoMol {
  //____ OneAtomPairNoExclusion

  /**
   * Computes the interaction for a given force between two atoms with the
   * template arguments defining the boundary conditions, switching function,
   * potential and optional constraint.
   */
  template<class TBoundaryConditions, class TSwitchingFunction,
           class TNonbondedForce, class TConstraint = NoConstraint>
  class OneAtomPairNoExclusion {
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Typedef & sub classes
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    typedef TBoundaryConditions BoundaryConditions;
    typedef SemiGenericTopology<TBoundaryConditions> TopologyType;
    // Make the boundary conditions visible

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    OneAtomPairNoExclusion() : switchingFunction(), nonbondedForceFunction() {};
    OneAtomPairNoExclusion(TNonbondedForce nF, TSwitchingFunction sF) :
      switchingFunction(sF), nonbondedForceFunction(nF),
      mySquaredCutoff(Cutoff<TNonbondedForce::CUTOFF>::cutoff(sF, nF)) {}

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class OneAtomPairNoExclusion
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    void initialize(const TopologyType *topo, const Vector3DBlock *pos,
                    Vector3DBlock *f, ScalarStructure *e) {
      realTopo = (TopologyType *)topo;
      positions = pos;
      forces = f;
      energies = e;
    }

    void initialize(TopologyType *topo, const Vector3DBlock *pos,
                    Vector3DBlock *f, ScalarStructure *e) {
      initialize(static_cast<const TopologyType *>(topo), pos, f, e);
    }

    // Computes the force and energy for atom i and j.
    void doOneAtomPair(const int i, const int j) {
      if (TConstraint::PRE_CHECK)
        if (!TConstraint::check(realTopo, i, j))
          return;

      // Get atom distance.
      Real distSquared;

      Vector3D diff = realTopo->
        boundaryConditions.minimalDifference((*positions)[i], (*positions)[j],
                                             distSquared);
      //      cout << "DIFF: " << diff << endl;
      // Do switching function rough test, if necessary.
      if (TSwitchingFunction::USE || TNonbondedForce::CUTOFF)
        if (distSquared > mySquaredCutoff)
          return;
      // Don't Check for an exclusion.
      int mi = realTopo->atoms[i].molecule;
      int mj = realTopo->atoms[j].molecule;
      bool same = (mi == mj);
      ExclusionClass excl =
        (same ? realTopo->exclusions.check(i, j) : EXCLUSION_NONE);
      //if (excl == EXCLUSION_FULL)
      //  return;

      // Calculate the force and energy.
      Real energy = 0, force = 0;
      Real rDistSquared = (TNonbondedForce::DIST_R2 ? 1.0 / distSquared : 1.0);
      nonbondedForceFunction(energy, force, distSquared, rDistSquared, diff,
                             realTopo, i, j, excl);
      //      cout << "EN: " << energy << " FO: " << force << endl;
      // Calculate the switched force and energy.
      if (TSwitchingFunction::MODIFY) {
        Real switchingValue, switchingDeriv;
        switchingFunction(switchingValue, switchingDeriv, distSquared);
        // This has a - sign because the force is the negative of the
        // derivative of the energy (divided by the distance between the atoms).
        force = force * switchingValue - energy * switchingDeriv;
        energy = energy * switchingValue;
      }

      // Add this energy into the total system energy.
      nonbondedForceFunction.accumulateEnergy(energies, energy);
      // Add this force into the atom forces.
      Vector3D fij(diff * force);
      (*forces)[i] -= fij;
      (*forces)[j] += fij;

      // compute the vector between molecular centers of mass
      if (!same && energies->molecularVirial())
        // Add to the atomic and molecular virials
        energies->
          addVirial(fij, diff, realTopo->boundaryConditions.
                    minimalDifference(realTopo->molecules[mi].position,
                                      realTopo->molecules[mj].position));
      else if (energies->virial())
        energies->addVirial(fij, diff);
      // End of force computation.
      if (TConstraint::POST_CHECK)
        TConstraint::check(realTopo, i, j, diff, energy, fij);
    }

    void getParameters(std::vector<Parameter> &parameters) const {
      nonbondedForceFunction.getParameters(parameters);
      switchingFunction.getParameters(parameters);
    }

    static OneAtomPairNoExclusion make(std::vector<Value> values) {
      unsigned int n = TNonbondedForce::getParameterSize();

      std::vector<Value> parmsNF(values.begin(), values.begin() + n);
      std::vector<Value> parmsSF(values.begin() + n, values.end());

      return OneAtomPairNoExclusion(TNonbondedForce::make(parmsNF),
                         TSwitchingFunction::make(parmsSF));
    }

    static std::string getId() {
      return TConstraint::getPrefixId() + TNonbondedForce::getId() +
        TConstraint::getPostfixId() +
        std::string((!TSwitchingFunction::USE) ? std::string("") :
                    std::string(" -switchingFunction " +
                                TSwitchingFunction::getId()));
    }

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    TNonbondedForce *getNonbondedForceFunction() {
      return &nonbondedForceFunction;
    }

    TSwitchingFunction &getSwitchingFunction() {return switchingFunction;}

  private:
    mutable TopologyType *realTopo;
    const Vector3DBlock *positions;
    Vector3DBlock *forces;
    ScalarStructure *energies;
    TSwitchingFunction switchingFunction;
    TNonbondedForce nonbondedForceFunction;
    Real mySquaredCutoff;
  };
}
#endif /* ONEATOMPAIRNOEXCLUSION_H */
