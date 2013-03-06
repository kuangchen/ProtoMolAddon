/* -*- c++ -*- */
#ifndef ONEATOMPAIRTHREE_H
#define ONEATOMPAIRTHREE_H

#include <protomol/topology/Topology.h>
#include <protomol/config/Parameter.h>
#include <protomol/force/OneAtomContraints.h>

namespace ProtoMol {
  //____ OneAtomPairThree

  /**
   * Computes the interaction for a given force between two atoms with the
   * template arguments defining the boundary conditions, three switching 
   * functions, three potentials and optional constraint.
   */
  template<class TBoundaryConditions, class TSwitchingFunctionFirst,
           class TNonbondedForceFirst, class TSwitchingFunctionSecond,
           class TNonbondedForceSecond, class TSwitchingFunctionThird,
           class TNonbondedForceThird, class TConstraint = NoConstraint>
  class OneAtomPairThree {
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Typedef & sub classes
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    typedef TBoundaryConditions BoundaryConditions;
    // Make the boundary conditions visible

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    OneAtomPairThree() {}

    OneAtomPairThree(TNonbondedForceFirst f1, TSwitchingFunctionFirst sF1,
                   TNonbondedForceSecond f2, TSwitchingFunctionSecond sF2,
                   TNonbondedForceThird f3, TSwitchingFunctionThird sF3 ) :
      switchingFunctionFirst(sF1), nonbondedForceFunctionFirst(f1),
      switchingFunctionSecond(sF2), nonbondedForceFunctionSecond(f2),
      switchingFunctionThird(sF3), nonbondedForceFunctionThird(f3),
      mySquaredCutoff(std::max(
                      std::max
                      (Cutoff<TNonbondedForceFirst::CUTOFF>::cutoff(sF1, f1),
                        Cutoff<TNonbondedForceSecond::CUTOFF>::cutoff(sF2, f2)),
                          Cutoff<TNonbondedForceThird::CUTOFF>::cutoff(sF3, f3))
                       )
    {}

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class OneAtomPairThree
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    void initialize(const SemiGenericTopology<TBoundaryConditions> *topo,
                    const Vector3DBlock *pos, Vector3DBlock *f,
                    ScalarStructure *e) {
      realTopo = topo;
      positions = pos;
      forces = f;
      energies = e;
    }

    // Computes the force and energy for atom i and j.
    void doOneAtomPair(const int i, const int j) {
      if (TConstraint::PRE_CHECK)
        if (!TConstraint::check(realTopo, i, j))
          return;

      // Get atom distance.
      Real distSquared;
      Vector3D diff(realTopo->boundaryConditions.
                    minimalDifference((*positions)[i], (*positions)[j],
                                      distSquared));
      // Do switching function rough test, if necessary.
      if ((TSwitchingFunctionFirst::USE || TSwitchingFunctionSecond::USE || TSwitchingFunctionThird::USE ||
           TNonbondedForceFirst::CUTOFF || TNonbondedForceSecond::CUTOFF || TNonbondedForceThird::CUTOFF ) 
           && distSquared > mySquaredCutoff)
        return;

      // Check for an exclusion.
      int mi = realTopo->atoms[i].molecule;
      int mj = realTopo->atoms[j].molecule;
      bool same = (mi == mj);
      ExclusionClass excl =
        (same ? realTopo->exclusions.check(i, j) : EXCLUSION_NONE);
      if (excl == EXCLUSION_FULL)
        return;

      // Calculate the force and energy.
      Real rDistSquared =
        ((TNonbondedForceFirst::DIST_R2 ||
            TNonbondedForceSecond::DIST_R2 ||
              TNonbondedForceThird::DIST_R2 ) ? 1.0 / distSquared : 1.0);
      Real energy1, force1, energy2 = 0, force2 = 0, energy3 = 0, force3 = 0;
      nonbondedForceFunctionFirst(energy1, force1, distSquared, rDistSquared,
                                  diff, realTopo, i, j, excl);
      nonbondedForceFunctionSecond(energy2, force2, distSquared, rDistSquared,
                                   diff, realTopo, i, j, excl);
      nonbondedForceFunctionThird(energy3, force3, distSquared, rDistSquared,
                                   diff, realTopo, i, j, excl);
      //Report::report << "\t"<<i << "\t"<<j<<Report::endr;
      // Calculate the switched force and energy.
      if (TSwitchingFunctionFirst::MODIFY || 
            TSwitchingFunctionSecond::MODIFY || 
              TSwitchingFunctionThird::MODIFY) {
        Real switchingValue, switchingDeriv;

        switchingFunctionFirst(switchingValue, switchingDeriv, distSquared);
        force1 = force1 * switchingValue - energy1 * switchingDeriv;
        energy1 = energy1 * switchingValue;

        switchingFunctionSecond(switchingValue, switchingDeriv, distSquared);
        force2 = force2 * switchingValue - energy2 * switchingDeriv;
        energy2 = energy2 * switchingValue;

        switchingFunctionThird(switchingValue, switchingDeriv, distSquared);
        force3 = force3 * switchingValue - energy3 * switchingDeriv;
        energy3 = energy3 * switchingValue;
      }

      // Add this energy into the total system energy.
      nonbondedForceFunctionFirst.accumulateEnergy(energies, energy1);
      nonbondedForceFunctionSecond.accumulateEnergy(energies, energy2);
      nonbondedForceFunctionThird.accumulateEnergy(energies, energy3);
      // Add this force into the atom forces.
      Vector3D fij(diff * (force1 + force2 + force3));
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
        TConstraint::check(realTopo, i, j, diff, energy1 + energy2 + energy3, fij);
    }

    void getParameters(std::vector<Parameter> &parameters) const {
      nonbondedForceFunctionFirst.getParameters(parameters);
      switchingFunctionFirst.getParameters(parameters);
      nonbondedForceFunctionSecond.getParameters(parameters);
      switchingFunctionSecond.getParameters(parameters);
      nonbondedForceFunctionThird.getParameters(parameters);
      switchingFunctionThird.getParameters(parameters);
    }

    static unsigned int getParameterSize() {
      return
        TNonbondedForceFirst::getParameterSize() +
        TSwitchingFunctionFirst::getParameterSize() +
        TNonbondedForceSecond::getParameterSize() +
        TSwitchingFunctionSecond::getParameterSize() +
        TNonbondedForceThird::getParameterSize() +
        TSwitchingFunctionThird::getParameterSize();
    }

    static OneAtomPairThree make(std::vector<Value> values) {
      unsigned int l1a = TNonbondedForceFirst::getParameterSize();
      unsigned int l1b = TSwitchingFunctionFirst::getParameterSize() + l1a;
      unsigned int l2a = TNonbondedForceSecond::getParameterSize() + l1b;
      unsigned int l2b = TSwitchingFunctionSecond::getParameterSize() + l2a;
      unsigned int l3a = TNonbondedForceThird::getParameterSize() + l2b;

      std::vector<Value> F1(values.begin(), values.begin() + l1a);
      std::vector<Value> S1(values.begin() + l1a, values.begin() + l1b);
      std::vector<Value> F2(values.begin() + l1b, values.begin() + l2a);
      std::vector<Value> S2(values.begin() + l2a, values.begin() + l2b);
      std::vector<Value> F3(values.begin() + l2b, values.begin() + l3a);
      std::vector<Value> S3(values.begin() + l3a, values.end());

      return OneAtomPairThree
        (TNonbondedForceFirst::make(F1), TSwitchingFunctionFirst::make(S1),
          TNonbondedForceSecond::make(F2), TSwitchingFunctionSecond::make(S2),
            TNonbondedForceThird::make(F3), TSwitchingFunctionThird::make(S3));
    }

    static std::string getId() {
      return
        TConstraint::getPrefixId() +
        headString(TNonbondedForceFirst::getId()) +
        TConstraint::getPostfixId() + " " + 
        TConstraint::getPrefixId() +
        headString(TNonbondedForceSecond::getId()) +
        TConstraint::getPostfixId() + " " + 
        TConstraint::getPrefixId() +
        headString(TNonbondedForceThird::getId()) +
        TConstraint::getPostfixId() +

        (tailString(TNonbondedForceFirst::getId()).empty() ? "" : " ") +
        tailString(TNonbondedForceFirst::getId()) +
        (tailString(TNonbondedForceSecond::getId()).empty() ? "" : " ") +
        tailString(TNonbondedForceSecond::getId()) +
        (tailString(TNonbondedForceThird::getId()).empty() ? "" : " ") +
        tailString(TNonbondedForceThird::getId()) +

        std::string((!TSwitchingFunctionFirst::USE) ?
                    std::string("") :
                    std::string(" -switchingFunction " +
                                TSwitchingFunctionFirst::getId())) +
        std::string((!TSwitchingFunctionSecond::USE) ?
                    std::string("") :
                    std::string(" -switchingFunction " +
                                TSwitchingFunctionSecond::getId())) +
        std::string((!TSwitchingFunctionThird::USE) ?
                    std::string("") :
                    std::string(" -switchingFunction " +
                                TSwitchingFunctionThird::getId()));

    }

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  private:
    const SemiGenericTopology<TBoundaryConditions> *realTopo;
    const Vector3DBlock *positions;
    Vector3DBlock *forces;
    ScalarStructure *energies;
    TSwitchingFunctionFirst switchingFunctionFirst;
    TNonbondedForceFirst nonbondedForceFunctionFirst;
    TSwitchingFunctionSecond switchingFunctionSecond;
    TNonbondedForceSecond nonbondedForceFunctionSecond;
    TSwitchingFunctionThird switchingFunctionThird;
    TNonbondedForceThird nonbondedForceFunctionThird;
    Real mySquaredCutoff;
  };
}

#endif /* ONEATOMPAIRTHREE_H */
