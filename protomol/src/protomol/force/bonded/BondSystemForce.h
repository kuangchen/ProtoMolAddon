/* -*- c++ -*- */
#ifndef BONDSYSTEMFORCE_H
#define BONDSYSTEMFORCE_H

#include <protomol/force/system/SystemForce.h>
#include <protomol/type/ScalarStructure.h>
#include <protomol/parallel/Parallel.h>

#include <string>

namespace ProtoMol {
  //____ BondSystemForce

  template<class TBoundaryConditions>
  class BondSystemForce : public SystemForce {
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    virtual ~BondSystemForce() {}

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class BondSystemForce
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    void calcBond(const TBoundaryConditions &boundary, const Bond &currentBond,
                  const Vector3DBlock *positions, Vector3DBlock *forces,
                  ScalarStructure *energies);

    Real calcBondEnergy(const TBoundaryConditions &boundary,
                        const Bond &currentBond,
                        const Vector3DBlock *positions);

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class SystemForce
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    using SystemForce::evaluate; // Avoid compiler warning/error
    virtual void evaluate(const GenericTopology *topo,
                          const Vector3DBlock *positions, Vector3DBlock *forces,
                          ScalarStructure *energies);

    virtual void parallelEvaluate(const GenericTopology *topo,
                                  const Vector3DBlock *positions,
                                  Vector3DBlock *forces,
                                  ScalarStructure *energies);

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class Force
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    virtual std::string getKeyword() const {return "Bond";}

    virtual unsigned int numberOfBlocks(const GenericTopology *topo,
                                        const Vector3DBlock *pos);

  private:
    virtual Force *doMake(const std::vector<Value> &) const {
      return new BondSystemForce();
    }

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class Makeable
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    virtual std::string getIdNoAlias() const {return getKeyword();}
    virtual void getParameters(std::vector<Parameter> &) const {}

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  private:
  };

  //____ INLINES

  template<class TBoundaryConditions>
  inline void BondSystemForce<TBoundaryConditions>::evaluate(
    const GenericTopology *topo, const Vector3DBlock *positions,
    Vector3DBlock *forces, ScalarStructure *energies) {
    const TBoundaryConditions &boundary =
      ((SemiGenericTopology<TBoundaryConditions> &)(*topo)).
        boundaryConditions;

    for (unsigned int i = 0; i < topo->bonds.size(); i++)
      calcBond(boundary, topo->bonds[i], positions, forces, energies);
  }

  template<class TBoundaryConditions>
  inline void BondSystemForce<TBoundaryConditions>::calcBond(
    const TBoundaryConditions &boundary, const Bond &currentBond,
    const Vector3DBlock *positions, Vector3DBlock *forces,
    ScalarStructure *energies) {
    int a1 = currentBond.atom1;
    int a2 = currentBond.atom2;

    Real restLength = currentBond.restLength;
    Real springConstant = currentBond.springConstant;

    Vector3D atom1((*positions)[a1]);
    Vector3D atom2((*positions)[a2]);

    // Vector from atom 1 to atom 2.
    Vector3D r12(boundary.minimalDifference(atom2, atom1));
    Real r = r12.norm();                      // Distance between atom 1 and 2.

    Real dpotdr = 2.0 * springConstant * (r - restLength);  // Calculate dpot/dr

    // Calculate force on atom1 due to atom2.
    Vector3D force1(r12 * (-dpotdr / r));

    // Add to the total force.
    (*forces)[a1] += force1;
    (*forces)[a2] -= force1;

    // Add energy
    (*energies)[ScalarStructure::BOND] += springConstant *
                                          (r - restLength) * (r - restLength);

    // Add virial
    if (energies->virial())
      energies->addVirial(force1, r12);
  }

  template<class TBoundaryConditions>
  inline Real BondSystemForce<TBoundaryConditions>::calcBondEnergy(
    const TBoundaryConditions &boundary, const Bond &currentBond,
    const Vector3DBlock *positions) {
    int a1 = currentBond.atom1;
    int a2 = currentBond.atom2;
    Real restLength = currentBond.restLength;
    Real springConstant = currentBond.springConstant;

    Vector3D atom1 = (*positions)[a1];
    Vector3D atom2 = (*positions)[a2];

    // Vector from atom 1 to atom 2.
    Vector3D r12 = boundary.minimalDifference(atom2, atom1);
    Real r = r12.norm();                      // Distance between atom 1 and 2.

    //Calculate energy.
    return springConstant * (r - restLength) * (r - restLength);
  }

  template<class TBoundaryConditions>
  inline void BondSystemForce<TBoundaryConditions>::parallelEvaluate(
    const GenericTopology *topo, const Vector3DBlock *positions,
    Vector3DBlock *forces, ScalarStructure *energies) {
    const TBoundaryConditions &boundary =
      (dynamic_cast<const SemiGenericTopology<TBoundaryConditions> &>(*topo)).
        boundaryConditions;

    unsigned int n = topo->bonds.size();
    unsigned int count = numberOfBlocks(topo, positions);

    for (unsigned int i = 0; i < count; i++)
      if (Parallel::next()) {
        int to = (n * (i + 1)) / count;
        if (to > static_cast<int>(n))
          to = n;
        int from = (n * i) / count;
        for (int j = from; j < to; j++)
          calcBond(boundary, topo->bonds[j], positions, forces, energies);
      }
  }

  template<class TBoundaryConditions>
  inline unsigned int BondSystemForce<TBoundaryConditions>::numberOfBlocks(
    const GenericTopology *topo, const Vector3DBlock *) {
    return std::min(Parallel::getAvailableNum(),
                    static_cast<int>(topo->bonds.size()));
  }
}

#endif /* BONDSYSTEMFORCE_H */
