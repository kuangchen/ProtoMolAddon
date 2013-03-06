/* -*- c++ -*- */
#ifndef DIHEDRALSYSTEMFORCE_H
#define DIHEDRALSYSTEMFORCE_H

#include <protomol/force/bonded/MTorsionSystemForce.h>
#include <protomol/type/ScalarStructure.h>
#include <protomol/parallel/Parallel.h>
#include <protomol/topology/SemiGenericTopology.h>

#include <string>

namespace ProtoMol {
  //____ DihedralSystemForce

  template<class TBoundaryConditions>
  class DihedralSystemForce : public MTorsionSystemForce<TBoundaryConditions> {
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    virtual ~DihedralSystemForce() {}

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class SystemForce
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    // Avoid compiler warning/error
    using MTorsionSystemForce<TBoundaryConditions>::evaluate;
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
    virtual std::string getKeyword() const {return "Dihedral";}

    virtual unsigned int numberOfBlocks(const GenericTopology *topo,
                                        const Vector3DBlock *pos);

  private:
    virtual Force *doMake(const std::vector<Value> &) const {
      return new DihedralSystemForce();
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
  inline void DihedralSystemForce<TBoundaryConditions>::
  evaluate(const GenericTopology *topo, const Vector3DBlock *positions,
           Vector3DBlock *forces, ScalarStructure *energies) {

    const TBoundaryConditions &boundary =
      ((SemiGenericTopology<TBoundaryConditions> &)(*topo)).
        boundaryConditions;
    for (unsigned int i = 0; i < topo->dihedrals.size(); i++)
      calcTorsion(boundary, topo->dihedrals[i], positions, forces,
                  (*energies)[ScalarStructure::DIHEDRAL], energies);
  }

  template<class TBoundaryConditions>
  inline void DihedralSystemForce<TBoundaryConditions>::
  parallelEvaluate(const GenericTopology *topo, const Vector3DBlock *positions,
                   Vector3DBlock *forces, ScalarStructure *energies) {
    const TBoundaryConditions &boundary =
      (dynamic_cast<const SemiGenericTopology<TBoundaryConditions> &>(*topo)).
        boundaryConditions;

    unsigned int n = topo->dihedrals.size();
    unsigned int count = numberOfBlocks(topo, positions);

    for (unsigned int i = 0; i < count; i++)
      if (Parallel::next()) {
        int to = (n * (i + 1)) / count;
        if (to > static_cast<int>(n))
          to = n;
        int from = (n * i) / count;
        for (int j = from; j < to; j++)
          calcTorsion(boundary, topo->dihedrals[j], positions, forces,
                      (*energies)[ScalarStructure::DIHEDRAL], energies);
      }
  }

  template<class TBoundaryConditions>
  inline unsigned int DihedralSystemForce<TBoundaryConditions>::
  numberOfBlocks(const GenericTopology *topo, const Vector3DBlock *) {
    return std::min(Parallel::getAvailableNum(),
                    static_cast<int>(topo->dihedrals.size()));
  }
}
#endif /* DIHEDRALSYSTEMFORCE_H */
