/* -*- c++ -*- */
#ifndef NONBONDEDCUTOFFSYSTEMFORCE_H
#define NONBONDEDCUTOFFSYSTEMFORCE_H

#include <protomol/force/system/SystemForce.h>
#include <protomol/force/nonbonded/NonbondedCutoffForce.h>

namespace ProtoMol {
  //____ NonbondedCutoffSystemForce

  template<class TCellManager, class TOneAtomPair>
  class NonbondedCutoffSystemForce :
    public NonbondedCutoffForce<TCellManager, TOneAtomPair, SystemForce,
                                NonbondedCutoffSystemForce<TCellManager,
                                                           TOneAtomPair> > {
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Typedef
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  protected:
    typedef typename TOneAtomPair::BoundaryConditions BoundaryConditions;
    typedef Topology<BoundaryConditions, TCellManager> RealTopologyType;
    typedef typename RealTopologyType::Enumerator EnumeratorType;
    typedef CellPair CellPairType;
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    NonbondedCutoffSystemForce() :
      NonbondedCutoffForce<TCellManager, TOneAtomPair, SystemForce,
                           NonbondedCutoffSystemForce>() {}

    NonbondedCutoffSystemForce(Real cutoff, TOneAtomPair oneAtomPair) :
      NonbondedCutoffForce<TCellManager, TOneAtomPair, SystemForce,
                           NonbondedCutoffSystemForce>(cutoff, oneAtomPair) {}

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class SystemForce
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    typedef NonbondedCutoffForce<TCellManager, TOneAtomPair, SystemForce,
                                 NonbondedCutoffSystemForce<TCellManager,
                                                            TOneAtomPair> >
    Super_T;
    using Super_T::evaluate; // Avoid compiler warning/error
    virtual void evaluate(const GenericTopology *topo,
                          const Vector3DBlock *positions,
                          Vector3DBlock *forces, ScalarStructure *energies) {
      const RealTopologyType *realTopo =
          (RealTopologyType*) topo;
      //  dynamic_cast<const RealTopologyType *>(topo);
      this->myOneAtomPair.initialize(realTopo, positions, forces, energies);
      realTopo->updateCellLists(positions);
      this->enumerator.initialize(realTopo, this->myCutoff);
      this->doEvaluate(topo, realTopo->cellLists.size());
    }

    virtual void parallelEvaluate(const GenericTopology *topo,
                                  const Vector3DBlock *positions,
                                  Vector3DBlock *forces,
                                  ScalarStructure *energies) {
      const RealTopologyType *realTopo =
        dynamic_cast<const RealTopologyType *>(topo);
      
      this->myOneAtomPair.initialize(realTopo, positions, forces, energies);
      realTopo->updateCellLists(positions);
      this->enumerator.initialize(realTopo, this->myCutoff);
      
      unsigned int n = realTopo->cellLists.size();
      unsigned int count = numberOfBlocks(realTopo, positions);
      
      for (unsigned int i = 0; i < count; i++) {
        unsigned int l = (n * (i + 1)) / count - (n * i) / count;

        if (Parallel::next()) this->doEvaluate(topo, l);
        else this->enumerator.nextNewPair(l);
      }
    }
  };
}
#endif /* NONBONDEDCUTOFFSYSTEMFORCE_H */
