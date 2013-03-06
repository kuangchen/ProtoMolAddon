/* -*- c++ -*- */
#ifndef NONBONDEDCUTOFFFORCE_H
#define NONBONDEDCUTOFFFORCE_H

#include <protomol/force/Force.h>
#include <protomol/parallel/Parallel.h>
#include <protomol/topology/Topology.h>

namespace ProtoMol {
  //____ NonbondedCutoffForce

  template<class TCellManager, class TOneAtomPair, class TForce,
           class TImplForce>
  class NonbondedCutoffForce : public TForce {
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Typedef
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  private:
    typedef typename TOneAtomPair::BoundaryConditions BoundaryConditions;
    typedef Topology<BoundaryConditions, TCellManager> RealTopologyType;
    typedef typename RealTopologyType::Enumerator EnumeratorType;
    typedef CellPair CellPairType;

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    NonbondedCutoffForce() : myCutoff(0.0) {}

    NonbondedCutoffForce(Real cutoff, TOneAtomPair oneAtomPair) :
      TForce(), myCutoff(cutoff), myOneAtomPair(oneAtomPair) {}

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class NonbondedCutoffForce
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  protected:
    void doEvaluate(const GenericTopology *topo, unsigned int n) {
      CellPairType thisPair;
      unsigned int count = 0;
      for (; !enumerator.done(); enumerator.next()) {
        enumerator.get(thisPair);
        bool notSameCell = enumerator.notSameCell();
        
        if (!notSameCell) {
          count++;
          if (count > n) break;
        }
        for (int i = thisPair.first; i != -1; i = topo->atoms[i].cellListNext)
          for (int j =
                 (notSameCell ? thisPair.second : topo->atoms[i].cellListNext);
               j != -1; j = topo->atoms[j].cellListNext)
            myOneAtomPair.doOneAtomPair(i, j);
      }
    }

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class Force
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    virtual unsigned int numberOfBlocks(const GenericTopology *topo,
                                        const Vector3DBlock *positions) {
      const RealTopologyType *realTopo =
        dynamic_cast<const RealTopologyType *>(topo);
      realTopo->updateCellLists(positions);
      return Parallel::getNumberOfPackages(realTopo->cellLists.size());
    }

    virtual std::string getKeyword() const {return "NonbondedCutoff";}

  private:
    virtual Force *doMake(const std::vector<Value> &values) const {
      std::vector<Value> atomPairValues(values.begin(), values.end() - 1);
      return new TImplForce(values[values.size() - 1], 
                            TOneAtomPair::make(atomPairValues));
    }

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class Makeable
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    virtual std::string getIdNoAlias() const {
      return TOneAtomPair::getId() + " -algorithm " + getKeyword();
    }

    virtual void getParameters(std::vector<Parameter> &parameters) const {
      myOneAtomPair.getParameters(parameters);
      parameters.push_back
        (Parameter("-cutoff",
                   Value(myCutoff, ConstraintValueType::Positive()),
                   Text("algorithm cutoff")));
    }

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    protected:
       Real myCutoff;
       TOneAtomPair myOneAtomPair;
       EnumeratorType enumerator;
  };
}
#endif /* NONBONDEDCUTOFFFORCE_H */
