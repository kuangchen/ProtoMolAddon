/* -*- c++ -*- */
#ifndef NONBONDEDINTERMITTENTFULLSYSTEMFORCE_H
#define NONBONDEDINTERMITTENTFULLSYSTEMFORCE_H

#include <protomol/force/system/SystemForce.h>
#include <protomol/parallel/Parallel.h>
#include <protomol/base/MathUtilities.h>
#include <protomol/base/Exception.h>

namespace ProtoMol {
  //____ NonbondedIntermittentFullSystemForce

  template<class TOneAtomPair>
  class NonbondedIntermittentFullSystemForce : public SystemForce {
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Typedef
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  private:
    typedef SemiGenericTopology<typename TOneAtomPair::BoundaryConditions>
    RealTopologyType;

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    NonbondedIntermittentFullSystemForce() :
      SystemForce(), myBlockSize(0), myRate(0), myCached(false), myStep(0) {};
    NonbondedIntermittentFullSystemForce(TOneAtomPair oneAtomPair,
                                   unsigned int blockSize = defaultBlockSize,
                                         unsigned int rate = 1) :
      SystemForce(), myOneAtomPair(oneAtomPair), myBlockSize(blockSize), myRate(rate),
      myCached(false), myStep(0) {}

    virtual ~NonbondedIntermittentFullSystemForce() {};

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class SystemForce
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    using SystemForce::evaluate; // Avoid compiler warning/error
    virtual void evaluate(const GenericTopology *topo, const Vector3DBlock *pos,
                          Vector3DBlock *f, ScalarStructure *e) {
      myCached = true;
      doEvaluate(topo, pos, f, e, 0, topo->atoms.size(), 0, topo->atoms.size());
    }

    virtual void parallelEvaluate(const GenericTopology *topo,
                                  const Vector3DBlock *pos, Vector3DBlock *f,
                                  ScalarStructure *e) {
      if (!myCached)
        splitRangeArea(static_cast<unsigned int>(Parallel::getAvailableNum()),
                       0, topo->atoms.size(), myFromRange, myToRange);
      myCached = true;

      for (int i = 0; i < Parallel::getAvailableNum(); i++)
        if (Parallel::next())
          doEvaluate(topo, pos, f, e, myFromRange[i].first,
                     myFromRange[i].second, myToRange[i].first,
                     myToRange[i].second);
    }

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class Force
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    virtual unsigned int numberOfBlocks(const GenericTopology *,
                                        const Vector3DBlock *) {
      return Parallel::getAvailableNum();
    }

    virtual std::string getKeyword() const {return "NonbondedIntermittentFull";}
    virtual void uncache() {myCached = false;};

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class Makeable
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    virtual void getParameters(std::vector<Parameter> &parameters) const {
      myOneAtomPair.getParameters(parameters);
      parameters.push_back
        (Parameter("-blocksize",
                   Value(myBlockSize, ConstraintValueType::Positive()),
                   defaultBlockSize));
      parameters.push_back
      (Parameter("-rate",
                 Value(myRate, ConstraintValueType::Positive()),
                 1));
    }

    virtual std::string getIdNoAlias() const {
      return TOneAtomPair::getId() + " -algorithm " + getKeyword();
    }

  private:
    virtual Force *doMake(const std::vector<Value> &values) const {
      unsigned int blockSize;
      int n = values.size() - 2;
      values[n].get(blockSize);
      if (!(values[n].valid()) || blockSize == 0)
        THROW(getKeyword() + " algorithm: 0 < blocksize (=" +
              values[n].getString() + ").");
      unsigned int rate;
      int nr = values.size() - 1;
      values[nr].get(rate);
      if (!(values[nr].valid()) || rate == 0)
        THROW(getKeyword() + " algorithm: 0 < rate (=" +
              values[nr].getString() + ").");
      
      std::vector<Value> atomPairValues(values.begin(), values.end() - 1);

      return
        new NonbondedIntermittentFullSystemForce(TOneAtomPair::make(atomPairValues),
                                           blockSize, rate);
    }

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class NonbondedIntermittentFullSystemForce
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  private:
    void doEvaluate(const GenericTopology *topo, const Vector3DBlock *positions,
                    Vector3DBlock *forces, ScalarStructure *energies, int i0,
                    int i1, int j0, int j1) {
      
      if(!(myStep % myRate)){  //do if step is rate * n

        const RealTopologyType *realTopo =
          (const RealTopologyType *)(topo);
      
        myOneAtomPair.initialize(realTopo, positions, forces, energies);
      
        for (int blocki = i0; blocki < i1; blocki += myBlockSize) {
          int blocki_max = blocki;
          if (blocki_max < j0) blocki_max = j0;
          for (int blockj = blocki_max; blockj < j1; blockj += myBlockSize) {
            int istart = blocki;
            int iend = blocki + myBlockSize;
            if (iend > i1) iend = i1;
            for (int i = istart; i < iend; i++) {
              int jstart = blockj;
              if (jstart <= i) jstart = i + 1;
              int jend = blockj + myBlockSize;
              if (jend > j1) jend = j1;
              for (int j = jstart; j < jend; j++)
                myOneAtomPair.doOneAtomPair(i, j);
            }
          }
        }
      }
      myStep++; //update step number
    }

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  private:
    TOneAtomPair myOneAtomPair;
    unsigned int myBlockSize;
    unsigned int myRate;
    std::vector<PairUInt> myFromRange;
    std::vector<PairUInt> myToRange;
    bool myCached;
    
    //step counter
    unsigned int myStep;

    static const unsigned int defaultBlockSize = 64;
  };
}
#endif /* NONBONDEDINTERMITTENTFULLSYSTEMFORCE_H */
