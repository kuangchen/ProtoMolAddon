/* -*- c++ -*- */
#ifndef COMPAREFORCE_H
#define COMPAREFORCE_H

#include <protomol/force/Force.h>
#include <vector>

namespace ProtoMol {
  //________________________________________ CompareForce

  class CompareForce : virtual public Force {
    // This class contains the definition of one force

    struct CompareError {
      CompareError() : absF2(0.0), rFavg(0.0), rFmax(0.0), rPE(0.0) {}
      CompareError(Real a, Real b, Real c,
                   Real d) : absF2(a), rFavg(b), rFmax(c), rPE(d) {}
      Real absF2, rFavg, rFmax, rPE;
    };

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    CompareForce(Force *actualForce, CompareForce *compareForce);
    virtual ~CompareForce();
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class CompareForce
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    const Vector3DBlock *getForces() const {return myForces;}
    const ScalarStructure *getEnergies() const {return myEnergies;}
    Force *getForceObject() const {return myActualForce;}
    unsigned int getIdNumber() const {return myIdNumber;}
  protected:
    void preprocess(unsigned int numAtoms);
    void postprocess(const GenericTopology *topo,
                     Vector3DBlock *forces,
                     ScalarStructure *energies);
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class Force
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    virtual std::string getKeyword() const {return keyword;}
    virtual unsigned int numberOfBlocks(const GenericTopology *topo,
                                        const Vector3DBlock *pos);
  private:
    virtual Force *doMake(const std::vector<Value> &values) const;
    virtual void uncache();

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class Makable
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    virtual std::string getIdNoAlias() const;
    virtual void getParameters(std::vector<Parameter> &parameters) const;

  private:
    virtual void doSetParameters(std::vector<Value> values);

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    static const std::string keyword;
  protected:
    Force *myActualForce;
    Vector3DBlock *myForces;
    ScalarStructure *myEnergies;
  private:
    CompareForce *myCompareForce;
    unsigned int myIdNumber;
    std::string myForcename;
    std::string myCompareForcename;
    static unsigned int myCounter;
    std::vector<CompareError> myErrors;
  };

  //________________________________________ INLINES
}
#endif /* COMPAREFORCE_H */
