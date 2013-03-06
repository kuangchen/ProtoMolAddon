/* -*- c++ -*- */
#ifndef HESSIAN_H
#define HESSIAN_H

#include <protomol/force/Force.h>
#include <protomol/type/Vector3DBlock.h>
#include <protomol/type/Matrix3By3.h>

namespace ProtoMol {
  /**
   *
   * Calculates the Hessian or mass re-weighted hessian
   * for the current force field.
   *
   */
  class Hessian {
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Types and Enums
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  protected:
    enum {COULOMB = 1};
    enum {COULOMBDIELEC = 2};
    enum {COULOMBSCPISM = 4};
    enum {LENNARDJONES = 8};
    enum {BORNCUTOFF2 = 25};

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    Hessian();
    Hessian(unsigned int szin);   //dimension of matrix (sz*sz)
    Hessian(const Hessian &hess); //copy const.
    ~Hessian();

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class Hessian
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    void initialData(unsigned int szin);
    void findForces(ForceGroup *overloadedForces);
    void evaluate(const Vector3DBlock *myPositions,   //positions
                  GenericTopology *myTopo,      //topology
                  const bool mrw);                    //mass re-weighted

    Matrix3By3 evaluatePairsMatrix(int i, int j, int pairType, const Vector3DBlock *myPositions,
                       const GenericTopology *myTopo, bool mrw);
    void evaluatePairs(int i, int j, int pairType, const Vector3DBlock *myPositions,
                       const GenericTopology *myTopo, bool mrw, 
                       int mat_i, int mat_j, int mat_sz, double * mat_array);
    void outputSparsePairMatrix(int i, int j, Real massi, Real massj, Matrix3By3 rha, 
                                bool mrw, int arrSz, double * basePoint);
    void outputSparseMatrix(int i, int j, Real massi, Real massj, Matrix3By3 rha, 
                            bool mrw, int arrSz, double * basePoint);

    void evaluateBornRadii(const Vector3DBlock *myPositions, GenericTopology *myTopo);

    bool setHessianColumn( const Vector3DBlock &hescol, const unsigned int columnNumber, 
                            const GenericTopology *myTopo, const bool massWeight );

    //GB specific code
    void evaluateGBBornBurialTerm(const Vector3DBlock *myPositions, GenericTopology *myTopo);
    void evaluateGBBornRadii(const Vector3DBlock *myPositions, GenericTopology *myTopo);

    Matrix3By3 evaluateBornSelfPair(int i, int j, const Vector3DBlock *myPositions,
                                       const GenericTopology *myTopo);

    Matrix3By3 evaluateGBACEPair(int i, int j, const Vector3DBlock *myPositions,
                                        const GenericTopology *myTopo);

    Matrix3By3 evaluateGBPair(int i, int j, const Vector3DBlock *myPositions,
                                  const GenericTopology *myTopo);

  public:
    void clear(); // clear the hessian matrix

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Friends of class Hessian
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // private data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  protected:
    bool myBond;    //force field components
    bool myAngle;
    bool myCoulomb;
    bool myCoulombDielec;
    bool myCoulombSCPISM;
    //bool myCoulombBornRadii;
    bool myBornRadii;
    bool myBornSelf;
    bool myLennardJones;
    bool myDihedral;
    bool myImproper;
    bool myRBDihedral;

    //GB stuff
    bool myGBBornBurialTerm, myGBBornRadii;
    bool myGBACEForce;
    Real solvationparam, watersphereradius;

    bool myGBForce;
    Real soluteDielec, solventDielec;

    Real cCutoff, cSwitchon, cSwitch, cOrder, cSwitchoff;   //switch data
    Real lCutoff, lSwitchon, lSwitch, lOrder, lSwitchoff;
    Real D, S, epsi;
    unsigned int sz; //size
    int myBornSwitch;
    Real myDielecConst;
  public:
    double *hessM;  //matrix
    Real cutOff;
  };
}
#endif
