/* -*- c++ -*- */
#ifndef HESSDIHEDRAL_H
#define HESSDIHEDRAL_H

#include <protomol/topology/GenericTopology.h>
#include <protomol/type/Matrix3By3.h>

namespace ProtoMol {
  /**
   *
   *  Calculates Hessian for Diherdal PE term.
   *  Rotates the dihedral into a bi-planar configuration
   *  to reduct the number of calculations. Then rotates
   *  the results back again.
   *  Extended to Rykart-Bellman Dihedrals
   */
  class HessDihedral {
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    HessDihedral();

    //torsion data including positions etc
    HessDihedral(const Torsion &currTorsion,
                 const Vector3D &a1, const Vector3D &a2,
                 const Vector3D &a3, const Vector3D &a4);

    //RBtorsion data including positions etc
    HessDihedral(const RBTorsion &currTorsion,
                 const Vector3D &a1, const Vector3D &a2,
                 const Vector3D &a3, const Vector3D &a4);

    ~HessDihedral() {}

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class HessDihedral
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    void evaluate(const Torsion &currTorsion, const Vector3D &a1,
                  const Vector3D &a2, const Vector3D &a3,
                  const Vector3D &a4);
    void evaluate(const RBTorsion &currRBTorsion, const Vector3D &a1,
                  const Vector3D &a2, const Vector3D &a3,
                  const Vector3D &a4);
    Matrix3By3 operator()(const unsigned int i, const unsigned int j) const;

  private:
    bool rotateAndCalcPhi(const Vector3D &a1, const Vector3D &a2, 
                            const Vector3D &a3, const Vector3D &a4,
                            Real &x21, Real &x43, Real &y43, Real &z21, Real &z43, Real &r23, 
                            double *aRot, Real &g);

    void outputHessian(const double fact1, const double fact2, 
                       const Real x21, const Real x43, const Real y43, const Real z21, const Real z43, const Real r23, 
                       const double *aRot);

    //Use aRot to rotate the vector back into real space
    double *rotateV3D(const double *aRot, double *mf);

    //Calculates determinate of 3x3 matrix
    double det(const double* a);

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Friends of class HessDihedral
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    double hessD[144];  //storage for 4x4 hessian matrices in 3D

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // private data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  private:
  };
}
#endif
