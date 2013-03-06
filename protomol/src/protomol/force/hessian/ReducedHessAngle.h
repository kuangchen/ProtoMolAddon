/* -*- c++ -*- */
#ifndef REDUCEDHESSANGLE_H
#define REDUCEDHESSANGLE_H

#include <protomol/type/Matrix3By3.h>


namespace ProtoMol {
  /**
   *
   *  Note that the U-B part of the angle energy is omitted because in the
   *  MOLLY averaging, we do not incorporate the U-B part.@n
   *
   *  The remaining computation is same as that in AngleHessian class except
   *  now a more compact form is used to compute each elements.
   *  See Appendix B in Qun Ma's dissertation for details.@n
   *
   *  It can compute the whole Hessian or only the blocks of (0,0), (0,2),
   *  (2,0), (2,2) for the "reduced" form when the Heavy atom is "anchored."
   */
  class ReducedHessAngle {
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    ReducedHessAngle();
    ReducedHessAngle(const Vector3D &atom_i, ///< position of atom i
                     const Vector3D &atom_j, ///< position of atom j
                     const Vector3D &atom_k, ///< position of atom k
                     const double k_t, ///< angluar spring constant
                     const double theta0, ///< rest angle
                     bool computeReduced = false);
    // spring constants, and rest length, rest angle: from Topology
    //ReducedHessAngle(const ReducedHessAngle &H);
    ~ReducedHessAngle() {}

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class ReducedHessAngle
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    void evaluate(const Vector3D &P1,
                  const Vector3D &P2,
                  const Vector3D &P3,
                  const double kt,
                  const double t0,
                  bool computeReduced = 0);

  public:
    /// Extract an element.
    Matrix3By3 operator()(int i, int j) const;

    /// Set an element.
    void operator()(int i, int j, Matrix3By3 x);

    void accumulateTo(int i, int j, Matrix3By3 x);
    void accumulateNegTo(int i, int j, Matrix3By3 x);

    ReducedHessAngle operator*(const ReducedHessAngle &tm);
    ReducedHessAngle operator*(const Real tm);
    ReducedHessAngle operator/(const Real tm);
    ReducedHessAngle &operator*=(const Real tm);
    ReducedHessAngle &operator/=(const Real tm);
    ReducedHessAngle operator+(const ReducedHessAngle &tm);
    ReducedHessAngle operator-(const ReducedHessAngle &tm);
    ReducedHessAngle &operator+=(const ReducedHessAngle &tm);
    ReducedHessAngle &operator-=(const ReducedHessAngle &tm);
    // operator overloading

    ReducedHessAngle transposed(); ///< the original object unchanged
    void convertFromJacobian(Real **jac, int n); ///< n is the dimension
    void identity();  ///< identity matrix
    void clear();     ///< clear the hessian matrix

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Friends of class ReducedHessAngle
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    friend std::ostream &operator<<(std::ostream &os,
                                    const ReducedHessAngle &tm);

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // private data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  private:
    Matrix3By3 H[3][3];
    //in some applications, we need to have a 3 X 3 block matrix,
    //which may not be a symmetric matrix. so we want to have the
    //lower half blocks too.
  };
}
#endif
