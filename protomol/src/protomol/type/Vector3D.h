/*  -*- c++ -*-  */
#ifndef VECTOR3D_H
#define VECTOR3D_H

#include <protomol/base/MathUtilities.h>
#include <protomol/base/Report.h>

#include <iostream>
using namespace std;


namespace ProtoMol {
  //_________________________________________________________________ Vector3D
  /**
   * Container to hold 3D vector/coordinate
   */


  template<class T>
  struct Vector3DImpl {
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
public:
    T c;

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
public:
    Vector3DImpl(Real* p) {c=p;c[0]=0;c[1]=0;c[2]=0;}

    Vector3DImpl(Real X, Real Y, Real Z, Real* p) {c=p;c[0]=X;c[1]=Y;c[2]=Z;}

    template<class U>
    Vector3DImpl(const Vector3DImpl<U> & rhs, Real* p) {
      c = p;
      c[0]=rhs.c[0];
      c[1]=rhs.c[1];
      c[2]=rhs.c[2];
    }


    Vector3DImpl() {c[0]=0;c[1]=0;c[2]=0;}  
    Vector3DImpl(Real X, Real Y, Real Z) {c[0]=X;c[1]=Y;c[2]=Z;}

    template<class U>
    Vector3DImpl(const Vector3DImpl<U> & rhs) {
      c[0] = rhs.c[0];
      c[1] = rhs.c[1];
      c[2] = rhs.c[2];
    }


    ~Vector3DImpl() {}

    template<class U>
    Vector3DImpl &operator=(const Vector3DImpl<U> &rhs) {
      c[0] = rhs.c[0];
      c[1] = rhs.c[1];
      c[2] = rhs.c[2];
      return *this;
    }



    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class Vector3DImpl
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
public:
    /// Index access
    Real operator[](int index) const {
      if (index < 0 || index > 2)
        Report::report << Report::error <<
        "[Vector3DImpl::operator[] const] index out of range" << Report::endr;
      return c[index];
      //return index == 0 ? c[0] : (index == 1 ? c[1] : c[2]);
    }

    /// Index access
    Real &operator[](int index) {
      if (index < 0 || index > 2)
        Report::report << Report::error <<
        "[Vector3DImpl::operator[] const] index out of range" << Report::endr;
      return c[index];
      //return index == 0 ? c[0] : (index == 1 ? c[1] : c[2]);
    }


    // Binary operators
    template<class U>
    Vector3DImpl<Real[3]> operator+(const Vector3DImpl<U> &b) const {
      return Vector3DImpl<Real[3]>(c[0] + b.c[0], c[1] + b.c[1], c[2] + b.c[2]);
    }

    /* No change */
    template<class U>
    Vector3DImpl<Real[3]> add(const Vector3DImpl<U> &b) const {
      return (*this) + b;
    }

    template<class U>
    Vector3DImpl<Real[3]> operator-(const Vector3DImpl<U> &b) const {
      return Vector3DImpl<Real[3]>(c[0] - b.c[0], c[1] - b.c[1], c[2] - b.c[2]);
    }

    /* No change */
    template<class U>
    Vector3DImpl<Real[3]> subtract(const Vector3DImpl<U> &b) const {
      return (*this) - b;
    }

    /// dot product
    template<class U>
    Real operator*(const Vector3DImpl<U> &b) const {
      return c[0] * b.c[0] + c[1] * b.c[1] + c[2] * b.c[2];
    }
    /// dot product
    /* No change */
    template<class U>
    Real dot(const Vector3DImpl<U> &b) const {
      return (*this) * b;
    }

    /// cross product
    template<class U>
    Vector3DImpl<Real[3]> operator^(const Vector3DImpl<U> &b) const {
      return Vector3DImpl<Real[3]>(c[1] * b.c[2] - c[2] * b.c[1], c[2] * b.c[0] - c[0] * b.c[2], c[0] * b.c[1] - c[1] * b.c[0]);
    }
    /// cross product
    /* No change */
    template<class U>
    Vector3DImpl<Real[3]> cross(const Vector3DImpl<U> &b) const {
      return (*this) ^ b;
    }


    Vector3DImpl<Real[3]> operator*(Real w) const {
      return Vector3DImpl<Real[3]>(c[0] * w, c[1] * w, c[2] * w);
    }
    /* No change */
    Vector3DImpl<Real[3]> multiply(Real w) const {
      return (*this) * w;
    }

    Vector3DImpl<Real[3]> operator/(Real w) const {
      return Vector3DImpl<Real[3]>(c[0] / w, c[1] / w, c[2] / w);
    }

    /* No change */
    Vector3DImpl<Real[3]> divide(Real w) const {
      return (*this) / w;
    }


    // Unary operators
    Vector3DImpl<Real[3]> operator-() const {
      return Vector3DImpl<Real[3]>(-c[0], -c[1], -c[2]);
    }


    // Comparison
    template<class U>
    bool operator==(const Vector3DImpl<U> &b) const {
      return c[0] == b.c[0] && c[1] == b.c[1] && c[2] == b.c[2];
    }


    template<class U>
    bool operator!=(const Vector3DImpl<U> &b) const {
      return c[0] != b.c[0] || c[1] != b.c[1] || c[2] != b.c[2];
    }

    // Assignment
    template<class U>
    Vector3DImpl<T> &operator+=(const Vector3DImpl<U> &b) {
      c[0] += b.c[0];
      c[1] += b.c[1];
      c[2] += b.c[2];
      return *this;
    }

    /* No change */
    template<class U>
    Vector3DImpl<T> &intoAdd(const Vector3DImpl<U> &b) {
      return (*this) += b;
    }

    template<class U>
    Vector3DImpl<T> &operator-=(const Vector3DImpl<U> &b) {
      c[0] -= b.c[0];
      c[1] -= b.c[1];
      c[2] -= b.c[2];
      return *this;
    }

    /* No change */
    template<class U>
    Vector3DImpl<T> &intoSubtract(const Vector3DImpl<U> &b) {
      return (*this) -= b;
    }


    Vector3DImpl<T> &operator*=(Real w) {
      c[0] *= w;
      c[1] *= w;
      c[2] *= w;
      return *this;
    }

    /* No change */
    Vector3DImpl<T> &intoMultiply(Real w) {
      return (*this) *= w;
    }

    Vector3DImpl<T> &operator/=(Real w) {
      c[0] /= w;
      c[1] /= w;
      c[2] /= w;
      return *this;
    }

    /* No change */
    Vector3DImpl<T> &intoDivide(Real w) {
      return (*this) /= w;
    }

    template<class U>
    Vector3DImpl<T> &intoWeightedAdd(Real w, const Vector3DImpl<U> &b) {
      c[0] += w * b.c[0];
      c[1] += w * b.c[1];
      c[2] += w * b.c[2];
      return *this;
    }

    template<class U>
    Vector3DImpl<T> &intoWeightedSubtract(Real w, const Vector3DImpl<U> &b) {
      c[0] -= w * b.c[0];
      c[1] -= w * b.c[1];
      c[2] -= w * b.c[2];
      return *this;
    }


    Real normSquared() const {
      return c[0] * c[0] + c[1] * c[1] + c[2] * c[2];
    }

    Real norm() const {
      return sqrt(c[0] * c[0] + c[1] * c[1] + c[2] * c[2]);
    }

    /// Normalize the Vector3DImpl and return the original length.
    Real normalize() {
      Real len = norm();

      if (len == 0.0)
        return 0.0;

      Real d = 1 / len;
      c[0] *= d; c[1] *= d; c[2] *= d;
      return len;
    }


    /// Return a normalized Vector3DImpl, leave the original Vector3DImpl unchanged.
    Vector3DImpl<Real[3]> normalized() const {
      Real len = norm();

      if (len == 0.0) {
        Report::report << Report::recoverable <<
        "[Vector3DImpl::normalized] length is zero." << Report::endr;
        return Vector3DImpl<Real[3]>(0.0, 0.0, 0.0);
      }

      return Vector3DImpl<Real[3]>(*this / len);
    }


    friend std::ostream &operator<<(std::ostream &OS, const Vector3DImpl &coords) {
      OS << coords.c[0] << " " << coords.c[1] << " " << coords.c[2] << " ";
      return OS;
    }

    friend std::istream &operator>>(std::istream &OS, Vector3DImpl &coords) {
      OS >> coords.c[0] >> coords.c[1] >> coords.c[2];
      return OS;
    }

    friend Report::MyStreamer &operator<<(Report::MyStreamer &OS,
                                          const Vector3DImpl &coords) {
      OS << "(" << coords.c[0] << "," << coords.c[1] << "," << coords.c[2] << ")";
      return OS;
    }
  };


  typedef Vector3DImpl<Real[3]> Vector3D;
  typedef Vector3DImpl<Real*> Vector3DB;
}

#endif

