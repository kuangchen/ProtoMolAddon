/*  -*- c++ -*-  */
#ifndef REDUCEDHESSCOULOMBDIELEC_H
#define REDUCEDHESSCOULOMBDIELEC_H

#include <protomol/type/Matrix3By3.h>
#include <protomol/topology/ExclusionTable.h>

namespace ProtoMol {
  class GenericTopology;

  class ReducedHessCoulombDielec {
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class ReducedHessCoulombDielec
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    Matrix3By3 operator()(const Real rawEnergy, const Real rawForce,
                          Real distSquared, Real rDistSquared,
                          const Vector3D &diff,
                          const GenericTopology *topo,
                          int atom1, int atom2, const Real switchingValue,
                          const Real switchingDeriv,
                          const Matrix3By3 &switchingHess,
                          ExclusionClass excl, const Real D, const Real S,
                          const Real epsi) const;
  };
}
#endif
