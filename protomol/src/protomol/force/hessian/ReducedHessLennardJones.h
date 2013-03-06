/*  -*- c++ -*-  */
#ifndef REDUCEDHESSLENNARDJONES_H
#define REDUCEDHESSLENNARDJONES_H

#include <protomol/type/Matrix3By3.h>
#include <protomol/topology/ExclusionTable.h>

/*
 *  LennardJonesEnergy = LennardJonesEnergy * switchingHess.
 *___________________________________________________ReducedHessLennardJones
 *  Now one of the atoms, heavy atom, in the nonbonded interaction is
 *  considered as fixed so that only the light atom position remain
 *  independent. This is used as the heavy-atom-fixed approximation in H-bond
 *  MOLLY integrator.
 *_____________________________________________________________________*/

namespace ProtoMol {
  class GenericTopology;

  class ReducedHessLennardJones {
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class ReducedHessLennardJones
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    Matrix3By3 operator()(const Real rawEnergy, const Real rawForce,
                          Real distSquared, Real rDistSquared,
                          const Vector3D &diff,
                          const GenericTopology *topo,
                          int atom1, int atom2, const Real switchingValue,
                          const Real switchingDeriv,
                          const Matrix3By3 &switchingHess,
                          ExclusionClass excl) const;
  };
}
#endif
