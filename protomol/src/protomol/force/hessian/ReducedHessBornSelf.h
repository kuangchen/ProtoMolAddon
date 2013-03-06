/*  -*- c++ -*-  */
#ifndef REDUCEDHESSBORNSELF_H
#define REDUCEDHESSBORNSELF_H

#include <protomol/type/Matrix3By3.h>
#include <protomol/topology/ExclusionTable.h>

namespace ProtoMol {
  class GenericTopology;

  class ReducedHessBornSelf {
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // New methods of class ReducedHessBornSelf
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    Matrix3By3 operator()(Real distSquared,
                          const Vector3D &diff,
                          const GenericTopology *topo,
                          int atom1, int atom2, 
                          int bornSwitch, Real dielecConst,
                          ExclusionClass excl) const;


  };
}
#endif
