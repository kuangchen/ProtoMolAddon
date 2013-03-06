/*  -*- c++ -*-  */
#ifndef REDUCEDHESSGBACE_H
#define REDUCEDHESSGBACE_H
#include <protomol/type/Matrix3By3.h>
#include <protomol/topology/ExclusionTable.h>

namespace ProtoMol {
  class GenericTopology;

  class ReducedHessGBACE {
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // New methods of class ReducedHessGBACE
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
      Matrix3By3 operator()(Real distSquared,
                          const Vector3D &diff,
                          const GenericTopology *topo,
                          int atom1, int atom2,
                          Real de, Real s_p,
                          ExclusionClass excl) const;


  };

}
#endif /* REDUCEDHESSGBACE_H */
