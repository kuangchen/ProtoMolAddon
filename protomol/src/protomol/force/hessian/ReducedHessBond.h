/*  -*- c++ -*-  */
#ifndef REDUCEDHESSBOND_H
#define REDUCEDHESSBOND_H

#include <protomol/type/Vector3D.h>
#include <protomol/type/Matrix3By3.h>

namespace ProtoMol {
  Matrix3By3 reducedHessBond(const Vector3D &atom_i, const Vector3D &atom_j,
                             const Real _k, const Real _r0);
}
#endif
