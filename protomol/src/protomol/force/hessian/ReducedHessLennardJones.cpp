#include <protomol/force/hessian/ReducedHessLennardJones.h>
#include <protomol/topology/GenericTopology.h>

using namespace ProtoMol;

//____ ReducedHessLennardJones
Matrix3By3 ReducedHessLennardJones::operator()(
  const Real rawEnergy, const Real rawForce, Real a, Real /*rDistSquared*/,
  const Vector3D &rij, const GenericTopology *topo, int atom1, int atom2,
  const Real switchingValue, const Real switchingDeriv,
  const Matrix3By3 &switchingHess, ExclusionClass excl) const {
  Matrix3By3 vec_rij_ij(rij, rij);    // outer products of vectors r_ij
  Real a4 = a * a * a * a;
  Real a7 = a4 * a * a * a;
  Real a5 = a4 * a;
  Real a8 = a4 * a4;

  const LennardJonesParameters &params =
    topo->lennardJonesParameters(topo->atoms[atom1].type,
                                 topo->atoms[atom2].type);
  Real const_A = params.A;
  Real const_B = params.B;
  if (excl == EXCLUSION_MODIFIED) {
    const_A = params.A14;
    const_B = params.B14;
  }

  Real tm1 = 6 * const_B / a4 - 12 * const_A / a7;
  Real tm2 = -48 * const_B / a5 + 168 * const_A / a8;

  Matrix3By3 H(tm1, 0, 0, 0, tm1, 0, 0, 0, tm1);
  // now  is tm1 times identity matrix

  H += vec_rij_ij * tm2;
  H *= switchingValue;
  // now finished the second derivative of the raw energy.

  H += (switchingHess * rawEnergy);
  H -= vec_rij_ij * (2 * rawForce * switchingDeriv);

  return H;
}
