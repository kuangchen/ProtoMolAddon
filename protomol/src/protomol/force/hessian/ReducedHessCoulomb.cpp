#include <protomol/force/hessian/ReducedHessCoulomb.h>
#include <protomol/topology/GenericTopology.h>

using namespace ProtoMol;
//____ ReducedHessCoulomb
Matrix3By3 ReducedHessCoulomb::operator()(const Real rawEnergy,
                                          const Real rawForce, Real a,
                                          Real /*rDistSquared*/,
                                          const Vector3D &rij,
                                          const GenericTopology *topo,
                                          int atom1, int atom2,
                                          const Real switchingValue,
                                          const Real switchingDeriv,
                                          const Matrix3By3 &switchingHess,
                                          ExclusionClass excl) const {
  Real scaledCharge_i = topo->atoms[atom1].scaledCharge;
  Real scaledCharge_j = topo->atoms[atom2].scaledCharge;
  Real coulombScalingFactor = 1;
  if (excl == EXCLUSION_MODIFIED)
    coulombScalingFactor = topo->coulombScalingFactor;

  Matrix3By3 vec_rij_ij(rij, rij);   // outer products of vectors r_ij

  Real na = sqrt(a);

  Matrix3By3 I(1, 0, 0, 0, 1, 0, 0, 0, 1);   // now I is identity matrix

  Real tm1 = -coulombScalingFactor * scaledCharge_i * scaledCharge_j / a / na;
  Matrix3By3 H((I - vec_rij_ij * (3 / a)) * (tm1 * switchingValue));
  // now done with the electrostatic (Coulomb) part.

  H += (switchingHess * rawEnergy);
  H -= vec_rij_ij * (2 * rawForce * switchingDeriv);

  return H;
}
