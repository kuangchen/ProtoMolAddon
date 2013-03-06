#include <protomol/force/hessian/ReducedHessCoulombDiElec.h>
#include <protomol/topology/GenericTopology.h>

using namespace ProtoMol;
//____ ReducedHessCoulombDielec
Matrix3By3 ReducedHessCoulombDielec::operator()(
  const Real rawEnergy, const Real rawForce, Real a, Real /*rDistSquared*/,
  const Vector3D &rij, const GenericTopology *topo, int atom1, int atom2,
  const Real switchingValue, const Real switchingDeriv,
  const Matrix3By3 &switchingHess, ExclusionClass excl, const Real D,
  const Real S,
  const Real epsi) const {
  //cout<<"Inside CoulombDiElec "<<endl;

  Real scaledCharge_i = topo->atoms[atom1].scaledCharge;
  Real scaledCharge_j = topo->atoms[atom2].scaledCharge;
  Real coulombScalingFactor = 1;
  if (excl == EXCLUSION_MODIFIED)
    coulombScalingFactor = topo->coulombScalingFactor;

  Real na = sqrt(a);

  Real expsr = exp(-S * na);
  Real epsilon_r =
    D - (D - epsi) * 0.5 * expsr * (S * S * na * na + 2 * S * na + 2);
  Real k_r = na * epsilon_r;
  Real z_r = D - (D - epsi) * 0.5 * expsr *
    (S * S * na * na + 2 * S * na + 2 - S * S * S * na * na * na);
  Real tm1 = (coulombScalingFactor * scaledCharge_i * scaledCharge_j);
  ///(k_r*k_r);

  Real Dervz_rVsr =
    (D - epsi) * 0.5 * expsr * (4 * S * S * S * a - S * S * S * S * a * na);

  Real P1r = -(z_r / na) / (k_r * k_r);

  Real P2r = (-1) * (k_r * Dervz_rVsr - 2 * z_r * z_r) / (k_r * k_r * k_r) /
    na - P1r / a;

  Matrix3By3 vec_rij_ij(rij, rij);   // outer products of vectors r_ij
  Matrix3By3 I(1, 0, 0, 0, 1, 0, 0, 0, 1);   // now I is identity matrix

  Matrix3By3 H((I * P1r + vec_rij_ij * (P2r / na)) * (tm1 * switchingValue));
  // now done with the electrostatic (Coulomb) part.

  H += (switchingHess * rawEnergy);
  H -= vec_rij_ij * (2 * rawForce * switchingDeriv);

  return H;
}
