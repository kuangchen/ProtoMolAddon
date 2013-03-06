#include <protomol/force/hessian/ReducedHessAngle.h>
#include <protomol/base/Report.h>

using namespace std;
using namespace ProtoMol::Report;
using namespace ProtoMol;

//____constructors
ReducedHessAngle::ReducedHessAngle() {
  // create ReducedHessAngle structure without init it to zero
  // user should take care of it by calling clear() to init it to zero.
}

ReducedHessAngle::ReducedHessAngle(const Vector3D &atom_i,  //position of atom i
                                   const Vector3D &atom_j,  //position of atom j
                                   const Vector3D &atom_k,  //position of atom k
                                   const Real k_t,   //angluar spring constant
                                   const Real theta0,   //rest angle
                                   bool computeReduced) {
  evaluate(atom_i, atom_j, atom_k, k_t, theta0, computeReduced);
  // call the evaluate function
}


void ReducedHessAngle::evaluate(const Vector3D &atom_i,   //position of atom i
                                const Vector3D &atom_j,   //position of atom j
                                const Vector3D &atom_k,   //position of atom k
                                const Real k_t,   //angluar spring constant
                                const Real theta0,   //rest angle
                                bool computeReduced) {
  Vector3D rij(atom_j - atom_i);
  Vector3D rkj(atom_j - atom_k);
  Vector3D rki(atom_i - atom_k);
  // the relative position vectors

  Real theta = atan2((rij.cross(rkj)).norm(), rij.dot(rkj));
  Real sinTheta = sin(theta);
  Real cosTheta = cos(theta);

  Real a = rij.normSquared();   // \alpha: square of the norm of r_ij
  Real b = rkj.normSquared();   // \beta: square of the norm of r_kj
  Real c = rki.normSquared();   // \gamma: square of the norm of r_ki
  Real na = sqrt(a);    // norm of r_ij
  Real nb = sqrt(b);    // norm of r_kj
  Real nanb = na * nb;   // temp constant

  Real E_c = -2.0 * k_t * (theta - theta0) / sinTheta;
  // derivative of E_theta w.r.t. C
  Real E_cc = 2.0 * k_t * (sinTheta - (theta - theta0) * cosTheta);
  E_cc /= sinTheta * sinTheta * sinTheta;
  // derivative of E_c w.r.t. C

  Real f = (a - b + c) / (4 * nanb * a);
  Real g = (-a + b + c) / (4 * nanb * b);
  Real h = -1. / (2 * nanb);
  //constants computed from \alpha, \beta, \gamma

  // derivative of f w.r.t. \alpha
  Real f_a = (-a + 3 * b - 3 * c) / (8 * a * a * nanb);
  // derivative of f w.r.t. \beta
  Real f_b = (-a - b - c) / (8 * a * b * nanb);
  Real f_c = 1. / (4 * a * nanb);    // derivative of f w.r.t. \gamma
  //derivative of g w.r.t. \beta
  Real g_b = (3 * a - b - 3 * c) / (8 * nanb * b * b);
  Real g_c = 1. / (4 * b * nanb);   //derivative of g w.r.t. \gamma

  Real t1 = 2 * (f + h) * E_c, t2 = -2 * h * E_c;
  Real t3 = 2 * (g + h) * E_c;

  H[0][0](t1, 0, 0, 0, t1, 0, 0, 0, t1);
  H[0][2](t2, 0, 0, 0, t2, 0, 0, 0, t2);
  H[2][2](t3, 0, 0, 0, t3, 0, 0, 0, t3);

  // first init the subblocks with const mult identity matrices.

  Matrix3By3 rij_ij(rij, rij), rij_ki(rij, rki), rij_kj(rij, rkj);
  Matrix3By3 rjk_jk(rkj, rkj), rki_kj(rki, rkj), rki_ki(rki, rki);
  // use the 2-coordinates constructor to get the outer-produce matrix

  Real FourEc = 4 * E_c;
  Real FourEcFc = FourEc * f_c, FourEcGc = FourEc * g_c;
  Real FourEcFa = FourEc * f_a, FourEcFb = FourEc * f_b;
  Real FourEcGb = FourEc * g_b;

  H[0][0] += rij_ij * FourEcFa - (rij_ki + rij_ki.transposed()) * FourEcFc;
  H[0][2] += rij_kj * FourEcFb + rij_ki * FourEcFc - rki_kj * FourEcGc;
  H[2][2] += rjk_jk * FourEcGb + (rki_kj + rki_kj.transposed()) * FourEcGc;

  Real FourEcc = 4 * E_cc;
  ;
  Real FourEccFF = FourEcc * f * f, FourEccFG = FourEcc * f * g;
  Real FourEccFH = FourEcc * f * h, FourEccGG = FourEcc * g * g;
  Real FourEccGH = FourEcc * g * h, FourEccHH = FourEcc * h * h;

  H[0][0] += rij_ij * FourEccFF -
             (rij_ki + rij_ki.transposed()) * FourEccFH + rki_ki * FourEccHH;
  H[0][2] += rij_kj * FourEccFG + rij_ki * FourEccFH - rki_kj * FourEccGH -
             rki_ki * FourEccHH;
  H[2][2] += rjk_jk * FourEccGG +
             (rki_kj + rki_kj.transposed()) * FourEccGH + rki_ki * FourEccHH;

  if (!computeReduced) {
    Real t4 = -2 * f * E_c, t5 = 2 * (g + f) * E_c;
    Real t6 = -2 * g * E_c;
    H[0][1](t4, 0, 0, 0, t4, 0, 0, 0, t4);
    H[1][1](t5, 0, 0, 0, t5, 0, 0, 0, t5);
    H[1][2](t6, 0, 0, 0, t6, 0, 0, 0, t6);

    H[0][1] += rij_ki.transposed() * FourEcFc + rki_kj * FourEcGc - rij_ij *
               FourEcFa - rij_kj * FourEcFb;
    H[1][1] += rij_ij * FourEcFa +
               (rij_kj + rij_kj.transposed()) * FourEcFb + rjk_jk * FourEcGb;
    H[1][2] -= rjk_jk * FourEcGb + rij_kj * FourEcFb + rij_ki * FourEcFc +
               rki_kj.transposed() * FourEcGc;

    H[0][1] += rki_kj * FourEccGH - rij_kj * FourEccFG +
               rij_ki.transposed() * FourEccFH - rij_ij * FourEccFF;
    H[1][1] += rij_ij * FourEccFF +
               (rij_kj +
                rij_kj.transposed()) * FourEccFG + rjk_jk * FourEccGG;
    H[1][2] -= rij_kj * FourEccFG + rij_ki * FourEccFH +
               rki_kj.transposed() * FourEccGH + rjk_jk * FourEccGG;

    H[1][0].transpose(H[0][1]);
    H[2][0].transpose(H[0][2]);
    H[2][1].transpose(H[1][2]);
  } else {
    H[0][1](0, 0, 0, 0, 0, 0, 0, 0, 0);
    H[1][1](0, 0, 0, 0, 0, 0, 0, 0, 0);
    H[1][2](0, 0, 0, 0, 0, 0, 0, 0, 0);
    H[1][0](0, 0, 0, 0, 0, 0, 0, 0, 0);
    H[2][0].transpose(H[0][2]);
    H[2][1](0, 0, 0, 0, 0, 0, 0, 0, 0);
  }
}

Matrix3By3 ReducedHessAngle::operator()(int i, int j) const {
  return H[i][j];
}

void ReducedHessAngle::operator()(int i, int j, Matrix3By3 x) {
  H[i][j] = x;
}

void ReducedHessAngle::accumulateTo(int i, int j, Matrix3By3 x) {
  H[i][j] += x;
}

void ReducedHessAngle::accumulateNegTo(int i, int j, Matrix3By3 x) {
  H[i][j] -= x;
}

//____ print the ReducedHessAngle matrix to tty
ostream &ProtoMol::operator<<(ostream &os, const ReducedHessAngle &tm) {
  for (unsigned int i = 0; i < 3; i++)
    for (unsigned int j = 0; j < 3; j++)
      os << "H[" << i << "][" << j << "]" << ": " << endl << tm.H[i][j] <<
      endl;

  return os;
}

ReducedHessAngle ReducedHessAngle::transposed() {
  ReducedHessAngle res;

  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      res.H[i][j] = H[j][i].transposed();

  return res;
}

//____makes *this identity matrix
void ReducedHessAngle::identity() {
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      if (i != j)
        H[i][j](0, 0, 0, 0, 0, 0, 0, 0, 0);
      else
        H[i][j](1, 0, 0, 0, 1, 0, 0, 0, 1);

}

void ReducedHessAngle::convertFromJacobian(Real **jac, int n) {
  if (n == 9)
    for (int i = 0; i < 9; i++)
      for (int j = 0; j < 9; j++)
        (H[i / 3][j / 3])(i % 3, j % 3, jac[i][j]);

  else {   // assume n=6 holding 0,0; 0,2; 2,0; 2,2
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++)
        H[0][0](i % 3, j % 3, jac[i][j]);

    for (int i = 0; i < 3; i++)
      for (int j = 3; j < 6; j++)
        H[0][2](i % 3, j % 3, jac[i][j]);

    for (int i = 3; i < 6; i++)
      for (int j = 0; j < 3; j++)
        H[2][0](i % 3, j % 3, jac[i][j]);

    for (int i = 3; i < 6; i++)
      for (int j = 3; j < 6; j++)
        H[2][2](i % 3, j % 3, jac[i][j]);

  }
}

//____clear the hessian matrix
void ReducedHessAngle::clear() {
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      H[i][j](0, 0, 0, 0, 0, 0, 0, 0, 0);

}

//____ operator *
ReducedHessAngle ReducedHessAngle::operator *(const ReducedHessAngle &tm) {
  ReducedHessAngle res;
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) {
      (res.H[i][j])(0, 0, 0, 0, 0, 0, 0, 0, 0);
      for (int k = 0; k < 3; k++)
        res.H[i][j] += H[i][k] * tm.H[k][j];
    }

  return res;
}

//____ operator +
ReducedHessAngle ReducedHessAngle::operator+(const ReducedHessAngle &tm) {
  ReducedHessAngle res;
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      res.H[i][j] = H[i][j] + tm.H[i][j];

  return res;
}

//____ operator -
ReducedHessAngle ReducedHessAngle::operator-(const ReducedHessAngle &tm) {
  ReducedHessAngle res;
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      res.H[i][j] = H[i][j] - tm.H[i][j];

  return res;
}

//____ operator +=
ReducedHessAngle &ReducedHessAngle::operator+=(const ReducedHessAngle &tm) {
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      H[i][j] += tm.H[i][j];

  return *this;
}

//____ operator -=
ReducedHessAngle &ReducedHessAngle::operator-=(const ReducedHessAngle &tm) {
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      H[i][j] -= tm.H[i][j];

  return *this;
}

//____ operator *
ReducedHessAngle ReducedHessAngle::operator *(const Real tm) {
  ReducedHessAngle res;
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      res.H[i][j] = H[i][j] * tm;

  return res;
}

//____ operator /
ReducedHessAngle ReducedHessAngle::operator/(const Real tm) {
  ReducedHessAngle res;
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      res.H[i][j] = H[i][j] / tm;

  return res;
}

ReducedHessAngle &ReducedHessAngle::operator*=(const Real tm) {
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      H[i][j] *= tm;

  return *this;
}

ReducedHessAngle &ReducedHessAngle::operator/=(const Real tm) {
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      H[i][j] /= tm;

  return *this;
}

