#include <protomol/addon/BufferGas.h>
#include <protomol/base/PMConstants.h>

using namespace ProtoMolAddon;
using namespace ProtoMol;
using namespace ProtoMol::Constant;

BufferGas::BufferGas() :
  mass(0),
  temp(0),
  vn(0),
  pos(),
  vel()
{ 
}

void BufferGas::Set(double m, double t) {
  mass = m;
  temp = t;
  vn = sqrt(SI::BOLTZMANN * temp / mass);
}
   
void BufferGas::Sample(gsl_rng *r) {
  double mag = gsl_ran_chisq(r, 3);
  mag = sqrt(mag) * vn;
  gsl_ran_dir_3d(r, &vel[0], &vel[1], &vel[2]);
  vel *= mag;
}

void BufferGas::Collide(double mi, Vector3D &pi, Vector3D &vi, gsl_rng *r) {
  double b1 = mi/(mass+mi);
  double b2 = 1-b1;

  Vector3D v_com, v_rel, v_rel_after;
  v_com = vi * b1 + vel * b2;
  v_rel = vi - vel;

  double norm = v_rel.norm();
  gsl_ran_dir_3d(r, &v_rel_after[0], &v_rel_after[1], &v_rel_after[2]);
  v_rel_after *= norm;
  vi = v_com + v_rel_after * b2;
}

