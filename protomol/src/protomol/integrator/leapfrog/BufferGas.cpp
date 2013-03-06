#include <protomol/integrator/leapfrog/BufferGas.h>
#include <protomol/base/PMConstants.h>
#include <cmath>

using namespace ProtoMol;

BufferGas::BufferGas() :
  mass(0),
  temperature(0),
  freq(0),
  vn(0),
  T(NULL),
  r(NULL)
{
}

BufferGas::~BufferGas() {
  if (r)
    gsl_rng_free(r);
}

        
BufferGas::BufferGas (LuaState::LuaState& L) {
  gsl_rng_env_setup();  
  T = gsl_rng_default;
  r = gsl_rng_alloc(T);
  
  mass = L.get<double>("neutral.mass") * ProtoMol::Constant::SI::AMU;
  freq = L.get<double>("neutral.collision_frequency");
  temperature = L.get<double>("neutral.temperature");
  vn = sqrt(ProtoMol::Constant::SI::BOLTZMANN * temperature / mass); 
}

void ProtoMol::BufferGas::collide(Real m, Vector3D& p, Vector3D& v) {
  Vector3D v_n, v_com, v_rel, v_rel_after;

  Real mag = gsl_ran_chisq(r, 3);
  mag = sqrt(mag) * vn;

  gsl_ran_dir_3d(r, &v_n[0], &v_n[1], &v_n[2]);
  v_n *= mag;
  
  Real b1 = m / (mass+m);
  Real b2 = 1-b1;

  v_com = v * b1 + v_n * b2;
  v_rel = v - v_n;
  
  double norm = v_rel.norm();
  gsl_ran_dir_3d(r, &v_rel_after[0], &v_rel_after[1], &v_rel_after[2]);
  v_rel_after *= norm;
  
  v = v_com + v_rel_after * b2;
}

bool compareFunction(const CollisionEvent p1, const CollisionEvent p2) {
  return p1.time < p2.time;
}

vector<CollisionEvent> BufferGas::createCollisionSchedule(Real start, Real end, int atomCount) {
  vector<CollisionEvent> s;

  for (int n=0; n<atomCount; n++) 
    for (double t = start; t<end; t += gsl_ran_exponential(r, 1.0/freq))
	s.push_back(CollisionEvent(n, t));

  sort(s.begin(), s.end(), compareFunction);
  
  return s;
}
