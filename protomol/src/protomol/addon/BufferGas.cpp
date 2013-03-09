#include <protomol/addon/BufferGas.h>
#include <protomol/base/PMConstants.h>
#include <cmath>

using namespace ProtoMol::Constant;
using namespace ProtoMolAddon::Collision;
using namespace ProtoMolAddon::Lua;

BufferGas::BufferGas() :
  _mass(0),
  _temperature(0),
  _freq(0),
  _vn(0),
  _schedule(0),
  _nextEvent(0),
  _T(NULL),
  _r(NULL)
{
}

BufferGas::~BufferGas() {
  if (_r)
    gsl_rng_free(_r);
}

        
BufferGas::BufferGas (LuaState& L) {
  gsl_rng_env_setup();  
  _T = gsl_rng_default;
  _r = gsl_rng_alloc(_T);
  
  _mass = L.get<double>("neutral.mass") * SI::AMU;
  _freq = L.get<double>("neutral.collision_frequency");
  _temperature = L.get<double>("neutral.temperature");
  _vn = sqrt(SI::BOLTZMANN * _temperature / _mass); 
}

void BufferGas::collide(double m, Vector3D& p, Vector3D& v) {
  Vector3D v_n, v_com, v_rel, v_rel_after;

  double mag = gsl_ran_chisq(_r, 3);
  mag = sqrt(mag) * _vn;

  gsl_ran_dir_3d(_r, &v_n[0], &v_n[1], &v_n[2]);
  v_n *= mag;
  
  double b1 = m / (_mass+m);
  double b2 = 1-b1;

  v_com = v * b1 + v_n * b2;
  v_rel = v - v_n;
  
  double norm = v_rel.norm();
  gsl_ran_dir_3d(_r, &v_rel_after[0], &v_rel_after[1], &v_rel_after[2]);
  v_rel_after *= norm;
  
  v = v_com + v_rel_after * b2;
}

bool Event::Compare(const Event p1, const Event p2) {
  return p1.time < p2.time;
}


void BufferGas::collide(ProtoMolApp *app){
  if (isCollisionFinished())
    return;

  int atomID = _schedule[_nextEvent].atomID;

  double p_conv = 1e-10;
  double v_conv = 1e-10 * SI::TIME_FS / TIMEFACTOR;

  Vector3D v_i(app->velocities[atomID]);
  Vector3D v_n, v_com, v_rel, v_rel_after;

  double mag = gsl_ran_chisq(_r, 3);
  mag = sqrt(mag) * _vn;

  gsl_ran_dir_3d(_r, &v_n[0], &v_n[1], &v_n[2]);
  v_n *= mag;

  double m = app->topology->atoms[atomID].scaledMass * SI::AMU;
  double b1 = m / (_mass+m);
  double b2 = 1-b1;

  v_com = v_i * b1 + v_n * b2;
  v_rel = v_i - v_n;
  
  double norm = v_rel.norm();
  gsl_ran_dir_3d(_r, &v_rel_after[0], &v_rel_after[1], &v_rel_after[2]);
  v_rel_after *= norm;
  
  app->velocities[atomID] = (v_com + v_rel_after * b2) / v_conv;
  _nextEvent++;
}

void BufferGas::scheduleCollision(double start, double end, int numAtom) {
  for (int n=0; n<numAtom; n++) 
    for (double t=start; t<end; t += gsl_ran_exponential(_r, 1.0/_freq))
      _schedule.push_back(Event(n, t));

  sort(_schedule.begin(), _schedule.end(), &Event::Compare);
}

vector<Event> BufferGas::createCollisionSchedule(double start, double end, int numAtom) {
  vector<Event> s(0);

  for (int n=0; n<numAtom; n++) 
    for (double t=start; t<end; t += gsl_ran_exponential(_r, 1.0/_freq))
	s.push_back(Event(n, t));

  sort(s.begin(), s.end(), &Event::Compare);
  
  return s;
}

bool BufferGas::isCollisionFinished() const {
  return _nextEvent==_schedule.size();
}

double BufferGas::getNextCollisionEventTime() const {
  return isCollisionFinished() ? 0 : _schedule[_nextEvent].time;
    
}
