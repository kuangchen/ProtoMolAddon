#include <protomol/addon/BufferGasManager.h>
#include <iostream>
#include <ctime>

using namespace ProtoMol;
using namespace ProtoMolAddon;

double BufferGasManager::p_conv = 1e-10;
double BufferGasManager::v_conv = 1e-10 * SI::TIME_FS / TIMEFACTOR;

BufferGasManager::BufferGasManager() :
  buffer_gas(),
  collision_schedule(),
  next_collision(NULL),
  T(gsl_rng_default),
  r(gsl_rng_alloc(T))
{
  gsl_rng_set(r, time(NULL));
}

BufferGasManager::~BufferGasManager() {
  if (r)
    gsl_rng_free(r);
}

void BufferGasManager::InitializeBufferGas(LuaConfigReader &reader) {
  double mass = reader.GetValue<double>("buffer_gas.mass") * SI::AMU;
  double temp = reader.GetValue<double>("buffer_gas.temp");

  buffer_gas.Set(mass, temp);
}

void BufferGasManager::InitializeCollisionSchedule(LuaConfigReader &reader, ProtoMolApp *app) {
  double freq = reader.GetValue<double>("collision.freq");
  bool always_on = reader.GetValue<bool>("collision.always_on");

  int atom_count = app->topology->atoms.size();
  double start, end;

  always_on = true;
  if (always_on) {				  
    double step = app->integrator->getTimestep() / SI::TIME_FS;
    start = app->currentStep * step;
    end = app->lastStep * step;
  }
  
  collision_schedule.Generate(start, end, atom_count, r, freq);
  next_collision = collision_schedule.begin();
}

void BufferGasManager::Collide(ProtoMolApp *app) {
  //std::cout << next_collision->GetTime() << "\t" 
  //	    << next_collision->GetAtomId() << "\n";
  if (next_collision == collision_schedule.end() )
    return;

  int atomId = next_collision->GetAtomId();
  double m = app->topology->atoms[atomId].scaledMass * SI::AMU;
  Vector3D pi(app->positions[atomId] * p_conv);
  Vector3D vi(app->velocities[atomId] * v_conv);
  
  buffer_gas.Sample(r);
  buffer_gas.Collide(m, pi, vi, r);
  
  app->positions[atomId] = pi / p_conv;
  app->velocities[atomId] = vi / v_conv;
  
  next_collision++;
}

bool BufferGasManager::IsCollisionFinished() const {
  return next_collision == collision_schedule.end();
}

double BufferGasManager::GetNextCollisionTime() const {
  return IsCollisionFinished() ? 0 : next_collision->GetTime();

}
