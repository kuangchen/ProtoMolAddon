#include <protomol/addon/Damping.h>
#include <protomol/addon/LuaState.h>

using namespace ProtoMolAddon;
using namespace ProtoMolAddon::Lua;

Damping::Damping() :
  coeff(0),
  start(0),
  end(0)
{}

Damping::Damping(const string &def) :
  coeff(0),
  start(0),
  end(0)
{
  LuaState ls(def);
  
  coeff = ls.get<double>("damping.coeff");
  start = ls.get<double>("damping.start");
  end = ls.get<double>("damping.end");
}

void Damping::GetForce(const Vector3D &vel, double time, Vector3D& force) {
  double c = ((start<time) && (end>time)) ? coeff : 0;
  force = - vel * c;
}
