#include <protomol/addon/Damping.h>
#include <protomol/addon/LuaConfigReader.h>

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
  LuaConfigReader reader(def);
  
  coeff = reader.GetValue<double>("damping.coeff");
  start = reader.GetValue<double>("damping.t_start");
  end = reader.GetValue<double>("damping.t_end");
}

void Damping::GetForce(const Vector3D& x, const Vector3D& v, double time, Vector3D& f) {
  double c = ((start<time) && (end>time)) ? coeff : 0;
  f = -(v+x*c)*2*c;
}
