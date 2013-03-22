#include <string>

#include <protomol/addon/HarmonicTrap.h>
#include <protomol/addon/LuaState.h>

using namespace std;
using namespace ProtoMolAddon;
using namespace ProtoMolAddon::Lua;

HarmonicTrap::HarmonicTrap() :
  freq(3)
{}

HarmonicTrap::HarmonicTrap(const string &def) : 
  freq(3) 
{
  LuaState ls(def);
  
  freq[0] = ls.get<double>("freq.x");
  freq[1] = ls.get<double>("freq.y");
  freq[2] = ls.get<double>("freq.z");

}

void HarmonicTrap::GetForce(const Vector3D &pos, double mass, Vector3D& force) {

  for (int i=0; i<3; i++)
    force[i] = -mass * freq[i] * freq[i] * pos[i];

}
