#include <string>

#include <protomol/addon/HarmonicTrap.h>
#include <protomol/addon/LuaConfigReader.h>

using namespace std;
using namespace ProtoMolAddon;
using namespace ProtoMolAddon::Lua;

HarmonicTrap::HarmonicTrap() :
  freq(3)
{}

HarmonicTrap::HarmonicTrap(const string &def) : 
  freq(3) 
{
  LuaConfigReader reader(def);
  
  freq[0] = reader.GetValue<double>("trap.freq.x");
  freq[1] = reader.GetValue<double>("trap.freq.y");
  freq[2] = reader.GetValue<double>("trap.freq.z");

}

void HarmonicTrap::GetForce(const Vector3D &pos, double mass, Vector3D& force) {
  for (int i=0; i<3; i++)
    force[i] = -mass * freq[i] * freq[i] * pos[i];
}
