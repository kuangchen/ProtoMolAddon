#ifndef _THERMAL_ATOM_H
#define _THERMAL_ATOM_H

#include <protomol/type/Vector3D.h>
#include <string>
#include <random>

using namespace ProtoMol;
using namespace std;

namespace ProtoMolAddon {
  namespace SympatheticCooling {

    using namespace std;

    class ThermalAtom {
    private:     
      double mass;
      Vector3D position;
      Vector3D velocity;
      double polarizability;
      double density;
      double temperature;
      double C4;
      random_device rd;
      normal_distribution<double> dice;
    
    public:
      typedef string initializer;

      ThermalAtom();
      ThermalAtom(const initializer &fname);
      void Resample();

      inline Vector3D GetVelocity() const { return velocity; }
      inline Vector3D GetPosition() const { return position; }
      inline double GetMass() const  { return mass; }
      inline double GetDensity() const { return density; }
      inline double GetC4() const { return C4; }
    };
  }
}

#endif
