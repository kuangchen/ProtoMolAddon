#ifndef __THERMAL_ATOM_H
#define __THERMAL_ATOM_H

#include <protomol/type/Vector3D.h>
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
      double density;
      double temperature;
      random_device rd;
      normal_distribution<double> dice;
    
    public:
      ThermalAtom(double m, double density, double temperature);
      void Resample();

    };
  }
}

#endif
