#ifndef _THERMAL_ATOM_H
#define _THERMAL_ATOM_H

#include <protomol/type/Vector3D.h>
#include <string>
#include <random>

namespace ProtoMolAddon {
  namespace SympatheticCooling {
    
    using namespace ProtoMol;

    class ThermalAtom {
    public:
      typedef struct Spec {
	std::string name;
	double mass;
	double density;
	double temperature;
	double polarizability;

	Spec();
	Spec(const std::string &fname);
      } Spec;

    private:
      Spec spec;
      double C4;
      Vector3D position;
      Vector3D velocity;
      std::random_device rd;
      std::normal_distribution<double> dice;
      
    public:
      ThermalAtom();
      ThermalAtom(const Spec &spec);

      void Resample();
      
      inline Vector3D GetVelocity() const { return velocity; }
      inline Vector3D GetPosition() const { return position; }

      inline double GetMass() const  { return spec.mass; }
      inline double GetDensity() const { return spec.density; }
      inline double GetC4() const { return C4; }
    };

  }
}

#endif
