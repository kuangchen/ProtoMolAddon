#ifndef __ISOTROPIC_SCATTERING_H
#define __ISOTROPIC_SCATTERING_H

#include <protomol/addon/sympathetic_cooling/GenericScattering.h>
#include <protomol/addon/sympathetic_cooling/ThermalAtom.h>
#include <protomol/type/Vector3D.h>
#include <random>

namespace ProtoMolAddon {
  namespace SympatheticCooling {

    using namespace std;
    
    class IsotropicScattering : public GenericScattering {
    private:
      random_device rd;
      mutable uniform_real_distribution<double> z_dice;
      mutable uniform_real_distribution<double> phi_dice;

    public:
      typedef int initializer;

      IsotropicScattering(const initializer &init = initializer());

      Vector3D Rotate(double mu, const Vector3D &v) const;
      double GetReactionRateConstant(const ThermalAtom &atom, const SIAtomProxy &ion) const;
    };
    
  }
}



#endif
