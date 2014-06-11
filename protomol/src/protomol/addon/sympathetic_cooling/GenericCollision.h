#ifndef __GENERIC_COLLISION_H
#define __GENERIC_COLLISION_H

#include <protomol/type/Vector3D.h>
#include <protomol/addon/util/SIAtomProxy.h>

namespace ProtoMolAddon {
  namespace SympatheticCooling {

    using namespace ProtoMol;
    using namespace std;

    class GenericScattering {
    public:
      GenericScattering() {};
      virtual Vector3D Rotate(double mu, const Vector3D &v) const = 0;
      virtual double GetReactionRateConstant(const SIAtomProxy &proxy) const = 0;
    };

  }
}




    // template <class CrossSection>
    // class IonNeutralCollision {
    // private:
    //   CrossSection cs;
      
    // public:
    //   IonNeutralCollision(const CrossSection::initializer &init);
    //   Vector3D Rotate(double mu, const Vector3D &v) const;
    //   double GetReactionRateConstant(const ProtoMolIonProxy &ion) const;
    // };

    // template <class CrossSection>
    // IonNeutralCollision<CrossSection>::IonNeutralCollision(const CrossSection::initializer &init) :
    //   cs(init) {

    // }

    // template <class CrossSection>
    // Vector3D IonNeutralCollision<CrossSection>::Rotate(double mu, const Vector3D &v_rel) const {
    //   double energy = v_rel.normSquared() * 0.5 * mu ;

    //   pair<double, double> solid_angle = CrossSection.ResampleSolidAngle(energy);

    //   Vector3D v_rel_after = BuildVector(v_rel.norm(), solid_angle.first, solid_angle.second);

    //   return v_rel_after;
    // }


  }
}

#endif
