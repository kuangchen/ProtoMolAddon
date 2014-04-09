#ifndef __ION_NEUTRAL_COLLISION_H
#define __ION_NEUTRAL_COLLISION_H

#include <protomol/type/Vector3D.h>
#include <vector>
#include <algorithm>
#include <unordered_map>
#include <random>
#include <stdexcept>

namespace ProtoMolAddon {
  namespace SympatheticCooling {

    using namespace ProtoMol;
    using namespace std;

    template <class CrossSection>
    class IonNeutralCollision {
    private:
      CrossSection cs;
      
    public:
      IonNeutralCollision(const CrossSection::initializer &init);
      Vector3D Rotate(double mu, const Vector3D &v) const;
      double GetReactionRateConstant(const ProtoMolIonProxy &ion) const;
    };

    template <class CrossSection>
    IonNeutralCollision<CrossSection>::IonNeutralCollision(const CrossSection::initializer &init) :
      cs(init) {

    }

    template <class CrossSection>
    Vector3D IonNeutralCollision<CrossSection>::Rotate(double mu, const Vector3D &v_rel) const {
      double energy = v_rel.normSquared() * 0.5 * mu ;
      pair<double, double> solid_angle = CrossSection.ResampleSolidAngle(energy);
      Vector3D v_rel_after = BuildVector(v_rel.norm(), solid_angle.first, solid_angle.second);
      return Rotate(v_rel, Vector3D &v_rel);
    }


  }
}

#endif
