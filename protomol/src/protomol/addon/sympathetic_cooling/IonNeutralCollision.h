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

    class IonNeutralCollision {
    private:

      class CrossSection {
	
      };

      struct energy_hash {
	vector<double> energy;

	energy_hash(vector<double> &e) : energy(e) {
	  if (energy.empty())
	    throw runtime_error("Nonzero length expected for energy cross-section");

	  sort(energy.begin(), energy.end());
	}

	size_t operator() (double e) {
	  vector<double>::iterator iter = lower_bound(energy.begin(), 
						      energy.end(), 
						      e);
	  // Move iter back one step if the energy is too high
	  if (iter==energy.end())
	    iter--;

	  return hash<double>()(*iter);
	}
      };

      unordered_map<double, CrossSection > cross_section;
      
    public:
      IonNeutralCollision();
      Vector3D Rotate(const Vector3D &v) const;
    };

  }
}

#endif
