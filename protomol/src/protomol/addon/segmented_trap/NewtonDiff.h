#ifndef _NEWTON_DIFF_H
#define _NEWTON_DIFF_H

#include <algorithm>
#include <vector>

namespace ProtoMolAddon {
  namespace SegmentedTrap {

    using namespace std;
    
    template <class Electrode>
    class NewtonDiff {
    private:
      const vector<Electrode> &elct;
      double epsilon;

    public:
      NewtonDiff(const vector<Electrode> &elct, double epsilon=1e-10);
      Vector3D Evaluate(const Vector3D &pos, double t) const;
    };

    template <class Electrode>    
    NewtonDiff<Electrode>::NewtonDiff(const vector<Electrode> &elct, double epsilon) 
      : elct(elct), epsilon(epsilon)
    {}
    
    template <class Electrode>    
    Vector3D NewtonDiff<Electrode>::Evaluate(const Vector3D &pos, double t) const {
      Vector3D g, a;
  
      for (int i=0; i<3; i++) {
	a[i] = epsilon;
	for (int j=0; j<elct.size(); j++)
	  g[i] += -(elct[j].GetPotential(pos+a, t) - elct[j].GetPotential(pos-a, t));
	a[i] = 0;
      }

      return g / (2*epsilon);
    }

  }
}

#endif
