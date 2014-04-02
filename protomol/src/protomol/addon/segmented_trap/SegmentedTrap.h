#ifndef _SEGMENTED_TRAP_H
#define _SEGMENTED_TRAP_H

#include <protomol/type/Vector3D.h>
#include <vector>

using namespace ProtoMol;

namespace ProtoMolAddon {
  namespace SegmentedTrap {

    using namespace std;

    template <class Electrode, class Diff>
    class SegmentedTrap {
    private:
      vector< Electrode > elct;
      Diff diff;
      void GetPotential(const Vector3D &pos, double t);

    public:
      SegmentedTrap();
      void GetForce(double charge, const Vector3D &pos, double t, Vector3D &force) const;
    };

    template <class Electrode, class Diff>
    SegmentedTrap<Electrode, Diff>::SegmentedTrap(): elct(12), diff(elct) {}

    template <class Electrode, class Diff>
    void SegmentedTrap<Electrode, Diff>::GetForce(double charge, const Vector3D &pos, double t, Vector3D& force) const {

      Vector3D gradient = diff.Evaluate(pos, t);
      force = gradient * charge;
    };

  }
}

#endif
