#ifndef _STRAY_FIELD_H
#define _STRAY_FIELD_H

#include <iosfwd>
#include <protomol/type/Vector3D.h>

using std::istream;
using ProtoMol::Vector3D;

namespace ProtoMolAddon {
  namespace StrayField {
    
    class StrayField {
    private:
      Vector3D field;

    public:
      StrayField(const Vector3D &f = Vector3D());
      Vector3D GetForce(double charge);

      friend istream& operator>> (istream &is, StrayField &f);
    };
  }
}

#endif

