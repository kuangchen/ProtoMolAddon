#include <protomol/addon/stray_field/StrayField.h>
#include <protomol/addon/Constants.h>
#include <iostream>

using namespace ProtoMolAddon::Constant;
using namespace ProtoMolAddon::StrayField;

StrayField::StrayField(const Vector3D &f) : field(f) {}

Vector3D StrayField::GetForce(double charge) {
  return field * charge;
}

namespace ProtoMolAddon {
  namespace StrayField {

    istream& operator>> (istream &is, StrayField &f) {
      is >> f.field[0] >> f.field[1] >> f.field[2];

      return is;
    }
  }
}


