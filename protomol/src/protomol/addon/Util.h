#ifndef __UTIL_H
#define __UTIL_H

#include <protomol/type/Vector3DBlock.h>

using namespace ProtoMol;

namespace ProtoMolAddon {
  namespace Util {

    Vector3D BuildVector(double r, double theta, double phi);

    // Apply the rotation, that transforms (0, 0, 1) to v1, to v2;
    Vector3D Rotate(const Vector3D &v1, const Vector3D &v2);

  }
}


#endif
