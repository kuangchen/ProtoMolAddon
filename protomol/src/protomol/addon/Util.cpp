#include <protomol/addon/Util.h>
#include <cmath>

namespace ProtoMolAddon {
  namespace Util {

    Vector3D BuildVector(double r, double theta, double phi) {
      return Vector3D(sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta)) * r; 
    }

    Vector3D Rotate(const Vector3D &v1, const Vector3D &v2) {
      Vector3D temp = v1 / v1.norm();

      double ct = temp[2];
      double st = sqrt(temp[0] * temp[0] + temp[1] * temp[1]);

      double sp = 1; 
      double cp = 0;
      if (fabs(st) > 1e-8) {
	sp = temp[1] / st;
	cp = temp[0] / st;
      }

      return Vector3D( cp*ct*v2[0] - sp*v2[1] + cp*st*v2[2],
		       sp*ct*v2[0] + cp*v2[1] + sp*st*v2[2],
		       -  st*v2[0] +        0 +    ct*v2[2]);
    }
  }
}
