#ifndef __DAMPING_H
#define __DAMPING_H

#include <protomol/type/Vector3D.h>
#include <string>

using namespace std;
using namespace ProtoMol;

namespace ProtoMolAddon {
  
  class Damping {
  private:
    double coeff;
    double start;
    double end;

  public:
    Damping();
    Damping(const string &def);
    void GetForce(const Vector3D &vel, double time, Vector3D& force);
  };
};


#endif
