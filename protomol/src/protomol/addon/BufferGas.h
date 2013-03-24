#ifndef __BUFFER_GAS_H_
#define __BUFFER_GAS_H_

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <protomol/type/Vector3D.h>

using namespace ProtoMol;
using namespace std;

namespace ProtoMolAddon {

  class BufferGas {
    double mass;
    double temp;
    double vn;
    Vector3D pos;
    Vector3D vel;

  public:
    BufferGas();
    void Set(double mass, double temp);
    void Sample(gsl_rng *r);
    void Collide(double mi, Vector3D &pi, Vector3D &vi, gsl_rng *r);
  };

  
}


#endif
