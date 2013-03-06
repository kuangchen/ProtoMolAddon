#ifndef __BUFFER_GAS_H
#define __BUFFER_GAS_H

#include <protomol/output/LuaState.h>
#include <protomol/type/Real.h>
#include <protomol/type/Vector3D.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include <vector>
#include <string>

using std::string;
using LuaState::LuaState;

namespace ProtoMol {

  struct CollisionEvent {
    unsigned int atom;
    Real time;
  public:
    CollisionEvent(unsigned int atom=0, Real time=0) : atom(atom), time(time) {}
  };

  class BufferGas {
  private:
    Real mass;
    Real temperature;
    Real freq;
    Real vn;

    const gsl_rng_type *T;
    gsl_rng *r;

  public:
    BufferGas();
    BufferGas(LuaState::LuaState& L);
    ~BufferGas();

    void collide(Real mass, Vector3D& pos, Vector3D& vel);
    vector<CollisionEvent> createCollisionSchedule(Real start, Real end, int atomCount);
  };
}

#endif
