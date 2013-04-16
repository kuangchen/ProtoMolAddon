#ifndef __BUFFERGAS_MANAGER_H
#define __BUFFERGAS_MANAGER_H

#include <protomol/ProtoMolApp.h>
#include <protomol/addon/BufferGas.h>
#include <protomol/addon/CollisionSchedule.h>
#include <protomol/addon/LuaConfigReader.h>
#include <protomol/base/PMConstants.h>

using namespace ProtoMol::Constant;
using namespace ProtoMolAddon::Lua;

namespace ProtoMolAddon {

  class BufferGasManager {

  private:
    BufferGas buffer_gas;
    CollisionSchedule collision_schedule;
    CollisionSchedule::iterator next_collision;
    const gsl_rng_type *T;
    gsl_rng *r;

  public:
    static double p_conv;
    static double v_conv;

    BufferGasManager();
    ~BufferGasManager();
    
    void InitializeBufferGas(LuaConfigReader &reader);
    void InitializeCollisionSchedule(LuaConfigReader &reader, ProtoMolApp *app);
    void Collide(ProtoMolApp *app);
    bool IsCollisionFinished() const;
    double GetNextCollisionTime() const ;
  };
  

};



#endif
