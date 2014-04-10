#ifndef _PROTOMOL_ION_PROXY_H
#define _PROTOMOL_ION_PROXY_H

#include <protomol/type/Vector3D.h>
#include <protomol/addon/Constants.h>
#include <protomol/ProtoMolApp.h>

namespace ProtoMolAddon {
  namespace SympatheticCooling {

    using ProtoMol::ProtoMolApp;
    using ProtoMol::Vector3DB;
    using ProtoMol::Vector3D;

    
    using namespace ProtoMolAddon::Constant;

    class ProtoMolIonProxy {
    private:
      const double &mass;
      const double &charge;
      Vector3DB &position;
      Vector3DB &velocity;

    public:
      ProtoMolIonProxy(ProtoMolApp *app, unsigned int i);

      double GetMass() const { return mass * ToSI::mass; }
      double GetCharge() const { return charge * ToSI::charge; }
      const Vector3D GetPosition() const { return position * ToSI::position; } 
      const Vector3D GetVelocity() const { return velocity * ToSI::velocity; }
      void SetVelocity(const Vector3D &v) { velocity = v / ToSI::velocity; }
    };

  }
}

#endif
