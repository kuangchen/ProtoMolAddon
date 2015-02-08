#ifndef __CONST_SI_ATOM_PROXY_H
#define __CONST_SI_ATOM_PROXY_H

#include <string>
#include <protomol/base/PMConstants.h>
#include <protomol/addon/Constants.h>
#include <protomol/ProtoMolApp.h>
#include <protomol/type/Vector3DBlock.h>

namespace ProtoMolAddon {
  namespace Util {
    
    using namespace ProtoMol;
    using namespace ProtoMolAddon::Constant;

    class ConstSIAtomProxy {
    private:
      unsigned int id;
      const std::string *name;
      const double *mass;
      const double *charge;
      const Vector3DB *pos;
      const Vector3DB *vel;

    public:
      ConstSIAtomProxy();
      ConstSIAtomProxy(const ConstSIAtomProxy &other);
      ConstSIAtomProxy(const ProtoMolApp *app, unsigned int i);
      ConstSIAtomProxy(const GenericTopology *topo, const Vector3DBlock *positions, const Vector3DBlock *velocities, unsigned int i);

      inline unsigned int GetID() const { return id; }
      inline Vector3D GetPosition() const { return *pos * ToSI::position; }
      inline Vector3D GetVelocity() const { return *vel * ToSI::velocity; }
      inline double GetMass() const { return *mass * ToSI::mass; }

      inline double GetCharge() const { return *charge * ToSI::charge; }
      inline int GetIntegerCharge() const { return int(GetCharge() / SI::ELECTRON_CHARGE+0.5); }

      inline std::string GetName() const { return *name; }

      bool operator< (const ConstSIAtomProxy &other) const { return this->id < other.id; }
      bool operator== (const std::string &name) const { return *(this->name)==name; }
      inline friend bool operator== (const ConstSIAtomProxy& ap, const std::string &name) { return ap.GetName() == name; }
    };
  }
}

#endif
