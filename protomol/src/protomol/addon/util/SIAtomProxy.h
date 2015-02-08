#ifndef _SI_ATOM_PROXY
#define _SI_ATOM_PROXY

#include <string>
#include <protomol/base/PMConstants.h>
#include <protomol/addon/Constants.h>
#include <protomol/ProtoMolApp.h>
#include <protomol/type/Vector3DBlock.h>

namespace ProtoMolAddon {
  namespace Util {
    
    using namespace ProtoMol;
    using namespace ProtoMolAddon::Constant;

    class SIAtomProxy {
    private:
      unsigned int id;
      std::string *name;
      double *mass;
      double *charge;
      Vector3DB *pos;
      Vector3DB *vel;

    public:
      SIAtomProxy();
      SIAtomProxy(const SIAtomProxy &other);
      SIAtomProxy(ProtoMolApp *app, unsigned int i);
      SIAtomProxy(GenericTopology *topo, Vector3DBlock *positions, Vector3DBlock *velocities, unsigned int i);

      inline unsigned int GetID() const { return id; }
      inline Vector3D GetPosition() const { return *pos * ToSI::position; }
      inline Vector3D GetVelocity() const { return *vel * ToSI::velocity; }

      inline void SetPosition(const Vector3D &p) { (*pos) = p / ToSI::position; }
      inline void SetVelocity(const Vector3D &v) { (*vel) = v / ToSI::velocity; }
      
      inline double GetMass() const { return *mass * ToSI::mass; }
      inline void SetMass(double m) const { *mass = m / ToSI::mass; }

      inline double GetCharge() const { return *charge * ToSI::charge; }
      inline int GetIntegerCharge() const { return int(GetCharge() / SI::ELECTRON_CHARGE+0.5); }
      inline void SetIntegerCharge(int c) const { *charge = c * SI::ELECTRON_CHARGE / ToSI::charge; }
      inline void SetCharge(double c) { *charge = c / ToSI::charge; }

      inline void SetName(const std::string &n) { *name = n; }
      inline std::string GetName() const { return *name; }

      bool operator< (SIAtomProxy &other) { return this->id < other.id; }
      bool operator== (const std::string &name) { return *(this->name)==name; }

      inline friend bool operator== (const SIAtomProxy& ap, const std::string &name) { return ap.GetName() == name; }
    };
  }
}

#endif
