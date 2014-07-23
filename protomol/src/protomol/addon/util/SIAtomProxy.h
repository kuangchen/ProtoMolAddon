#ifndef _SI_ATOM_PROXY
#define _SI_ATOM_PROXY

#include <string>
#include <protomol/addon/Constants.h>
#include <protomol/ProtoMolApp.h>
#include <protomol/type/Vector3DBlock.h>

namespace ProtoMolAddon {
  namespace Util {
    
    using namespace ProtoMol;
    using namespace ProtoMolAddon::Constant;

    class ConstSIAtomProxy {
    private:
      const string *name;
      const double *mass;
      const double *charge;
      const Vector3DB *pos;
      const Vector3DB *vel;

    public:
      ConstSIAtomProxy(const ProtoMolApp *app, unsigned int i);

      ConstSIAtomProxy(const GenericTopology *topo, 
		       const Vector3DBlock *positions, 
		       const Vector3DBlock *velocities, 
		       const unsigned int i);

      inline Vector3D GetPosition() const { return *pos * ToSI::position; }
      inline Vector3D GetVelocity() const { return *vel * ToSI::velocity; }

      inline double GetMass() const { return *mass * ToSI::mass; }
      inline double GetCharge() const { return *charge * ToSI::charge; }
      inline string GetName() const { return *name; }
    };


    class SIAtomProxy {
    private:
      unsigned int id;
      const string *name;
      const double *mass;
      double *charge;
      Vector3DB *pos;
      Vector3DB *vel;

    public:
      SIAtomProxy();
      SIAtomProxy(ProtoMolApp *app, unsigned int i);
      SIAtomProxy(GenericTopology *topo, Vector3DBlock *positions, Vector3DBlock *velocities, unsigned int i);

      inline unsigned int GetID() const { return id; }
      inline Vector3D GetPosition() const { return *pos * ToSI::position; }
      inline Vector3D GetVelocity() const { return *vel * ToSI::velocity; }

      inline void SetPosition(const Vector3D &p) { (*pos) = p / ToSI::position; }
      inline void SetVelocity(const Vector3D &v) { (*vel) = v / ToSI::velocity; }
      
      inline double GetMass() const { return *mass * ToSI::mass; }
      inline double GetCharge() const { return *charge * ToSI::charge; }
      
      inline void SetCharge(double c) { *charge = c / ToSI::charge; }

      inline string GetName() const { return *name; }
      bool operator< (SIAtomProxy &other) { return this->id < other.id; }
    };
  }
}

#endif
