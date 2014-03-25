#ifndef _ATOM_RECORD_H
#define _ATOM_RECORD_H

#include <protomol/type/Vector3D.h>
#include <iosfwd>

using namespace ProtoMol;

namespace ProtoMolAddon {
  
  namespace ToF {
    
    class AtomRecord {
    public:
      const std::string &name;
      const Vector3DB &pos;
      const Vector3DB &vel;

      Vector3D pos_hit;
      Vector3D vel_hit;

      double time_hit;
      int status;

    public:
      AtomRecord(const string &name, const Vector3DB &pos, const Vector3DB &vel);
      friend std::ostream& operator<< (std::ostream &os, const AtomRecord &r);
    };

  }
}

#endif
