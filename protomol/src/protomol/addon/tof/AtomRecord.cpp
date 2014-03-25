#include <protomol/addon/tof/AtomRecord.h>
#include <protomol/addon/Constants.h>
#include <protomol/type/Vector3D.h>
#include <iostream>

using std::ostream;

using namespace ProtoMolAddon::ToF;
using namespace ProtoMolAddon::Constant;

AtomRecord::AtomRecord(const string &name,
		       const Vector3DB &pos, 
		       const Vector3DB &vel) :
  name(name), pos(pos), vel(vel), pos_hit(), vel_hit(), time_hit(-1), status(0) {
}

namespace ProtoMolAddon {
  namespace ToF {

    ostream& operator<< (ostream &os, const AtomRecord &r) {
      os << r.name << "\t" 
	 << r.status << "\t"
	 << r.time_hit << "\t" 
	 << r.pos_hit * POSITION_CONV << "\t" 
	 << r.vel_hit * VELOCITY_CONV << std::endl;

      return os;
    }
  }
}

