#include <protomol/addon/tof/CEMDetector.h>
#include <protomol/addon/Constants.h>


using namespace ProtoMolAddon::ToF;
using namespace ProtoMolAddon::Constant;
using namespace ProtoMol;


CEMDetector::CEMDetector(const init_type &init):
  pos(init) {}

void CEMDetector::Process(AtomRecord *r, double now) const {
  if (r->status == 0 && r->pos[0] * POSITION_CONV > pos) {
    r->status = 1;
    r->pos_hit = r->pos;
    r->vel_hit = r->vel;
    r->time_hit = now;
  }
}

namespace ProtoMolAddon {
  namespace ToF {

    ostream& operator << (ostream &os, const CEMDetector &detector) {
      os << detector.pos << std::endl;

      return os;
    }

  }
}
