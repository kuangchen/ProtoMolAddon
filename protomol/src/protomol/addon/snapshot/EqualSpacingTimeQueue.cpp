#include <protomol/addon/snapshot/EqualSpacingTimeQueue.h>
#include <stdexcept>
#include <iostream>

using namespace ProtoMolAddon::Snapshot;
using std::invalid_argument;
using std::runtime_error;

EqualSpacingTimeQueue::EqualSpacingTimeQueue(double t0, double dt, size_t n): 
  t0(t0), dt(dt), n(n), current(0) {

  if (n<0)
    throw invalid_argument("n < 0");
  if (dt<0)
    throw invalid_argument("dt < 0");

}

double EqualSpacingTimeQueue::PopFront() { 
  return t0 + dt * current++;
}

size_t EqualSpacingTimeQueue::Size() const {
  return n;
}

bool EqualSpacingTimeQueue::IsDue(double now) {
  return ((current < n) && (now > t0 + current * dt));
}

namespace ProtoMolAddon {	
  namespace Snapshot {

    istream& operator>> (istream &is, EqualSpacingTimeQueue &spec) {
      is >> spec.t0 >> spec.dt >> spec.n;
      spec.current = 0;
  
      return is;
    }
  }
}
