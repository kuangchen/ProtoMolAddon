#include <string>
#include <sstream>
#include <iostream>

#include <protomol/addon/snapshot/SnapshotTimeQueue.h>

using namespace ProtoMolAddon::Snapshot;
using std::istream;

void SnapshotTimeQueue::Reset(double t0, double dt, size_t n) {

  for (unsigned i=0; i<n; i++)
    push(t0 + dt * i);
}


bool SnapshotTimeQueue::IsDue(double t) {
  return (!empty() && front() < t);
}


namespace ProtoMolAddon {
  namespace Snapshot {

    istream& operator>> (istream &is, SnapshotTimeQueue &tq) {
      double t0, dt;
      size_t n;
      is >> t0 >> dt >> n;

      tq.Reset(t0, dt, n);

      return is;
    }
  }
}
